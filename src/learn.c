/*    ceif - categorized extended isolation forest
 *
 *    Copyright (C) 2019 Timo Savinen
 *    This file is part of ceif.
 * 
 *    ceif is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    ceif is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with ceif; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *    F607480034
 *    HJ9004-2
 *
 */
#include "ceif.h"
#include "cmap.h"
#include <math.h>
#include <regex.h>
#include <time.h>

static double fast_n_cache[FAST_N_SAMPLES];
static double fast_c_cache[FAST_C_SAMPLES];

static char input_line[INPUT_LEN_MAX];
static time_t now;
static double centroid_tresshold = CENTROID_TRESSHOLD;

/* hash function for hash table
 * calculates hash for string s
 */
static 
size_t hash(char *s)
{
    register unsigned long h = 5381;
    int c;

    while ((c = *s++) != 0)
    {
        h = ((h << 5) + h) + c;
    }

    return (size_t) (h % HASH_MAX);
}

/* search forest hash 
 * return -1 if not found, else index to forest table
 */
int search_forest_hash(char *category_string)
{
    int i;
    size_t h = hash(category_string);

    for(i = 0;i < fhash[h].idx_count;i++)
    {
        if (strcmp(category_string,forest[fhash[h].idx[i]].category) == 0) return (int) fhash[h].idx[i];
    }

    return -1;
}

/*
 * Add new entry to foretst hash 
 * idx is index to forest table
 */
void add_forest_hash(int idx, char *category_string)
{
    size_t h = hash(category_string);

    if(fhash[h].idx_count >= fhash[h].idx_cap)
    {
        fhash[h].idx_cap += 10;
        fhash[h].idx = xrealloc(fhash[h].idx,fhash[h].idx_cap * sizeof(size_t));
    }

    fhash[h].idx[fhash[h].idx_count] = idx;

    fhash[h].idx_count++;
}

/* generate a random integer from range min...max
 */ 
inline int ri(int min, int max)
{
    return ((rand() % (max - min + 1)) + min);
}

/* generate a random double from range min...max
 */ 
static 
inline double rd(double min, double max)
{
    return ((double) rand() * (max - min) / (double) (RAND_MAX)) + min;
}
    
/* Select forest using sample data and category dimensions. 
   If category dimensions are not used select first forest

   Allocates new forest if category dim is not found from existing forests

   Making string compare here, for large number of forests hash based search would be faster

   returns index to forest table
 */
static
int select_forest(int value_count,char **values)
{
    int i;
    char *category_string;


    category_string = make_category_string(value_count,values);

    i = search_forest_hash(category_string);
   
    if(i >= 0) {
        if(forest[i].last_updated != now) forest[i].last_updated = now;
        return i;
    }

    if(forest_count >= forest_cap) {
       if(forest_cap == 0) forest_cap = 64;
       forest_cap *= 2;
       forest = xrealloc(forest,forest_cap * sizeof(struct forest));
    }

   forest[forest_count].category = xstrdup(category_string);
   forest[forest_count].filter = 0;
   forest[forest_count].X_count = 0;
   forest[forest_count].X_current = 0;
   forest[forest_count].X_cap = 0;   
   forest[forest_count].X_summary = -1;   
   forest[forest_count].X = NULL;   
   forest[forest_count].t = NULL;
   forest[forest_count].min = NULL;
   forest[forest_count].max = NULL;
   forest[forest_count].scale_range_idx = -1;
   forest[forest_count].summary = NULL;
   forest[forest_count].c = 0;
   forest[forest_count].heigth_limit = 0;
   forest[forest_count].analyzed = 0;
   forest[forest_count].last_updated = now;
   forest[forest_count].avg = xmalloc(dimensions * sizeof(double));
   forest[forest_count].dim_density = xmalloc(dimensions * sizeof(double));
   forest[forest_count].total_rows = 0;
   forest[forest_count].analyzed_rows = 0;
   forest[forest_count].high_analyzed_rows = 0;
   forest[forest_count].extra_rows = 0;
   forest[forest_count].percentage_score = 0.0;
   forest[forest_count].min_score = 1.0;
   forest[forest_count].test_average_score = 0.0;

   add_forest_hash(forest_count,category_string);

   forest_count++;

   return forest_count - 1;
}

   
/* copy a vector
 */

void v_copy(double *t,double *s)
{
    memcpy(t,s,dimensions * sizeof(double));
}

/* duplicate a vector
 */
static
double *v_dup(double *v)
{
    double *new = xmalloc(dimensions * sizeof(double));

    v_copy(new,v);
    
    return new;
}

/* Compare vectors
 * return 0 if vectors are equal
 */
static
int v_cmp(double *t,double *s)
{
    return memcmp(t,s,dimensions * sizeof(double));
}


/* parse single numeric sample attribute. data expression is called to evaluate possible expression
   return 0 in case number cannot be parsed
   */
double parse_dim_attribute_expr(int data_idx, int value_count, char **values)
{
    double d;

    d = atof(evaluate_data_expression(data_idx,value_count,values));

    return isnan(d) ? 0.0 : d;
}

/* parse single numeric sample attribute. */
double parse_dim_attribute(char *value)
{
    double d;

    d = atof(value);

    return isnan(d) ? 0.0 : d;
}

/* parse single text attribute 
 * a hash is calculated based on string.
 * NOTE! collisions may occur
 *
 * Probably this should be used only for simple classifications like "yes/no" or "Male/Female/Unknown"
 * or text codes like A-Z or A1,A2,A3,...
 * 
 * Note also that "Yes" and "yes" will produce same value (upper and lower case are mapped to same value)
 *
 * Algorith is based on nothing, it just tries to produce different value for different strings, but probably not 100% sure
 * Also value for similar texts should be close, e.g.
 * values for 'attribute' and 'attricute' are relative close
 *
 * Values start app. 700 
 */
double parse_dim_hash_attribute(char *value)
{
    unsigned char *c,count=1,prev=85,pprev=60;
    unsigned long long int weight = 346;

    c = (unsigned char *) value;

    while(*c)
    {
        weight += cmap[*c].map_value * ((prev ^ count) + 1) * ((pprev ^ count) + 1);
        pprev = prev;
        prev = cmap[*c].map_value;
        if(count == 255) count = 0;
        count++;
        c++;
    }

    return (double) weight / 345.6789;
}

/* parse ascii value table to dimension table of type double
 */
void parse_values(double *dim,char **values, int value_count, int saved)
{
    int i;

    for(i = 0;i < dimensions;i++)
    {
        // silently ignore missing input dimension values
        if((saved ? i : dim_idx[i]) < value_count)
        {
            if(saved)
            {
                dim[i] = parse_dim_attribute(values[i]);
            } else
            {
                if(check_idx(dim_idx[i],text_idx_count,text_idx))
                {
                    dim[i] = parse_dim_hash_attribute(values[dim_idx[i]]);
                } else
                {
                    dim[i] = parse_dim_attribute_expr(dim_idx[i],value_count,values);
                }
            }
        } else
        {
            dim[i] = 0.0;   // use default value for missing dimension values
        }
    }
}

/* check if dimension values are allready in sample table X
 */
static
int duplicate_sample(struct forest *f,double *new_dim)
{
    int i;

    for(i = 0;i < f->X_count;i++)
    {
        if(v_cmp(new_dim,f->X[i].dimension) == 0) return 1;
    }
    return 0;
}

/* return scaled or non scaled sample dimension. Select dim using auto_weigth
 */
inline double *
sample_dimension(struct sample *s)
{
    return auto_weigth ? s->scaled_dimension : s->dimension;
}


/* add one sample to X in selected forest.
  if X is full (== samples_total) implement  reservoir sampling 
 * and check min/max values
 *
 * saved == 1 means that values are from saved forest file and == 0 means that values are from user given data file.
 * indices are different in those two cases
 */
void add_to_X(struct forest *f,double *new, int value_count,int saved)
{
    int sample_idx;

    DEBUG("Adding %sdimension to forest %s: ",saved ? "saved " : "",f->category);
    DEBUG_ARRAY(value_count,new);

    // check if this is a duplicate sample. Unique_samples is a value between 0..100. 
    // if value is .e.g 10 then 10 percent of samples are checked for uniqueness
    // check is not done for saved values, only new values are checked

    if(!saved && unique_samples > 0 && ri(0,100) <= unique_samples && duplicate_sample(f,new))
    {
        DEBUG(" Not added to sample table due to the duplicate check\n");
        return;
    }

    DEBUG(", Not a duplicate,");

    if(f->X_count >= f->X_cap) {
        if(f->X_cap == 0) f->X_cap = 32;
        f->X_cap *= 2;
        f->X = xrealloc(f->X,f->X_cap * sizeof(struct sample));
        DEBUG(" Reallocating sample table size to %u, ",f->X_cap);
    }

    if(f->X_count < samples_total)    // check the samples table size
    {
        DEBUG(" Adding as a new item to sample table");
        if(f->X_count == 0)
        {
            sample_idx = f->X_count;
            f->X[sample_idx].dimension = xmalloc(dimensions * sizeof(double));
            f->X[sample_idx].scaled_dimension = NULL;                            
        } else
        {
            sample_idx = ri(0,f->X_count - 1);             // ceif does not work well with sorted samples, make sure that  samples are shuffled
            f->X[f->X_count].dimension = v_dup(f->X[sample_idx].dimension);
            f->X[f->X_count].scaled_dimension = NULL;
        }
        f->X[f->X_count].cluster_center_idx = -1;
        f->X_count++;
    } else
    {
        DEBUG(" Replacing an existing item in sample table");
        sample_idx = ri(0,f->X_count + f->extra_rows -1); // minus 1 in order to get first sample to be saved

        if(!saved) f->extra_rows++;                  // Number of extra rows for this forest read from train file

        if(sample_idx >= f->X_count) return;         // check if old sample should be replaced with this or not
    }

    v_copy(f->X[sample_idx].dimension,new);

    DEBUG("\n");
}

/* Aggregate new values for a certain sample item in a forest.
 * if forest -> X_summary == -1, this is the first data to be aggregated
 * if forest -> X_summary > -1, this is the index to X where data is aggregated
 * in case we have max number of samples, a new data is allways added, no reservoir sampling is implemented
 */
static void
add_aggregate(struct forest *f,char **values, int value_count)
{
    int i,sample_idx;
    double *s;
    static double new[DIM_MAX];

    parse_values(new,values,value_count,0); 
    DEBUG("Adding dimension to be aggregated to forest %s: ",f->category);
    DEBUG_ARRAY(value_count,new);

    if(f->X_count >= f->X_cap) {
        if(f->X_cap == 0) f->X_cap = 32;
        f->X_cap *= 2;
        f->X = xrealloc(f->X,f->X_cap * sizeof(struct sample));
        DEBUG(" Reallocating sample table size to %u,",f->X_cap);
    }

    if(f->X_summary > -1)
    {
        sample_idx = f->X_summary;
    } else
    {
        if(f->X_count < samples_total)    // check the samples table size
        {
            DEBUG(" Adding as a new item to sample table,");
            sample_idx = f->X_count;
            f->X[sample_idx].dimension = xmalloc(dimensions * sizeof(double));
            f->X[sample_idx].scaled_dimension = NULL;
            f->X_count++;
        } else                                          // max number of samples in X
        {
            DEBUG(" Replacing an existing item in sample table,");
            sample_idx = ri(0,f->X_count - 1);          // replace a random sample with this new one
        }

        f->X_summary = sample_idx;

        for(i = 0;i < dimensions;i++) f->X[sample_idx].dimension[i] = 0.0; // init summary
    }
    
    s = f->X[sample_idx].dimension;

    for(i = 0;i < dimensions;i++) s[i] += new[i];

    DEBUG(" Aggegated values so far: ");
    DEBUG_ARRAY(dimensions,f->X[sample_idx].dimension);
    DEBUG("\n");
}

/* calculate min,max and avg values after is loaded
 */
static
void calculate_stats(struct forest *f)
{
    int i,j;
    double *s;

    if(f->min == NULL)
    {
        f->min = xmalloc(dimensions * sizeof(double));
        f->max = xmalloc(dimensions * sizeof(double));
    }

    if(f->avg == NULL) f->avg = xmalloc(dimensions * sizeof(double));

    for(i = 0;i < f->X_count;i++)
    {
        s = f->X[i].dimension;

        if(i == 0)
        {
            for(j = 0;j < dimensions;j++)
            {
                f->min[j] = s[j];
                f->max[j] = s[j];
                f->avg[j] = s[j];
            }
        } else
        {
            for(j = 0;j < dimensions;j++)
            {
                if(s[j] < f->min[j]) f->min[j] = s[j];
                if(s[j] > f->max[j]) f->max[j] = s[j];
                f->avg[j] += s[j];
            }
        }
    }
    for(j = 0;j < dimensions;j++) f->avg[j] /= (double) f->X_count;        // turn to average
}

/* Populate sample table with indices to X table
 * If X_sixe is less than samples, take all
 *
 * return  number of samples actually populated
 */
static
int populate_sample(int *s, struct forest *f)
{
    int i,n;

    n = (f->X_count < samples_max) ? f->X_count : samples_max;

    for(i = 0;i < n;i++)
    {
        s[i] = f->X_current;
        f->X_current++;
        if(f->X_current == f->X_count) f->X_current = 0;    // reset counter in case the sample count in X is less than total samples need for a forest
    }
    return n;
}


/*
 * Random normal distribution generator,
 * copied from http://c-faq.com/lib/gaussian.html
 */
static inline
double gaussrand()
{
    static double U, V;
    static int phase = 0;
    double Z;

    if(phase == 0) {
        U = (rand() + 1.) / (RAND_MAX + 2.);
        V = rand() / (RAND_MAX + 1.);
        Z = sqrt(-2. * log(U)) * sin(2. * M_PI * V);
    } else
        Z = sqrt(-2. * log(U)) * cos(2. * M_PI * V);

    phase = 1 - phase;

    return Z;
}

/* Init N cache  */
void init_fast_n_cache()
{
    int i;

    for(i = 0;i < FAST_N_SAMPLES;i++) fast_n_cache[i] = gaussrand();
}

/* calculate random normal distributed number from [0,1]
 */
static inline
double N()
{
    static int next_n=0;

    if(next_n == FAST_N_SAMPLES) next_n = 0;
    return fast_n_cache[next_n++];
}

/* calculate n vector, vector should have dimensions nmber of coordinates
 * returns pointer to vector
 * */
static 
double *calculate_n()
{
    int i;
    static double n[DIM_MAX];

    for(i = 0;i < dimensions;i++)
    {
        n[i] = N();
    }

    return n;
}

/* vector add
   adds b to a
   */
static inline
void v_add(double *a, double *b)
{
    int i;

    for(i = 0;i < dimensions;i++) a[i] += b[i];
}

/* calculate the distance of two samples
 * The sqrt is left out for performance reason
 * */
double v_dist_nosqrt(double *a, double *b)
{
    int i;
    double d = 0.0;

    for(i = 0;i < dimensions;i++) d += POW2(a[i] - b[i]);

    return d;
}

/* calculate the distance of two samples  */
double v_dist(double *a, double *b)
{
    return sqrt(v_dist_nosqrt(a,b));
}


/* vector subtract
   subtracts b from a
   */
static inline
void v_subt(double *a, double *b)
{
    int i;

    for(i = 0;i < dimensions;i++) a[i] -= b[i];
}

void set_centroid_tresshold(double new)
{
    centroid_tresshold = new;
}

/* generate p from sample data.
 * returns pointer to p array.
 * In first nodes are taken from random sample point centered n-sphere having random diameter. Diameter length is proportional to tree heigth (larger at root) and dimension value range / 2.
 * In more deeper nodes the sample cetroid is used as p, this ensures more balanced tree (hopefully)
 */
static
double *generate_p(int sample_count,int *samples,struct sample *X,double heigth_ratio, double *max, double *min)
{
    int i,j;
    int random_sample;
    static int start = 1;
    double *n_vector;
    static double p[DIM_MAX];

    if(heigth_ratio < centroid_tresshold)  // In deeper nodes of tree use sample centroid as p, take only every other sample, speeds things and adds ramdomness
    {
        DEBUG("(centroid)");
    
        for(i = 0;i < dimensions;i++) p[i] = 0.0;

        start = 1 - start;

        for(i = start;i < sample_count;i += 2)
        {
            for(j = 0;j < dimensions;j++)  p[j] += X[samples[i]].dimension[j];
        }

        for(i = 0;i < dimensions;i++) p[i] /= (double) (sample_count >> 1);  // turn to average, divide by sample_count / 2
    } else
    {
        DEBUG("(random)");
        // get a random sample point
        random_sample = ri(0,sample_count - 1);

        // copy random sample to p vector 
        v_copy(p,X[samples[random_sample]].dimension);
    
    
        n_vector = calculate_n();

        // Add adjustment vector

        for(i = 0;i < dimensions;i++) {
            p[i] += n_vector[i] * heigth_ratio * (max[i] - min[i] > 0.0 ? (max[i] - min[i]) / 2.0 : 0.5);   // move sample by adjustment
        }
    }
    

    return p;
}


/* calculate dot from two arrays with size of dimensions
 */

double dot(double *a, double *b)
{
    int i;
    double d = 0.0;

    for(i = 0;i < dimensions;i++)      
    {
        d += a[i]*b[i];
    }
 
    return d;
} 

/* scale a double value. Scaling is done using scale_min and scale_max values
 * with values min...max range
 */
double
scale_double(double value, double range, double scale_min, double min, double max)
{
    if(max == min) return value;

    return range * (value - min) / (max - min) + scale_min;
}


/* calculate average tree height for given sample size n
 *  */
double _c(int n)
{
    if (n < 2) return 0.0;   
    if (n == 2) return 1.0;
    return 2.0 * (log(n - 1) + 0.5772156649) - (2.0 * (double) (n - 1) / (double) n);
}

/* get the average tree height for given sample size n
 *  */
double c(int n)
{
    if(n < FAST_C_SAMPLES) return fast_c_cache[n];
    return _c(n);
}

/*  Init the fast_c_cache table */
void init_fast_c_cache()
{
    int i;

    for(i = 0;i < FAST_C_SAMPLES;i++) fast_c_cache[i] = _c(i);
}

/* make an n vector, 
*/

static
double *make_n_vector()
{
    return calculate_n();
}


/* scale a dim and return pointer to that
 */
double *scale_dimension(double *dim,struct forest *f)
{
    int i;
    double range;
    static double sd[DIM_MAX];

    if(f->scale_range_idx == -1)
    {
        v_copy(sd,dim);
    } else
    {
        range = f->max[f->scale_range_idx] - f->min[f->scale_range_idx];

        for(i = 0;i < dimensions;i++) sd[i] = scale_double(dim[i],range,f->min[f->scale_range_idx],f->min[i],f->max[i]);
    }

    return sd;
}

/*  copy samples array for leaf node
 */
static 
int *copy_samples(int sample_count,int *samples)
{
    int *new = xmalloc(sample_count * sizeof(int));

    memcpy(new,samples,sample_count * sizeof(int));

    return new;
}


/* add nodes to tree
 * returns the index of this node, -1 if end of tree
 */
static 
int add_node(struct forest *f,struct tree *t,int sample_count,int *samples,struct sample *X,int heigth, int heigth_limit)
{
    struct node *this;
    double *p;
    int i,node_index;
    int left_count = 0, rigth_count = 0,new;
    int *left_samples;
    int *rigth_samples;

    if(heigth >= heigth_limit || sample_count < NODE_MIN_SAMPLE) return -1;
    
    DEBUG("    Adding a node with %d samples at heigth %d,",sample_count,heigth);
    
    left_samples = xmalloc(sample_count*sizeof(int));
    rigth_samples = xmalloc(sample_count*sizeof(int));

    if(t->node_count >= t->node_cap) {
        if(t->node_cap == 0) t->node_cap = 16;
        t->node_cap *= 2;
        t->n = xrealloc(t->n,t->node_cap*sizeof(struct node));
    }

    this = &t->n[t->node_count];
    node_index = t->node_count;

    t->node_count++;

    this->sample_count = sample_count;
    this->samples = NULL;

    this->n = v_dup(make_n_vector());

    this->left = -1;
    this->rigth = -1;

    DEBUG(" interception point ");
    p = generate_p(sample_count,samples,X,1.0 - ((double) heigth / (double) heigth_limit),f->max,f->min);

    if(auto_weigth) p = scale_dimension(p,f);

    DEBUG(" p: ");
    DEBUG_ARRAY(dimensions,p);

    this->pdotn = dot(p,this->n);

    for(i = 0;i < sample_count;i++)
    {
        if(dot(sample_dimension(&X[samples[i]]),this->n) < this->pdotn)
        {
            left_samples[left_count] = samples[i];
            left_count++;
        } else
        {
            rigth_samples[rigth_count] = samples[i];
            rigth_count++;
        }
    }

    DEBUG(" left samples: %d, rigth samples: %d\n",left_count,rigth_count);

    if(left_count > 1) 
    {
        new = add_node(f,t,left_count,left_samples,X,heigth + 1,heigth_limit);
        this = &t->n[node_index];
        this->left = new;
    }

    if(rigth_count > 1) 
    {
        new = add_node(f,t,rigth_count,rigth_samples,X,heigth + 1,heigth_limit);
        this = &t->n[node_index];
        this->rigth = new;
    }


    // copy leaf node sample indices for 1-nearest analysis  
    if(this->left == -1 && this->rigth == -1)
    {
        DEBUG("\n    Reached a leaf node at heigth %d",heigth);
        if (nearest && f->avg_sample_dist > 0.0)
        {
            DEBUG(", copying %d samples for nearest distance analysis",sample_count);
            this->samples = copy_samples(sample_count,samples);  // copy samples for nearest distance calculation, c is calculated using sample count adjusted by the distance to nearest sample
        } 
        DEBUG("\n");
    }

    free(left_samples);
    free(rigth_samples);

    return node_index;
}

/* populate one tree
 */
static 
void populate_tree(struct forest *f, struct tree *t,int sample_count,int *samples,struct sample *X,int heigth_limit)
{
    t->first = add_node(f,t,sample_count,samples,X,0,heigth_limit);
}

/* find min...max range to be used in auto scale of dimension attributes 
 * min...max range is the attribute dimension having the smallest maximun absolute value
 *
 * Returns the index of min/max tables. Index points to largest min...max range attribute
 * Returns -1, if no range can be found. In this case all attributes have same value
 */
int find_weigth_scale_idx(double *min,double *max)
{
    int i;
    int idx = -1;
    double max_range = 0.0;
    double s;

    for(i = 0;i < dimensions;i++)
    { 
        s = max[i] - min[i];                         

        if(s > max_range)
        {
            max_range = s;
            idx = i;
        }
    }
    return idx;
}

/* calculate scaled dimension values to be used later
 */
static
void calculate_scaled_dimensions(struct forest *f)
{
    int i;

    for(i = 0;i < f->X_count;i++)
    {
        f->X[i].scaled_dimension = v_dup(scale_dimension(f->X[i].dimension,f));
        DEBUG("   Dimension: ");
        DEBUG_ARRAY(dimensions,f->X[i].dimension);
        DEBUG(" scaled values are: ");
        DEBUG_ARRAY(dimensions,f->X[i].scaled_dimension);
        DEBUG("\n");
    }
}



/* train one forest
 */
#define DIST_AVG(d) ((d) / 1.5 + 1.0 / (2.4 * (d)) - 1.0 / 12.0)
static 
void train_one_forest(int forest_idx)
{
    int i = 0;
    struct forest *f = &forest[forest_idx];
    static int *s = NULL; 
    int sample_count,total_samples = 0;
    double volume;

    DEBUG(" *Training forest %s\n",f->category);
    
    if(s == NULL) s = xmalloc(samples_max * sizeof(int));

    if(!tree_count) f->filter = 1;

    if(f->X_count < SAMPLES_MIN) f->filter = 1;  /*  check the resonable amount of samples */

    if(!f->X_count) return;

    if(f->dim_density == NULL) f->dim_density = xmalloc(dimensions * sizeof(double));
    
    if(auto_weigth)
    {
        f->scale_range_idx = find_weigth_scale_idx(f->min,f->max);              // find larges attribute range
        DEBUG("  Dimension %d has largest value range to be used in value scaling\n",f->scale_range_idx + 1);
        calculate_scaled_dimensions(f);
    }

    volume = 1.0;

    for(i = 0;i < dimensions;i++) 
    {
        f->dim_density[i] = (f->max[i] - f->min[i]) / (double) f->X_count;              // calculate avg dimension density
        if(f->dim_density[i] == 0.0) f->dim_density[i] = 1.0 / (double) f->X_count;     // make sure that density is not zero

        if(!auto_weigth || f->scale_range_idx == -1)
        {
            if(f->max[i] > f->min[i]) volume *= f->max[i] - f->min[i];
        }
    }

    // Calculate the average distance from evenly distributed point to closest points in a hypercube. 
    // This is estimated by dividing the volume by sample count and taking dimensions root, which yields the side length of a cube around 
    // evently distributed points.
    // Side length is multiplyed by sqrt(dimensions / 1.5 + 1 / (2.4 * dimensions) - 1/12)  in order to get app. average distance to all touching (nearest) points.
    // This equation is found to be quite good approximation when comparing real avg. distances to square root of the dimension (tested dimensions 1-19)
    
    if(!auto_weigth || f->scale_range_idx == -1)
    {
        f->avg_sample_dist = sqrt(DIST_AVG((double) dimensions)) *
            pow(volume / (double) ((f->X_count < samples_max) ? f->X_count : samples_max),1.0 / (double) dimensions);    
    } else // if autoscaling the hypercube side is the same as f->max[f->scale_range_idx] - f->min[f->scale_range_idx]
    {
        f->avg_sample_dist = sqrt(DIST_AVG((double) dimensions)) *
            ((f->max[f->scale_range_idx] - f->min[f->scale_range_idx]) / pow((double) ((f->X_count < samples_max) ? f->X_count : samples_max),1.0 / (double) dimensions));
    }

    if(f->filter) return;
    
    if(f->t == NULL) f->t = xmalloc(tree_count * sizeof(struct tree));

    f->X_current = ri(0,f->X_count - 1);           // start at random point

    for(i = 0;i < tree_count;i++)
    {
         sample_count = populate_sample(s,f);
         total_samples += sample_count;

         f->t[i].node_count = 0;
         f->t[i].node_cap = 0;
         f->t[i].n = NULL;
         f->t[i].sample_count = sample_count;
         DEBUG("\n   Populating tree %d with %d samples\n",i,sample_count);
         populate_tree(f,&f->t[i],sample_count,s,f->X,ceil(log2(sample_count)) + 1);
    }

    f->heigth_limit = ceil(log2(total_samples / tree_count)) + 2;
    f->c = c(total_samples / tree_count);    
}

/* Add a regular expression to category filter list
 */
void add_category_filter(char *regexp)
{
    if(cat_filter_count < FILTER_MAX)
    {
        cat_filter[cat_filter_count] = xstrdup(regexp);
        cat_filter_count++;
    }
}

/*  mark forests filtered using
    regular expressions in cat_filter 

    if regex starts with "-v ", then remove it and invert the result
*/
static
void filter_forests()
{
    int i,j,s;
    regex_t reg;

    for(i = 0;i < cat_filter_count;i++)
    {
        s = strncmp(cat_filter[i],"-v ",3) == 0 ? 3 : 0;
         
        if(regcomp(&reg,&cat_filter[i][s],REG_EXTENDED | REG_NOSUB))
        {
            panic("Error in regular expression",&cat_filter[i][s],NULL);
        }

        for(j = 0;j < forest_count;j++)
        {
            if(forest[j].category[0] != '\000')
            {
                if(regexec(&reg,forest[j].category,0,NULL,0) == 0)
                {
                    if(s == 0) forest[j].filter = 1;
                } else
                {
                    if(s != 0) forest[j].filter = 1;
                }
            }
        }
        regfree(&reg);
    }
}


/* build a new forest structure (new=1) or add new samples to existing forest (new=0)
 * make tree in acording make_tree. Tree is needed only if analysing/categorizing or making test data
   */
void
train_forest(FILE *in_stream,int new,int make_tree)
{
    int i,first;
    int value_count;
    int lines = 0;
    int forest_idx;
    static char *values[DIM_MAX];
    static double numval[DIM_MAX];

    

    first = 1;
    now = time(NULL);

    if(in_stream != NULL)
    {
        while(fgets(input_line,INPUT_LEN_MAX,in_stream) != NULL)  // Read data to  memory
        {
            lines++;

            if(header && lines == 1) continue;

            value_count = parse_csv_line(values,DIM_MAX,input_line,input_separator);

            if(first)
            {
                init_dims(value_count);   // init dimensions tables based on first line
                first = 0;
            } 

            forest_idx = select_forest(value_count,values);

            // If we are adding lines to allready loded samples then adjust the line count accordingly
            if(aggregate)
            {
                add_aggregate(&forest[forest_idx],values,value_count);
            } else
            {
                parse_values(numval,values,value_count,0);
                add_to_X(&forest[forest_idx],numval,value_count,0);
            }
        }
    }

    filter_forests();

    // train only once, if new data is only added (!make_tree) no training is run and only new samples are collected
    if(new && make_tree)  
    {
        DEBUG("\n **Starting forest training\n");
        for(i = 0;i < forest_count;i++)
        {
            calculate_stats(&forest[i]);
            train_one_forest(i);
            find_cluster_centers(i);
        }
    }
}



/* Make a test run through forests using points between each dimension min..max range
 * Range is adjusted by test_extension_factor (larger value means larger space)
 * and number of points between max--min is test_sample_interval
 * values are printed as outlier values
 *
 * The sample points as they are are printed too, but with score 0. 
 * So the samples can be plotted with black color (score 0 gives black)
 *
 */
void
test2(FILE *outs,double test_extension_factor,int test_sample_interval)
{
    int i;
    int forest_idx;
    int samples = 0;
    double score, forest_score;
    double zero_len_divisor;
    double *test_dimension;
    double *prev_dimension;
    double *len;             // precalculated max - min value
    double *sidx;            // contains sample number (between 0 - test_sample_interval) for each dimension value, all combinations are processed
    struct forest *f;

    test_dimension = xmalloc(dimensions * sizeof(double));
    prev_dimension = xmalloc(dimensions * sizeof(double));
    len = xmalloc(dimensions * sizeof(double));
    sidx = xmalloc(dimensions * sizeof(double));
    
    for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
    {
        if(!forest[forest_idx].filter)
        {
            calculate_forest_score(forest_idx);
        }
    }

    for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
    {
        if(!forest[forest_idx].filter)
        {
            f = &forest[forest_idx];

            forest_score = get_forest_score(forest_idx);

            zero_len_divisor = 0.0;         // This is used in test range calculation for dimensions having no variation (max == min), the narrowest range in
                                            // data is used in case the narrowest range is below 1 (max - min < 1)

            for(i = 0;i < dimensions;i++) 
            {
                sidx[i] = 0;
                prev_dimension[i] = (double) rand();
                len[i] = f->max[i] - f->min[i];
                if(len[i] > 0.0 && 2.0 / len[i] > zero_len_divisor) zero_len_divisor = 2.0 / len[i];
            }

            if(zero_len_divisor < 2.0) zero_len_divisor = 2.0;

            while(sidx[0] <= test_sample_interval)
            {
                for(i = 0;i < dimensions;i++) 
                {
                    if(len[i] == 0.0)
                    {
                        test_dimension[i] = (1.0 + test_extension_factor) * (2.0 * (double) sidx[i] / (double) test_sample_interval / zero_len_divisor) +
                                            (f->min[i] - ((1.0 + test_extension_factor) / zero_len_divisor));
                    } else
                    {
                        test_dimension[i] = (1.0 + test_extension_factor) * ((double) sidx[i] / (double) test_sample_interval) * len[i] + 
                                            (f->min[i] - (test_extension_factor * len[i]) / 2.0);
                    }
                }

                if(v_cmp(test_dimension,prev_dimension) != 0)
                {
                    score = calculate_score(forest_idx,test_dimension);

                    if(score > forest_score && get_dim_score(forest_idx,test_dimension) > forest_score) print_(outs,score,0,forest_idx,0,NULL,test_dimension,print_string,"sduaxCemgX");

                    v_copy(prev_dimension, test_dimension);
                }

                // next sample
                for(i = dimensions - 1;i >= 0;i--)              
                {
                    if(i && sidx[i] == test_sample_interval)
                    {
                        sidx[i] = 0;
                        sidx[i - 1]++;
                    } else
                    {
                        if(i == dimensions - 1) sidx[i]++;
                    }
                }
            }
                           
            for(samples = 0;samples < TEST_SAMPLES && samples < f->X_count;samples++)
            {
                 print_(outs,0.0,0,forest_idx,0,NULL,f->X[f->X_count <= TEST_SAMPLES ? samples : ri(0,f->X_count - 1)].dimension,print_string,"sduaxCX");
            }
        }
    }
    free(test_dimension);
    free(prev_dimension);
    free(len);
    free(sidx);
}



