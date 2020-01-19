/*  
 *    ceif - categorized extended isolation forest
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
#include <math.h>
#include <regex.h>


static char input_line[INPUT_LEN_MAX];

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



/* make category sring using values from file and category_idx
 * values are concatenad with semicolon
 * return pointer to string
 */
char *make_category_string(int value_count,char **values)
{
    int i;
    static char c[10240];

    c[0] = '\000';

    for(i = 0;i < category_idx_count;i++)
    {
        if(category_idx[i] < value_count)
        {
            strcat(c,values[category_idx[i]]);
            if(i < category_idx_count - 1) strcat(c,":");
        }
    }
    return c;
}


/* generate a random integer from range min...max
 */ 
static 
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

    if(i >= 0) return i;

    if(forest_count >= forest_cap) {
       if(forest_cap == 0) forest_cap = 128;
       forest_cap *= 2;
       forest = xrealloc(forest,forest_cap * sizeof(struct forest));
    }

   forest[forest_count].category = xstrdup(category_string);
   forest[forest_count].filter = 0;
   forest[forest_count].X_count = 0;
   forest[forest_count].X_current = 0;
   forest[forest_count].X_cap = 0;   
   forest[forest_count].X = NULL;   
   forest[forest_count].t = NULL;
   forest[forest_count].min = NULL;
   forest[forest_count].max = NULL;
   forest[forest_count].c = 0;
   forest[forest_count].heigth_limit = 0;
   forest[forest_count].avg = xmalloc(dimensions * sizeof(double));
   forest[forest_count].dim_density = xmalloc(dimensions * sizeof(double));

   add_forest_hash(forest_count,category_string);

   forest_count++;

   return forest_count - 1;
}

   
/* copy a vector
 */
static
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


/* parse single sample attribute 
   Possible calculate hash for text attributes later
   now assume all data being float
   */
double parse_dim_attribute(char *value)
{
    return atof(value);
}

/* parse ascii value table to dimension table of type double
 */
void parse_values(double *dim,char **values, int value_count, int saved)
{
    int i;

    for(i = 0;i < dimensions;i++)
    {
        // silently ignore missing input row dimension values
        if(dim_idx[i] < value_count || saved)
        {
            dim[i] = parse_dim_attribute(values[saved ? i : dim_idx[i]]);
        } else
        {
            dim[i] = 0.0;   // use default value for missing dimension vlaue
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


/* add one sample to X in selected forest.
  if X is full (== samples_total) implement  reservoir sampling 
 * and check min/max values
 *
 * saved == 1 means that values are from saved forest file and == 0 means that values are from user given data file.
 * indices are different in those two cases
   */
void add_to_X(struct forest *f,char **values, int value_count,int sample_number,int saved)
{
    int i,sample_idx;
    int first = 0;
    double *s;
    static double new[DIM_MAX];

    parse_values(new,values,value_count,saved); 

    // check if this is duplicate sample. Unique_samples is a value bweteen 0..100. 
    // if value is .e.g 10 then 10 percent of samples are checked for uniqueness
    // check is not done for saved values, only new values are checked

    if(!saved && unique_samples > 0 && ri(0,100) <= unique_samples && duplicate_sample(f,new)) return;

    if(f->X_count >= f->X_cap) {
        if(f->X_cap == 0) f->X_cap = 32;
        f->X_cap *= 2;
        f->X = xrealloc(f->X,f->X_cap * sizeof(struct sample));
    }

    if(f->X_count < samples_total)    // check the samples table size
    {
        if(f->X_count == 0)
        {
            sample_idx = f->X_count;
            f->X[sample_idx].dimension = xmalloc(dimensions * sizeof(double));
        } else
        {
            sample_idx = ri(0,f->X_count - 1);             // ceif does not work well with sorted samples, make sure that  samples are shuffled
            f->X[f->X_count].dimension = v_dup(f->X[sample_idx].dimension);
        }
        f->X_count++;
    } else
    {
        sample_idx = ri(0,sample_number);
        if(sample_idx >= samples_total) return;     // check if old sample should be replaced with this or not
    }

    if(f->min == NULL)
    {
        f->min = xmalloc(dimensions * sizeof(double));
        f->max = xmalloc(dimensions * sizeof(double));
        f->avg = xmalloc(dimensions * sizeof(double));
        first = 1;
    }

    s = f->X[sample_idx].dimension;

    v_copy(s,new);

    for(i = 0;i < dimensions;i++)
    {
        // update min/max dimensions and average
        if(first)
        {
            f->min[i] = s[i];
            f->max[i] = s[i];
            f->avg[i] = s[i];
        } else
        {
            if(s[i] < f->min[i]) f->min[i] = s[i];
            if(s[i] > f->max[i]) f->max[i] = s[i];
            f->avg[i] += s[i];
        }
    }

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

/* calcula random normal distributed number from [0,1]
 */
static inline
double N()
{
    return(gaussrand());
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


/* vector subtract
   subtracts b from a
   */
static inline
void v_subt(double *a, double *b)
{
    int i;

    for(i = 0;i < dimensions;i++) a[i] -= b[i];
}

/* generate p from sample data.
 * returns pointer to p array.
 * p are taken from random sample point centered n-sphere having random diameter
 *
 * diamenter length is proportional to avaerage dimension density, tree heigth (larger at root) and user given prange_extension_factor
 *
 */
static
double *generate_p(int sample_count,int *samples,struct sample *X,double heigth_ratio, double *dimension_density)
{
    int i;
    int random_sample;
    double *n_vector;
    static double p[DIM_MAX];

    // get a random sample point
    random_sample = ri(0,sample_count - 1);
    n_vector = calculate_n();

    // copy random sample to p vector and make adjustment vector using random unit vector n_vector, dimension density, current heigth (= max_heigth at start) and user given factor
    v_copy(p,X[samples[random_sample]].dimension);

    for(i = 0;i < dimensions;i++) {
        p[i] += n_vector[i] * dimension_density[i] * heigth_ratio *  prange_extension_factor;   // move sample by adjust 
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


/* calculate average tree height for given sample size n
 *  */
double c(int n)
{
    if(n <= 2) return 1;   
    return 2*(log(n - 1) + 0.5772156649) - (2*(n - 1)/n);
}

/* make an n vector, If n_vector_adjust is set then try to find n
   which is perpendicular to vector adjust.

   Adjust should be parallel to main data set.

   Adjustment can be used if dimension value ranges have big difference compared to each other
*/

static
double *make_n_vector(double *adjust)
{
    int i;
    double d,min_dot;
    double *n;
    static double best[DIM_MAX];
    
    n = calculate_n();

    if(adjust == NULL) return n;

    v_copy(best,n);
    min_dot = fabs(dot(best,adjust));

    for(i = 1;i < N_ADJUST_COUNT;i++)
    {
        n = calculate_n();
        d = fabs(dot(n,adjust));

        if(d < min_dot)
        {
            min_dot = d;
            v_copy(best,n);
        }
    }

    return best;
}


/* add nodes to tree
 * returns the index of this node, -1 if end of tree
 */
static 
int add_node(struct tree *t,int sample_count,int *samples,struct sample *X,int heigth, int heigth_limit, double *dimension_density)
{
    struct node *this;
    double *p,*adjust = NULL;
    int i,node_index;
    int left_count = 0, rigth_count = 0,new;
    int *left_samples;
    int *rigth_samples;

    if(heigth >= heigth_limit || sample_count < 3) return -1;
    
    left_samples = xmalloc(sample_count*sizeof(int));
    rigth_samples = xmalloc(sample_count*sizeof(int));

    if(t->node_count >= t->node_cap) {
        if(t->node_cap == 0) t->node_cap = 32;
        t->node_cap *= 2;
        t->n = xrealloc(t->n,t->node_cap*sizeof(struct node));
    }

    this = &t->n[t->node_count];
    node_index = t->node_count;

    t->node_count++;

    this->level = heigth;
    this->sample_count = sample_count;

    if(n_vector_adjust)
    {
        adjust = v_dup(generate_p(sample_count,samples,X,1.0 - ((double) heigth / (double) heigth_limit),dimension_density));
        v_subt(adjust,generate_p(sample_count,samples,X,1.0 - ((double) heigth / (double) heigth_limit),dimension_density));
    }

    this->n = v_dup(make_n_vector(adjust));

    if(adjust != NULL) free(adjust);

    this->left = -1;
    this->rigth = -1;
    p = generate_p(sample_count,samples,X,1.0 - ((double) heigth / (double) heigth_limit),dimension_density);
    this->pdotn = dot(p,this->n);

    for(i = 0;i < sample_count;i++)
    {
        if(dot(X[samples[i]].dimension,this->n) < this->pdotn)
        {
            left_samples[left_count] = samples[i];
            left_count++;
        } else
        {
            rigth_samples[rigth_count] = samples[i];
            rigth_count++;
        }
    }

    if(left_count > 1) 
    {
        new = add_node(t,left_count,left_samples,X,heigth + 1,heigth_limit,dimension_density);
        this = &t->n[node_index];
        this->left = new;
    }

    if(rigth_count > 1) 
    {
        new = add_node(t,rigth_count,rigth_samples,X,heigth + 1,heigth_limit,dimension_density);
        this = &t->n[node_index];
        this->rigth = new;
    }

    free(left_samples);
    free(rigth_samples);

    return node_index;
}

/* populate one tree
 */
static 
void populate_tree(struct tree *t,int sample_count,int *samples,struct sample *X,int heigth_limit, double *dimension_density)
{
    t->first = add_node(t,sample_count,samples,X,0,heigth_limit,dimension_density);
}

/* train one forest
 */
static 
void train_one_forest(int forest_idx)
{
    int i = 0;
    struct forest *f = &forest[forest_idx];
    int *s = xmalloc(samples_max * sizeof(int));
    int sample_count,total_samples = 0;

    if(!tree_count) return;

    if(f->X_count < SAMPLES_MIN)  /*  check the resonable amount of samples */
    {
        f->filter = 1;
        return; 
    }

    if(f->t == NULL) f->t = xmalloc(tree_count * sizeof(struct tree));
    if(f->avg == NULL) f->avg = xmalloc(dimensions * sizeof(double));
    if(f->dim_density == NULL) f->dim_density = xmalloc(dimensions * sizeof(double));

    for(i = 0;i < dimensions;i++) 
    {
        f->dim_density[i] = (f->max[i] - f->min[i]) / (double) f->X_count;              // calculate avg dimension density
        if(f->dim_density[i] == 0.0) f->dim_density[i] = 1.0 / (double) f->X_count;     // make sure that density is not zero

        f->avg[i] /= (double) f->X_count;              // turn to average
    }

    f->X_current = ri(0,f->X_count - 1);           // start at random point

    for(i = 0;i < tree_count;i++)
    {
         sample_count = populate_sample(s,f);
         total_samples += sample_count;

         f->t[i].node_count = 0;
         f->t[i].node_cap = 0;
         f->t[i].n = NULL;
         populate_tree(&f->t[i],sample_count,s,f->X,ceil(log2(sample_count + 1) + log(prange_extension_factor + 1)),f->dim_density);
    }

    f->heigth_limit = ceil(log2(total_samples / tree_count + 1) + log(prange_extension_factor + 1));
    f->c = c(total_samples / tree_count);    
    free(s);
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
   */
void
train_forest(FILE *in_stream,int new)
{
    int i,first;
    int value_count;
    int lines = 0;
    int forest_idx;
    char *values[DIM_MAX];

    first = 1;

    if(in_stream != NULL)
    {
        while(fgets(input_line,INPUT_LEN_MAX,in_stream) != NULL)  // Read data to  memory
        {
            lines++;

            if(header && lines == 1) continue;

            value_count = parse_csv_line(values,DIM_MAX,input_line,input_separator);

            if(first)
            {
                first = 0;
                init_dims(value_count);   // init dimensions tables based on first line
            } 

            forest_idx = select_forest(value_count,values);

            // If we are adding lines to allready loded samples then adjust the line count accordingly
            add_to_X(&forest[forest_idx],values,value_count,(header ? (lines - 1) : lines) + (new ? 0 : forest[forest_idx].X_count),0);
        }
    }

    if(!forest_count || !forest[0].X_count) panic("Can't process data with given parameters",NULL,NULL);
    
    for(i = 0;i < forest_count;i++)
    {
        train_one_forest(i);
    }

    filter_forests();
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
    double score;
    double *test_dimension;
    static double len[DIM_MAX];             // precalculated max - min value
    static int sidx[DIM_MAX];               // contains sample number (between 0 - test_sample_interval) for each dimension value, all combinations are processed
    struct forest *f;

    test_dimension = xmalloc(dimensions * sizeof(double));

    for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
    {
        if(!forest[forest_idx].filter)
        {
            f = &forest[forest_idx];

            for(i = 0;i < dimensions;i++) 
            {
                sidx[i] = 0;
                len[i] = f->max[i] - f->min[i];
            }

            while(sidx[0] <= test_sample_interval)
            {
                for(i = 0;i < dimensions;i++) 
                {
                    test_dimension[i] = (1.0 + test_extension_factor) * ((double) sidx[i] / (double) test_sample_interval) * len[i] + (f->min[i] - (test_extension_factor * len[i]) / 2.0);
                }

                score = calculate_score(forest_idx,test_dimension);

                if(score >= outlier_score) print_test(outs,score,forest_idx,test_dimension);

                // next sample
                for(i = dimensions - 1;i >= 0;i--)              
                {
                    if(sidx[i] == test_sample_interval)
                    {
                        if(i)
                        {
                            sidx[i] = 0;
                            sidx[i - 1]++;
                        }
                    } else
                    {
                        if(i == dimensions - 1) sidx[i]++;
                    }
                }
            }
                           
            for(samples = 0;samples < TEST_SAMPLES && samples < f->X_count;samples++)
            {
                 print_test(outs,0.0,forest_idx,f->X[ri(0,f->X_count - 1)].dimension);
            }
        }
    }
    free(test_dimension);
}



