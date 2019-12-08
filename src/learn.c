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


static char input_line[INPUT_LEN_MAX];


/* make category sring using values from file and category_idx
 * values are concatenad with semicolon
 * return pointer to string
 */
char *make_category_string(int value_count,char **values)
{
    int i;
    static char c[10240];

    c[0] = '\000';

    for(i = 0;i < categories;i++)
    {
        if(category_idx[i] < value_count)
        {
            strcat(c,values[category_idx[i]]);
            if(i < categories - 1) strcat(c,":");
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
    register int i;
    char *category_string;


    category_string = make_category_string(value_count,values);

    for(i = 0;i < forest_count;i++)
    {
        if (strcmp(category_string,forest[i].category) == 0) return i;
    }

    if(forest_count >= forest_cap) {
       forest_cap += 100;
       forest = xrealloc(forest,forest_cap * sizeof(struct forest));
    }

   forest[forest_count].category = xstrdup(category_string);
   forest[forest_count].X_count = 0;
   forest[forest_count].X_current = 0;
   forest[forest_count].X_cap = 0;   
   forest[forest_count].X = NULL;   
   forest[forest_count].t = NULL;
   forest[forest_count].min = xmalloc(dimensions * sizeof(double));
   forest[forest_count].max = xmalloc(dimensions * sizeof(double));
   forest[forest_count].dim_density = xmalloc(dimensions * sizeof(double));

   forest_count++;

   return forest_count - 1;
}

   

/* parse single sample attribute 
   Possible calculate hash for text attributes later
   now assume all data being float
   */
double parse_dim_attribute(char *value)
{
    return atof(value);
}

/* add one sample to X in selected forest.
  if X is full (== samples_total) implement  reservoir sampling 
 * and check min/max values
   */
void add_to_X(struct forest *f,char **values, int value_count,int sample_number,int saved)
{
    int i,j,sample_idx;

    if(f->X_count >= f->X_cap) {
        f->X_cap += 25000;
        f->X = xrealloc(f->X,f->X_cap * sizeof(struct sample));
    }

    if(f->X_count < samples_total)    // check the samples table size
    {
        sample_idx = f->X_count;
        f->X[sample_idx].dimension = xmalloc(dimensions * sizeof(double));
        f->X_count++;
    } else
    {
        sample_idx = ri(0,sample_number);
        if(sample_idx >= samples_total) return;     // check if old sample should be replaced with this or not
        printf("%d %d\n",sample_number,sample_idx);
    }

    j = 0;

    for(i = 0;i < dimensions;i++)
    {
           // silently ignore missing input row dimensions 
           if(dim_idx[i] < value_count || saved)
           {
               f->X[sample_idx].dimension[j] = parse_dim_attribute(values[saved ? i : dim_idx[i]]);

               // update min/max dimensions
               if(f->X_count == 0)
               {
                   f->min[j] = f->X[sample_idx].dimension[j];
                   f->max[j] = f->X[sample_idx].dimension[j];
               } else
               {
                   if(f->X[sample_idx].dimension[j] < f->min[j]) f->min[j] = f->X[sample_idx].dimension[j];
                   if(f->X[sample_idx].dimension[j] > f->max[j]) f->max[j] = f->X[sample_idx].dimension[j];
               }
               j++;
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

/* duplicate a vector
 */
static
double *v_dup(double *v)
{
    double *new = xmalloc(dimensions * sizeof(double));

    memcpy(new,v,dimensions * sizeof(double));

    return new;
}

/* copy a vector
 */
static
void v_copy(double *t,double *s)
{
    memcpy(t,s,dimensions * sizeof(double));
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
    if(n <= 1) return 0;   
    return 2*(log(n - 1) + 0.5772156649) - (2*(n - 1)/n);
}

        

/* add nodes to tree
 * returns the index of this node, -1 if end of tree
 */
static 
int add_node(struct tree *t,int sample_count,int *samples,struct sample *X,int heigth, int heigth_limit, double *dimension_density)
{
    struct node *this;
    double *p;
    int i,node_index;
    int left_count = 0, rigth_count = 0,new;
    int *left_samples;
    int *rigth_samples;

    if(heigth >= heigth_limit || sample_count < 3) return -1;
    
    left_samples = xmalloc(sample_count*sizeof(int));
    rigth_samples = xmalloc(sample_count*sizeof(int));

    if(t->node_count >= t->node_cap) {
        t->node_cap += 256;
        t->n = xrealloc(t->n,t->node_cap*sizeof(struct node));
    }

    this = &t->n[t->node_count];
    node_index = t->node_count;

    t->node_count++;

    this->level = heigth;
    this->sample_count = sample_count;
    this->n = v_dup(calculate_n());
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

    if(f->X_count < SAMPLES_MIN) return;  /*  check the resonable amount of samples */

    if(f->t == NULL) f->t = xmalloc(tree_count * sizeof(struct tree));

    for(i = 0;i < dimensions;i++) f->dim_density[i] = (f->max[i] - f->min[i]) / (double) f->X_count;   // calculate avg dimension density

    for(i = 0;i < tree_count;i++)
    {
         sample_count = populate_sample(s,f);
         total_samples += sample_count;

         f->t[i].node_count = 0;
         f->t[i].node_cap = 0;
         f->t[i].n = NULL;
         f->heigth_limit = ceil(log2(sample_count + 1) + log(prange_extension_factor + 1));
         populate_tree(&f->t[i],sample_count,s,f->X,f->heigth_limit,f->dim_density);
    }
    f->c = c(total_samples / tree_count);    
    free(s);
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
}
