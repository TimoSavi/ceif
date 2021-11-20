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
#include <time.h>
#include <float.h>

#define INPUT_LEN_MAX 1048576

static char input_line[INPUT_LEN_MAX];
static int first = 1;
static char *float_format = "%.*f";

/*struct for finding sample cluster centers 
 */
struct sample_score
{
    size_t idx;         // index to X array
    double score;       // sample X[idx] score
};

/* find a forest for a data row
 *
 * check filter and update analyzed only if filter_on is set
 * 
 * returns index to forest table, -1 if not found
 */
static 
int find_forest(int value_count,char **values, int filter_on)
{
    int i;
    char *category_string;

    category_string = make_category_string(value_count,values);

    i = search_forest_hash(category_string);

    if(filter_on)
    {
        if(i == -1 || forest[i].filter) return -1;

        if(!forest[i].analyzed) forest[i].analyzed = 1;
    }

    return i;
} 


/* Search the nearest training sample for a analyzed point a
   returns the shortest relative distance

   relative distance is 1 if the actual distance is the same as forest average sample distance
   relative distance < 1 if the actual distance is smaller than forest average sample distance but never smaller than MIN_REL_DIST
   relative distance > 1 if the actual distance is larger than forest average sample distance

   a is assumed be scaled in case auto scaling (auto_weigth)
 */
#define MIN_REL_DIST 0.05
double nearest_rel_distance(double *a, int sample_count,int *samples,struct forest *f)
{
    int i;
    double distance,d; 

    distance = v_dist_nosqrt(a,sample_dimension(&f->X[samples[0]]));

    for(i = 1;i < sample_count;i++)
    {
        d = v_dist_nosqrt(a,sample_dimension(&f->X[samples[i]]));

        if(d < distance) distance = d;
    }

    distance = sqrt(distance) / f->avg_sample_dist + MIN_REL_DIST;

    return distance;
}


/* search through nodes, 
 * returns the heigth from last node
 *
 * if nearest is true the shortes distance to node samples is calculated.
 * The distance is used to adjust node sample count and node sample c value is recalculed
 * Shorter distance compared to forest average, the larger adjusted sample_count is used 
 *
 */
static 
double search_last_node(struct forest *f,int this_idx,struct node *n,double *dimension,int heigth)
{
    struct node *this = &n[this_idx];

    if(this->left == -1 && this->rigth == -1)
    {
        DEBUG("\n    Reached a leaf node at heigth %d with %d samples",heigth,this->sample_count);
        if(nearest && f->avg_sample_dist > 0.0)
        {
            double rel_dist = nearest_rel_distance(dimension,this->sample_count,this->samples,f);

            DEBUG(", Calculated nearest relative distance to be: %f\n",rel_dist); 
            return (double) heigth + c((double) this->sample_count / rel_dist);
        }
        DEBUG("\n");
        return (double) heigth + c(this->sample_count);
    }

    DEBUG("    Reached a node at heigth %d with %d samples\n",heigth,this->sample_count);

    if(dot(dimension,this->n) < this->pdotn)
    {
        if(this->left == -1) return (double) heigth;
        return search_last_node(f,this->left,n,dimension,heigth + 1);
    } else
    {
        if(this->rigth == -1) return (double) heigth;
        return search_last_node(f,this->rigth,n,dimension,heigth + 1);
    }
}
    

/*
 * travel through a tree with a test case.
 * return the path length
 */
static 
double calculate_path_length(struct forest *f,struct tree *t,double *dimension)
{
    return search_last_node(f,t->first,t->n,dimension,0);
}


/* calculates score for a data row in given forest
 * returns score
 */
double _score(int forest_idx,double *dimension)
{
    int i;
    struct forest *f = &forest[forest_idx];
    double path_length = 0.0;

    DEBUG("\n Calculating score in forest %s for values: ",f->category);
    DEBUG_ARRAY(dimensions,dimension);
    DEBUG("\n");

    if(f->t != NULL)
    {
        for(i = 0;i < tree_count;i++)
        {
            DEBUG("\n    Scan tree %d\n",i + 1);
            path_length += calculate_path_length(f,&f->t[i],dimension);
            DEBUG("    Average path length now: %f\n",path_length / (i + 1));
        }
    }

    path_length = path_length / tree_count;  // turn to average 

    return (1.0/pow(2,path_length/f->c));
}


/* Calculates max score for a forest
 * This is done making 3^dimensions combinations of +MAX_DIM,-MAX_DIM and 0 
 * The  largest score is returned
 * Number of tested dimensions is limited by LIMIT_DIM for performance reasons
 */
#define MAX_DIM_VALUE (1e+100)
#define pwrtwo(x) ((unsigned int) 1 << (x))
#define LIMIT_DIM 8

double calculate_max_score(int forest_idx)
{
    unsigned int i,j,k;
    unsigned int bitmap1,bitmap2,state,lim_dim; 
    double *dim;
    double score,max_score = 0.0;
    int l;

    dim = xmalloc(dimensions * sizeof(double));

    lim_dim = (dimensions > LIMIT_DIM) ? LIMIT_DIM : dimensions;
    
    for(l = lim_dim;l < dimensions;l++) dim[l] = MAX_DIM_VALUE;   // Init possible rest values with +max

    for(i = 0;i < pwrtwo(lim_dim);i++)
    {
        for(j = 0;j <  pwrtwo(lim_dim);j++)
        {
            bitmap1 = i;        // Two bitmaps in order to get state values 0..2, for 0, +max and -max
            bitmap2 = j;

            if(!(bitmap1 & bitmap2))  // skip cases where there is both bits set in the same position (meaning state value 3, which is not allowed here)
            {
                for(k = 0;k < lim_dim;k++)
                {
                    state = ((bitmap1 & 1) << 1) | (bitmap2 & 1);
                    switch(state)
                    {
                        case 0:
                            dim[k] = 0.0;
                            break;
                        case 1:
                            dim[k] = MAX_DIM_VALUE;
                            break;
                        case 2:
                            dim[k] = -MAX_DIM_VALUE;
                            break;
                    }
                    bitmap1 >>= 1;
                    bitmap2 >>= 1;
                }

                score = _score(forest_idx,dim);
                if(score > max_score) max_score = score;
            }
        }
    }

    free(dim);

    return max_score;
}

/* calculate scaled score. Forest min (from sample having lowest score) and max range (found using a dim with "big" values) 
 * is used to scale score to range 0...1
 */
double calculate_score_scale(int forest_idx,double *dimension)
{
    double score;

    score = scale_double(_score(forest_idx,dimension),1.0,0.0,forest[forest_idx].min_score,forest[forest_idx].max_score);

    if(score < 0.0) score = 0.0;
    if(score > 1.0) score = 1.0;

    return score;
}


/* calculate score, scale if scale_score is set
 */
double calculate_score(int forest_idx,double *dimension)
{
    struct forest *f = &forest[forest_idx];
    double *dim = auto_weigth ? scale_dimension(dimension,f) : dimension;

    return scale_score ? calculate_score_scale(forest_idx,dim) : _score(forest_idx,dim);
}

/* Try to find out how each dimension value effects to outlier score
 * This is done by assingning each dimension value to each cluster center dim. and
 * calculating the score. If returned score is high for all clusters, then it can be assumed that this particular dim
 * is an outlier value. 
 * Minimum scores for each dim is stored to array
 *
 * Returns the result array
 */
static
double * find_dimension_score(int forest_idx,double *dimension)
{
    static double result[DIM_MAX];
    static double test[DIM_MAX];
    double min,score;
    struct forest *f;
    int i,j;

    f = &forest[forest_idx];

    for(i = 0;i < dimensions;i++)
    {
        min = 1.0;
        for(j = 0;j < f->cluster_count;j++)
        {
            v_copy(test,f->X[f->cluster_center[j]].dimension);

            test[i] = dimension[i];
            score = calculate_score(forest_idx,test);

            if(score < min) min = score;
        }
        result[i] = min;
    }
    return result;
}


/* Returns the score for dimensions attributes given by option -G. 
 * This can be used to check only certain attributes or if certain attributes needs to be left out from 
 * outlier analysis
 *
 * This is done by assigning those attribute values to each cluster center and analysing the score.
 * If several cluster centers then the smallest score is taken
 * This tries to find out which score those attributes cause alone.
 *
 * Returns score value or 2.0 if cluster centers or attribute indices given with option -G are not available
 */
double get_dim_score(int forest_idx,double *dimension)
{
    struct forest *f = &forest[forest_idx];
    static double test[DIM_MAX];
    double score,min;
    int i,j;

    min = 2.0;

    if(score_idx_count && cluster_relative_size > 0.0)
    {
        for(i = 0;i < f->cluster_count;i++)
        {
            v_copy(test,f->X[f->cluster_center[i]].dimension);

            for(j = 0;j < dimensions;j++) if(check_idx(j,score_idx_count,score_idx)) test[j] = dimension[j];

            score = calculate_score(forest_idx,test);
            if(score < min) min = score;
        }
    }

    return min;
}
            

/* calculate sample score. Select right dimension using auto_weigth and scale score
 */
double sample_score_scale(int forest_idx,struct sample *s)
{
    double *dim = sample_dimension(s);

    return scale_score ? calculate_score_scale(forest_idx,dim) :  _score(forest_idx,dim);
}


/* calculate sample score. Select right dimension using auto_weigth. Do not implement score scaling
 */
double sample_score(int forest_idx,struct sample *s)
{
    return _score(forest_idx,sample_dimension(s));
}


/* make a RGB value using score.
 * can be used for plotting the results
 * value mapping:
 * 1 : red
 * 0.5: green
 * 0: blue
 * and colors between them
 * if score == 0,return 0 (black)
 */
static 
unsigned int score_to_rgb(double score)
{
    unsigned int red,blue,green;

    if(score == 0.0) return 0;

    green = score < 0.5 ? (unsigned int) (2 * 0xff * score) : (unsigned int) (0xff * 2 * (1 - score)); 
    red = score > 0.5 ? (unsigned int) (2 * 0xff * (score - 0.5)) : 0; 
    blue = score < 0.5 ? (unsigned int) (2 * 0xff * (0.5 - score)) : 0; 
    
    return (red << 16) + (green << 8) + blue;
}                      

/* make label sring using values from file and label_idx
 *   values are concatenad with semicolon
 *   return pointer to string
 */
static 
char *make_label_string(int value_count,char **values)
{
    int i;
    char *label;

    label = make_separated_string(NULL,0);

    if(label_idx_count)
    {
        for(i = 0;i < label_idx_count - 1;i++)
        {
            if(label_idx[i] < value_count)
            {
                make_separated_string(values[label_idx[i]],label_separator);
            }
        }
        make_separated_string(values[label_idx[label_idx_count - 1]],0);
    }

    return label;
}

/* make category sring using values from file and category_idx
 * values are concatenad with category_separator
 * return pointer to string
 */
char *make_category_string(int value_count,char **values)
{
    int i;
    char *category;

    category = make_separated_string(NULL,0);

    if(category_idx_count)
    {
        for(i = 0;i < category_idx_count - 1;i++)
        {
            if(category_idx[i] < value_count)
            {
                make_separated_string(values[category_idx[i]],category_separator);
            }
        }
        make_separated_string(values[category_idx[category_idx_count - 1]],0);
    }

    return category;
}

/* prints escaped char
 * returns number of characters consumed
 */
static 
size_t print_escaped(FILE *outs, char *esc)
{
    size_t consumed = 0;

    if(*esc == '\\' && esc[1] != '\000')
    {
        esc++;
        consumed = 2;
        switch(*esc)
        {
            case 't':
                fprintf(outs,"%c",'\t');
                break;
            case 'n':
                fprintf(outs,"%c",'\n');
                break;
            case '\\':
                fprintf(outs,"%c",'\\');
                break;
            case '"':
                fprintf(outs,"%c",'"');
                break;
            case '\'':
                fprintf(outs,"%c",'\'');
                break;
            default:
                consumed = 0;
                break;
        }
    }
    return consumed;
}

/* print list separator for dimension list
 */
static 
void print_dim_list_separator(FILE *outs,int index)
{
    if(index < dimensions - 1) fputc(list_separator,outs);
}


/* Print something
 * printing is done using printf style string and %-directives
 */
void print_(FILE *outs, double score, int lines,int forest_idx,int value_count,char **values,double *dimension,char *format,char *directives)
{
    int i;
    size_t consumed;
    char *c = format;
    char *d;
    char outstr[100];
    struct tm *tmp;
    double *earray = NULL;

    while(*c != '\000')
    {
        if(*c == '%' && c[1] != '\000' && (strchr(directives,c[1]) != NULL || strchr(":.%",c[1]) != NULL))
        {
            c++;
            switch(*c)
            {
                case 'r':
                    fprintf(outs,"%d",lines);
                    break;
                case 'n':
                    fprintf(outs,"%d",forest[forest_idx].total_rows);
                    break;
                case 'o':
                    fprintf(outs,"%d",forest[forest_idx].analyzed_rows);
                    break;
                case 'h':
                    fprintf(outs,"%d",forest[forest_idx].high_analyzed_rows);
                    break;
                case 's':
                    fprintf(outs,"%f",score);
                    break;
                case 'S':
                    fprintf(outs,"%f",forest[forest_idx].test_average_score);
                    break;
                case 'c':
                    fprintf(outs,"%s",make_category_string(value_count,values));
                    break;
                case 'l':
                    fprintf(outs,"%s",make_label_string(value_count,values));
                    break;
                case 'm':
                    if(print_dimension != NULL)
                    {
                        for(i = 0;i < dimensions;i++)
                        {
                            d = print_dimension;
                            while(*d != '\000')
                            {
                                if(*d == '%' && d[1] != '\000' && strchr("daei",d[1]) != NULL)
                                {
                                    d++;
                                    switch(*d)
                                    {
                                        case 'd':
                                            if(check_idx(dim_idx[i],text_idx_count,text_idx))
                                            {
                                                if(values != NULL) fprintf(outs,"%s",values[dim_idx[i]]);
                                            } else
                                            {
                                                if(dimension != NULL) *printf_format != '\000' ? fprintf(outs,printf_format,dimension[i]) : fprintf(outs,float_format,decimals,dimension[i]);
                                            }
                                            break;
                                        case 'a':
                                            if(forest_idx > -1) *printf_format != '\000' ? fprintf(outs,printf_format,forest[forest_idx].avg[i]) : fprintf(outs,float_format,decimals,forest[forest_idx].avg[i]);
                                            break;
                                        case 'e':
                                            if(forest_idx > -1 && !forest[forest_idx].filter &&  cluster_relative_size > 0.0 && dimension != NULL)
                                            {
                                                if(earray == NULL) earray = find_dimension_score(forest_idx,dimension);
                                                fprintf(outs,"%f",earray[i]);
                                            }
                                            break;
                                        case 'i':
                                            fprintf(outs,"%d",i + 1);
                                            break;
                                    }
                                    d++;
                                } else if((consumed = print_escaped(outs,d))) // yes, it's an assignment
                                {
                                    d += consumed;
                                }  else
                                {
                                    fprintf(outs,"%c",*d);
                                    d++;
                                }
                            }
                            print_dim_list_separator(outs,i);
                        }
                    }
                    break;
                case 'd':
                    for(i = 0;i < dimensions;i++)
                    {
                        if(check_idx(dim_idx[i],text_idx_count,text_idx) && values != NULL)
                        {
                            fprintf(outs,"%s",values[dim_idx[i]]);
                        } else
                        {
                            *printf_format != '\000' ? fprintf(outs,printf_format,dimension[i]) : fprintf(outs,float_format,decimals,dimension[i]);
                        }
                        print_dim_list_separator(outs,i);
                    }
                    break;
                case 'e':
                    if(cluster_relative_size > 0.0)
                    {
                        if(earray == NULL) earray = find_dimension_score(forest_idx,dimension);
                        for(i = 0;i < dimensions;i++)
                        {
                            fprintf(outs,"%f",earray[i]);
                            print_dim_list_separator(outs,i);
                        }
                    }
                    break;
                case 'a':
                    for(i = 0;i < dimensions;i++)
                    {
                        *printf_format != '\000' ? fprintf(outs,printf_format,forest[forest_idx].avg[i]) : fprintf(outs,float_format,decimals,forest[forest_idx].avg[i]);
                        print_dim_list_separator(outs,i);
                    }
                    break;
                case 'v':
                    for(i = 0;i < value_count;i++)
                    {
                        fprintf(outs,"%s",values[i]);
                        if(i < value_count - 1) fputc(list_separator,outs);
                    }
                    break;
                case 'x':
                    fprintf(outs,"%06X",score_to_rgb(score));
                    break;
                case 'C':
                    fprintf(outs,"%s",forest[forest_idx].category);
                    break;
                case 't':
                    tmp = localtime(&forest[forest_idx].last_updated);

                    if(tmp != NULL)
                    {
                        outstr[0] = '\000';
                        strftime(outstr, sizeof(outstr), "%c", tmp);
                        fprintf(outs,"%s",outstr);
                    }
                    break;
                case ':':
                    fprintf(outs,"%c",category_separator);
                    break;
                case '.':
                    fprintf(outs,"%c",label_separator);
                    break;
                case '%':
                    fprintf(outs,"%c",'%');
                    break;
            }
            c++;
        } else if((consumed = print_escaped(outs,c))) // yes, it's an assignment
        {
            c += consumed;
        } else
        {
            fprintf(outs,"%c",*c);
            c++;
        }
    }
    if(*format) fprintf(outs,"%c",'\n');
}


/* check if index is in a index table.
 */
int check_idx(int index,int index_count, int *idx_table)
{
    register int i;

    for(i = 0;i < index_count;i++) if(index == idx_table[i]) return 1;

    return 0;
}


/*  mark category and label dims as non dimensions dims  and populate dim_idx and
 *  category_idx tables
 *  */
void
init_dims(int value_count)
{
    int i,d;

    d = 0;

    for(i = 0;i < value_count;i++) 
    {
        if((!check_idx(i,ignore_idx_count,ignore_idx) || check_idx(i,include_idx_count,include_idx)) &&
            !check_idx(i,category_idx_count,category_idx) && 
            !check_idx(i,label_idx_count,label_idx) &&
            d < DIM_MAX)
        {
            dim_idx[d] = i;
            d++;
        }
    }

    if(!dimensions) dimensions = d;   // If number of dims allready read from saved file, dont mess that
}

/* aggregate values to forest summary
 */
static 
void aggregate_values(int forest_idx,double *sample)
{
    int i;
    struct forest *f = &forest[forest_idx];

    if(f->summary == NULL)
    {
        f->summary = xmalloc(dimensions * sizeof(double));
        for(i = 0;i < dimensions;i++) f->summary[i] = sample[i];
    } else
    {
        for(i = 0;i < dimensions;i++) f->summary[i] += sample[i];
    }
}

/* comparison function for score table sorting
 */
static
int pscore_cmp(const void * a, const void * b)
{
    if(*(double*)a < *(double*)b) return -1;
    if(*(double*)a > *(double*)b) return 1;
    return 0;
}

/* calculate percentage based score for a forest
 * first all sample scores are calculated and sorted
 * the forest score is the score at point x% of sorted scores, 
 */
void calculate_forest_percentage_score(int forest_idx)
{
    double *all_scores;
    int i;
    struct forest *f;

    f = &forest[forest_idx];
        
    if(f->filter) return;

    all_scores = xmalloc(sizeof(double) *  f->X_count);

    for(i = 0;i < f->X_count;i++) all_scores[i] = sample_score(forest_idx,&f->X[i]);

    qsort(all_scores,f->X_count,sizeof(double),pscore_cmp);

    f->percentage_score = all_scores[(size_t) ((double) (f->X_count - 1) * (outlier_score / 100.0))];

    free(all_scores);
}
    
    


/* remove outlier
 * For each non filterd forest the sample with the highest  score is removed from 
 * samples.
 *
 * Forest must be saved with option -w if the if cleaned forest data is need afterwards
 */
void remove_outlier()
{
    int forest_idx,i,outlier_idx;
    double score,max_score;
    struct forest *f;

    for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
    {
        f = &forest[forest_idx];

        if(f->filter || f->X_count <= SAMPLES_MIN) continue;    // Do not remove samples below limit, renders forest filtered

        outlier_idx = -1;
        max_score = 0.0;

        for(i = 0;i < f->X_count;i++)
        {
            score = sample_score(forest_idx,&f->X[i]);
            if(score > max_score)
            {
                max_score = score;
                outlier_idx = i;
            }
        }

        if(outlier_idx >= 0)
        {
            free(f->X[outlier_idx].dimension);

            for(i = outlier_idx;i < f->X_count - 1;i++)
            {
                f->X[i] = f->X[i + 1];
            }
            f->X_count--;
        }
    }
}

/* calculate the range of scores for a forest this is used to scale scores to 0..1 range in analysis
 * Min score is considered to be the smallest score among samples
 * Max score is tested using huge values for each dimension attributes (calculate_max_score)
 * and in order to make sure that max score is really max adjust it by MAX_SCORE_ADJUST
 *
 */
#define MAX_SCORE_ADJUST 1.01
void calculate_sample_score_range(int forest_idx)
{
    int i;
    struct forest *f;
    double score;

    f = &forest[forest_idx];

    if(f->filter || f->min_score < 1.0) return;   // Range is calculated if f->min_score < 1.0
     
    f->min_score = 1.0;

    for(i = 0;i < f->X_count;i++)
    {
        score = sample_score(forest_idx,&f->X[i]);

        if(score < f->min_score) f->min_score = score;
    }

    f->max_score = calculate_max_score(forest_idx) * MAX_SCORE_ADJUST;   

    if(f->max_score > 1.0) f->max_score = 1.0;
}

/* Return score for a forest, special scores are handled here too
 */
inline 
double get_forest_score(int forest_idx)
{
    if(percentage_score) return forest[forest_idx].percentage_score;
    return outlier_score;
}

/* Calculate appropriate automatic score for a forest
 */
inline 
void calculate_forest_score(int forest_idx)
{
    if(!forest[forest_idx].filter) DEBUG("\n** Calculating forest score\n");

    if(scale_score) 
    {
        calculate_sample_score_range(forest_idx);
    } else if(forest[forest_idx].percentage_score == 0.0 && percentage_score)
    {
        calculate_forest_percentage_score(forest_idx);
    }
}

/* qsort comparison function for cluster center analysis
 */
static 
int cluster_score_cmp(const void *a, const void *b)
{
    const struct sample_score *e1 = a;
    const struct sample_score *e2 = b;
    
    if(e1->score < e2->score)
    {
        return -1;
    } else if(e1->score > e2->score)
    {
        return 1;
    } else
    {
        return 0;
    }
}



/* Find forest cluster centers
 * Centers are found using method:
 * - sort all samples by score
 * - take first n scores (having the smallest score -> cluster centers are among these). n = 97.5% of samples
 * - calculate largest distance from sample having smallest score to farthest sample among first n scores
 * - Take first sample and remove all samples being near (samples in the same cluster) to first sample
 * - Take next sample as next cluster center and remove nearby samples
 * - continue until all n samples are done
 */
#define CLUSTER_SAMPLE_DIV 0.975
#define CLUSTER_CHECK 0
#define CLUSTER_DONE 1
void find_cluster_centers(int forest_idx)
{
    struct forest *f;
    struct sample_score *samples;
    static double min[DIM_MAX],max[DIM_MAX];
    double dist,longest_dist = 0.0,same_cluster_dist;
    size_t *status;
    int samples_to_analyze,cluster_samples = 0;
    int i,j,current,next,min_cluster_sample_count;
    int cluster_sample_count[CLUSTER_MAX];

    f = &forest[forest_idx];
    
    f->cluster_count = 0;

    if(f->filter || cluster_relative_size == 0.0) return;

    samples_to_analyze = (int) (CLUSTER_SAMPLE_DIV * (double) f->X_count);

    samples = xmalloc(sizeof(struct sample_score) *  f->X_count);
    status = xmalloc(sizeof(size_t) * samples_to_analyze);

    for(i = 0;i < CLUSTER_MAX;i++) cluster_sample_count[i] = 0;

    for(i = 0;i < f->X_count;i++)
    {
        samples[i].idx = i;
        samples[i].score = sample_score(forest_idx,&f->X[i]);
    }

    qsort(samples,f->X_count,sizeof(struct sample_score),cluster_score_cmp);

    for(i = 0;i < samples_to_analyze;i++) status[i] = CLUSTER_CHECK;

    // We use nosqrt distances for performance reasons
    //
    // First find longest distance. Here is assumed that the sample with lowest score (first one) is the center of all samples
    // distance to the farthest sample is measured from that
    for(i = 0;i < dimensions;i++) 
    {
        min[i] = f->X[samples[0].idx].dimension[i];
        max[i] = f->X[samples[0].idx].dimension[i];
    }

    for(i = 1;i < samples_to_analyze;i++)
    {
        for(j = 0;j < dimensions;j++)
        {
            if(f->X[samples[i].idx].dimension[j] < min[j]) min[j] = f->X[samples[i].idx].dimension[j];
            if(f->X[samples[i].idx].dimension[j] > max[j]) max[j] = f->X[samples[i].idx].dimension[j];
        }
    }

    longest_dist =  v_dist_nosqrt(max,min);

    // First cluster is around the sample having lowest score
    f->cluster_center[f->cluster_count] = samples[0].idx;
    f->cluster_count++;

    // samples are considered to be in same cluster if their distance is smaller than same_cluster_dist calculated here
    // cluster_relative_size default value is in cluster_relative_size, use power of two because square root is not used in distance calculation
    same_cluster_dist = POW2(cluster_relative_size) * longest_dist;
    f->cluster_radius = cluster_relative_size * sqrt(longest_dist);


    current = 0;   // cluster center point, index to X
    status[0] = CLUSTER_DONE;

    while(f->cluster_count < CLUSTER_MAX)
    {
        next = -1;

        for(i = 0;i < samples_to_analyze;i++)
        {
            if(status[i] == CLUSTER_CHECK)
            {
                dist = v_dist_nosqrt(f->X[samples[current].idx].dimension,f->X[samples[i].idx].dimension);    // check distance
                if(dist <= same_cluster_dist)
                {
                    status[i] = CLUSTER_DONE;                     // if in cluster then  mark as done
                    cluster_sample_count[f->cluster_count - 1]++;
                    cluster_samples++;
                }
                else if(next == -1)  // add new cluster, but check that the distance for existing clusters is longer that 2 * f->cluster_radius
                {
                    next = i;
                    for(j = 0;j < f->cluster_count && next > -1;j++)
                        if(v_dist_nosqrt(f->X[f->cluster_center[j]].dimension,f->X[samples[next].idx].dimension) < POW2(2.0 * f->cluster_radius)) next = -1;
                    if(next > -1) status[i] = CLUSTER_DONE;                     // mark as done
                }
            }
        }

        if(next == -1) break;
        current = next;
        f->cluster_center[f->cluster_count] = samples[current].idx;
        f->cluster_count++;
    }

    f->cluster_coverage = (double) cluster_samples / (double) samples_to_analyze;

    // Remove clusters having small number of samples
    // Minimum number is app. the number of samples / 2 in clusters when they are evenly distributed for each cluster
    min_cluster_sample_count = (cluster_samples / f->cluster_count) / 2;

    for(i = 0;i < f->cluster_count;i++)
    {
        while(f->cluster_count > i && cluster_sample_count[i] < min_cluster_sample_count)
        {
            for(j = i;j < f->cluster_count - 1;j++)
            {
                f->cluster_center[j] = f->cluster_center[j + 1];
                cluster_sample_count[j] = cluster_sample_count[j + 1];
            }
            f->cluster_count--;
        }
    }

    free(samples);
    free(status);
}

/* Check if analyzed rows should be sampled and 
 * implement reservoir sampling if number of rows read is
 * larger than analyze_sampling_count
 * return true if row should be analyzed
 */
static
int take_this_row(int total_rows)
{
    if(analyze_sampling_count && total_rows > analyze_sampling_count && ri(1,total_rows) > analyze_sampling_count) return 0;
    return 1;
}
/* analyze data from file. 
 * All lines are analyzed against loaded forest/tree data
 * and print anomalies (having score > outlier_score) using printing mask
 */
void
analyze(FILE *in_stream, FILE *outs,char *not_found_format,char *average_format)
{
    int value_count;
    int lines = 0;
    int forest_idx;
    char *values[DIM_MAX];
    double *dimension = NULL;
    double score,forest_score;
    
    if(!first) dimension =  xmalloc(dimensions * sizeof(double));

    DEBUG("*** Starting analysis\n");
        
    while(fgets(input_line,INPUT_LEN_MAX,in_stream) != NULL) 
    {
        lines++;

        if(header && lines == 1) continue;

        value_count = parse_csv_line(values,DIM_MAX,input_line,input_separator);

        if(first && value_count) 
        {
            init_dims(value_count);
            dimension =  xmalloc(dimensions * sizeof(double));
            first = 0;
        }

        if(value_count)
        { 
            parse_values(dimension,values,value_count,0);

            forest_idx = find_forest(value_count,values,1);

            if(forest_idx >= 0)
            {
                if(!forest[forest_idx].total_rows) calculate_forest_score(forest_idx);

                forest[forest_idx].total_rows++;

                if(aggregate)
                {
                    aggregate_values(forest_idx,dimension);
                } else
                {
                    if(take_this_row(forest[forest_idx].total_rows))   // check if analyzed rows are reservoir sampled
                    {
                        forest[forest_idx].analyzed_rows++;
                        
                        DEBUG("\n *Calculate score for a dimension\n");

                        score = calculate_score(forest_idx,dimension);

                        if(average_format != NULL) forest[forest_idx].test_average_score += score;

                        forest_score =  get_forest_score(forest_idx);

                        if(score > forest_score && get_dim_score(forest_idx,dimension) > forest_score)
                        {
                            forest[forest_idx].high_analyzed_rows++;
                            print_(outs,score,lines,forest_idx,value_count,values,dimension,print_string,"rscldavxCtnohem");
                        }
                    }
                }
            } else
            {
                if(not_found_format != NULL && find_forest(value_count,values,0) == -1) print_(outs,0,lines,-1,value_count,values,dimension,not_found_format,"dvclm");
            }

        }
    }

    if(aggregate)
    {
        for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
        {
            if(forest[forest_idx].summary != NULL && !forest[forest_idx].filter)
            {
                forest[forest_idx].analyzed_rows = 1;
                        
                score = calculate_score(forest_idx,forest[forest_idx].summary);
                    
                if(average_format != NULL) forest[forest_idx].test_average_score = score;

                 forest_score =  get_forest_score(forest_idx);
                
                if(score > forest_score && get_dim_score(forest_idx,forest[forest_idx].summary) > forest_score)
                {
                    forest[forest_idx].high_analyzed_rows++;
                    print_(outs,score,0,forest_idx,0,NULL,forest[forest_idx].summary,print_string,"rsdaxCtnohem");
                }
            }
        }
    }

    if(average_format != NULL)
    {
        for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
        {
            if(!forest[forest_idx].filter && forest[forest_idx].analyzed_rows > 0)
            {
                forest[forest_idx].test_average_score /= forest[forest_idx].analyzed_rows;
                print_(outs,get_forest_score(forest_idx),lines,forest_idx,0,NULL,NULL,average_format,"sraxCthSno");
            }
        }
    }
    
    if(dimension != NULL) free(dimension);
}


/* Categorize dimensions
 * All lines are analyzed against loaded forest/tree data
 * All forests are analyzed and a forest having lowest anomaly score is selected as category forest
 * If score_limit then do not print cases where lowest score is higher than forest outlier score
 * Note that scaled score is used in order to get more comparable scores between forests
 */
void
categorize(FILE *in_stream, int score_limit, FILE *outs)
{
    int i,j;
    int value_count;
    int lines = 0;
    int forest_idx;
    int best_forest_idx;
    int save_scale_score = scale_score;
    char *values[DIM_MAX];
    double *dimension = NULL;
    double score,min_score;

    if(!first) dimension =  xmalloc(dimensions * sizeof(double));
    
    DEBUG("*** Starting categorizing\n");

    scale_score = 1;

    for(forest_idx = 0;forest_idx < forest_count;forest_idx++)  calculate_sample_score_range(forest_idx); // calculate score range for socre scaling

    while(fgets(input_line,INPUT_LEN_MAX,in_stream) != NULL) 
    {
        lines++;

        if(header && lines == 1) continue;

        value_count = parse_csv_line(values,DIM_MAX,input_line,input_separator);
        
        if(first && value_count) 
        {
            init_dims(value_count);
            dimension =  xmalloc(dimensions * sizeof(double));
            first = 0;
        }

        if(value_count)
        { 
            parse_values(dimension,values,value_count,0);

            if(aggregate)
            {
                forest_idx = find_forest(value_count,values,0);
                if(forest_idx > -1) 
                {
                    aggregate_values(forest_idx,dimension);
                    forest[forest_idx].total_rows++;
                }
            } else
            {
                best_forest_idx = -1;

                for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
                {
                    if(!forest[forest_idx].filter)
                    {
                            score = calculate_score(forest_idx,dimension);

                            if(best_forest_idx == -1 || score <= min_score)
                            {
                                min_score = score;
                                best_forest_idx = forest_idx;
                                forest[best_forest_idx].total_rows++;
                            }
                    }
                }

                if(best_forest_idx >= 0 && (!score_limit || (score_limit && min_score <= get_forest_score(best_forest_idx))))
                    print_(outs,min_score,lines,best_forest_idx,value_count,values,dimension,print_string,"rscldavxCtnem");
            }
        }
    }

    if(aggregate)
    {
        for(i = 0;i < forest_count;i++)
        {
            best_forest_idx = -1;
            if(forest[i].summary != NULL)
            {
                for(j = 0;j < forest_count;j++)
                {
                    if(!forest[j].filter)
                    {
                        score = calculate_score(j,forest[i].summary);
                        if(best_forest_idx == -1 || score < min_score)
                        {
                            min_score = score;
                            best_forest_idx = j;
                        }
                    }
                }
            }
            if(best_forest_idx >= 0 && (!score_limit || (score_limit && min_score <= get_forest_score(best_forest_idx))))
                print_(outs,min_score,0,best_forest_idx,0,NULL,forest[i].summary,print_string,"sdaxCtnem");
        }
    }

    scale_score = save_scale_score;

    if(dimension != NULL) free(dimension);
}


/* print category label of those forests which are not
 * filtered and which are not used in analysis
 *
 * This can be used to indicate which probably expected data is missing from analyzed file
 */
void print_missing_categories(FILE *outs,char *format)
{
    int i;

    for(i = 0;i < forest_count;i++)
    {
        if(!forest[i].filter && !forest[i].analyzed) print_(outs,0.0,0,i,0,NULL,NULL,format,"Catm");
    }
}

/* Reset sample count for a forest given by parameter
 */
void remove_samples(char *forest_string)
{
    int forest_idx;

    forest_idx = search_forest_hash(forest_string);

    if(forest_idx >= 0) 
    {
        forest[forest_idx].X_count = 0;
    } else
    {
        info("No forest having string",forest_string,NULL);
    }
}

