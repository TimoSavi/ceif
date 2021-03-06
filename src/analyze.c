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

   relative distance is 1 if the actual distance is the same as forersst average sample distance
   relative distance < 1 if the actual distance is smaller as forerst average sample distance but never larger than MIN_REL_DIST
   relative distance > 1 if the actual distance is larger as forersst average sample distance

   Scale the a if auto scaling is in use
 */
#define MIN_REL_DIST 0.05
double nearest_rel_distance(double *a, int sample_count,struct sample *samples,struct forest *f)
{
    int i;
    double *dim;
    double distance,d; 

    if(auto_weigth)
    {
        dim = scale_dimension(a,f);
    } else
    {
        dim = a;
    }

    distance = v_dist_nosqrt(dim,samples[0].dimension);

    for(i = 1;i < sample_count;i++)
    {
        d = v_dist_nosqrt(dim,samples[i].dimension);

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
        if(nearest && f->avg_sample_dist > 0.0)
        {
            return (double) heigth + c((double) this->sample_count / nearest_rel_distance(dimension,this->sample_count,this->samples,f));
        }
        return (double) heigth + c(this->sample_count);
    }

    if(wdot(dimension,this->n,f->scale_range_idx,f->min,f->max) < this->pdotn)
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
    double path_length = 0;

    if(f->t != NULL)
    {
        for(i = 0;i < tree_count;i++)
        {
            path_length += calculate_path_length(f,&f->t[i],dimension);
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
    int l,save_auto_weigth = auto_weigth;

    auto_weigth = 0;  // no need for scaling

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

    auto_weigth = save_auto_weigth;

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
    if(scale_score) return calculate_score_scale(forest_idx,dimension);

    return _score(forest_idx,dimension);
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


/* Print something
 * printing is done using printf style string and %-directives
 */
void print_(FILE *outs, double score, int lines,int forest_idx,int value_count,char **values,double *dimension,char *format,char *directives)
{
    int i;
    char *c = format;
    char outstr[100];
    struct tm *tmp;

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
                case'n':
                    fprintf(outs,"%d",forest[forest_idx].total_rows);
                    break;
                case'o':
                    fprintf(outs,"%d",forest[forest_idx].analyzed_rows);
                    break;
                case'h':
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
                        if(i < dimensions - 1) fputc(list_separator,outs);
                    }
                   break;
               case 'a':
                   for(i = 0;i < dimensions;i++)
                   {
                       *printf_format != '\000' ? fprintf(outs,printf_format,forest[forest_idx].avg[i]) : fprintf(outs,float_format,decimals,forest[forest_idx].avg[i]);
                       if(i < dimensions - 1) fputc(list_separator,outs);
                   }
                   break;
               case 'v':
                   for(i = 0;i < value_count;i++) i < value_count - 1 ? fprintf(outs,"%s%c",values[i],list_separator) : fprintf(outs,"%s",values[i]);
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
       } else
       {
           if(*c == '\\' && c[1] != '\000')
           {
               c++;
               switch(*c)
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
               }
           } else
           {
               fprintf(outs,"%c",*c);
           }
       }
       c++;
    }
    if(*format) fprintf(outs,"%c",'\n');
}


/* check if index is in index table.
 */
int check_idx(int index,int index_count, int *idx_table)
{
    int i;

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

/* Init auto (max) score for a forest 
 * Auto score is initialized by the maximun sample score value.
 * Sample values are sligthly randomly moved using v_expand, 
 * */
void calculate_forest_auto_score(int forest_idx)
{
    double score;
    int i;
    struct forest *f;

    f = &forest[forest_idx];
        
    f->auto_score = 0.0;

    if(f->filter) return;

    for(i = 0;i < f->X_count;i++)
    {
        score = _score(forest_idx,v_expand(f->X[i].dimension,f->avg,f->dim_density,f->X_count));
        if(score > f->auto_score) f->auto_score = score;
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

    for(i = 0;i < f->X_count;i++) all_scores[i] = _score(forest_idx,f->X[i].dimension);

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
            score = calculate_score(forest_idx,f->X[i].dimension);
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

/* calculate average sample score for a forest
 * Average is counted once, this is chekked using f->average_score (it is practically never zero)
 *
 */
void calculate_average_sample_score(int forest_idx)
{
    int i;
    struct forest *f;
    double stddev = 0.0;
    double *scores;

    f = &forest[forest_idx];

    if(f->filter) return;
     
    f->average_score = 0.0;
    f->min_score = 1.0;

    scores = xmalloc(f->X_count * sizeof(double));

    for(i = 0;i < f->X_count;i++)
    {
        scores[i] = _score(forest_idx,f->X[i].dimension);

        if(scores[i] < f->min_score) f->min_score = scores[i];

        f->average_score += scores[i];
    }

    f->average_score /= (double) f->X_count;

    for(i = 0;i < f->X_count;i++)
    {
        stddev += POW2(f->average_score - scores[i]);
    }

    stddev = sqrt(stddev / (double) f->X_count);

    f->average_score += stddev * average_score_factor;

    f->max_score = calculate_max_score(forest_idx) + stddev / 2.0;   // add stddev/2 in order to make sure that scaled score remains always < 1.0

    if(f->max_score > 1.0) f->max_score = 1.0;

    free(scores);
}

/* Return score for a forest, special scores for auto and average socres are handled here
 */
inline 
double get_forest_score(int forest_idx)
{
    if(outlier_score == AUTO_SCORE) return forest[forest_idx].auto_score;
    if(outlier_score == AVERAGE_SCORE) return forest[forest_idx].average_score;
    if(percentage_score) return forest[forest_idx].percentage_score;
    return outlier_score;
}

/* Calculate appropriate automatic score for a forest
 */
inline 
void calculate_forest_score(int forest_idx)
{
    if(forest[forest_idx].average_score == 0.0 && (outlier_score == AVERAGE_SCORE || scale_score)) 
    {
        calculate_average_sample_score(forest_idx);
    } else if(forest[forest_idx].auto_score == 0.0 && outlier_score == AUTO_SCORE)
    {
        calculate_forest_auto_score(forest_idx);
    } if(forest[forest_idx].percentage_score == 0.0 && percentage_score)
    {
        calculate_forest_percentage_score(forest_idx);
    }
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
    double score;
    
    if(!first) dimension =  xmalloc(dimensions * sizeof(double));
        
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
                calculate_forest_score(forest_idx);

                forest[forest_idx].total_rows++;

                if(aggregate)
                {
                    aggregate_values(forest_idx,dimension);
                } else
                {
                    if(take_this_row(forest[forest_idx].total_rows))   // check if analyzed rows are reservoir sampled
                    {
                        forest[forest_idx].analyzed_rows++;

                        score = calculate_score(forest_idx,dimension);

                        if(average_format != NULL) forest[forest_idx].test_average_score += score;

                        if(score > get_forest_score(forest_idx))
                        {
                            forest[forest_idx].high_analyzed_rows++;
                            print_(outs,score,lines,forest_idx,value_count,values,dimension,print_string,"rscldavxCtnoh");
                        }
                    }
                }
            } else
            {
                if(not_found_format != NULL && find_forest(value_count,values,0) == -1) print_(outs,0,lines,0,value_count,values,dimension,not_found_format,"dvcl");
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
                
                if(score > get_forest_score(forest_idx))
                {
                    forest[forest_idx].high_analyzed_rows++;
                    print_(outs,score,0,forest_idx,0,NULL,forest[forest_idx].summary,print_string,"rsdaxCtnoh");
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
 */
void
categorize(FILE *in_stream, int score_limit, FILE *outs)
{
    int i,j;
    int value_count;
    int lines = 0;
    int forest_idx;
    int best_forest_idx;
    char *values[DIM_MAX];
    double *dimension = NULL;
    double score,min_score;

    if(!first) dimension =  xmalloc(dimensions * sizeof(double));

    for(forest_idx = 0;forest_idx < forest_count;forest_idx++)  /* Get forest sample score min ... max range */
    {
        if(!forest[forest_idx].filter)
        {
            if(forest[forest_idx].average_score == 0.0) calculate_average_sample_score(forest_idx);
        }
    }

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
                            score = calculate_score_scale(forest_idx,dimension);

                            if(best_forest_idx == -1 || score <= min_score)
                            {
                                min_score = score;
                                best_forest_idx = forest_idx;
                                forest[best_forest_idx].total_rows++;
                            }
                    }
                }

                if(best_forest_idx >= 0 && (!score_limit || (score_limit && min_score <= get_forest_score(best_forest_idx))))
                    print_(outs,min_score,lines,best_forest_idx,value_count,values,dimension,print_string,"rscldavxCtn");
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
                        score = calculate_score_scale(j,forest[i].summary);
                        if(best_forest_idx == -1 || score < min_score)
                        {
                            min_score = score;
                            best_forest_idx = j;
                        }
                    }
                }
            }
            if(best_forest_idx >= 0 && (!score_limit || (score_limit && min_score <= get_forest_score(best_forest_idx))))
                print_(outs,min_score,0,best_forest_idx,0,NULL,forest[i].summary,print_string,"sdaxCtn");
        }
    }

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
        if(!forest[i].filter && !forest[i].analyzed) print_(outs,0.0,0,i,0,NULL,NULL,format,"Cat");
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

