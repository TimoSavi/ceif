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


/* Populates values from parsed input row
 */
static
void populate_dimension(double *d,char **values,int value_count)
{
    int i;

    for(i = 0;i < dimensions;i++)
    {
        if(dim_idx[i] < value_count) d[i] = parse_dim_attribute(values[dim_idx[i]]);
    }
}

/* search through nodes, 
 * returns the length from last node
 */
static 
double search_last_node(int this_idx,struct node *n,double *dimension,int heigth)
{
    struct node *this = &n[this_idx];

    if(this->left == -1 && this->rigth == -1)
    {
        return (double) heigth + c(this->sample_count);
    }

    if(dot(dimension,this->n) < this->pdotn)
    {
        if(this->left == -1) return (double) heigth;
        return search_last_node(this->left,n,dimension,heigth + 1);
    } else
    {
        if(this->rigth == -1) return (double) heigth;
        return search_last_node(this->rigth,n,dimension,heigth + 1);
    }
}
    

/*
 * travel through a tree with a test case.
 * return the path length
 */
static 
double calculate_path_length(struct tree *t,double *dimension)
{
    return search_last_node(t->first,t->n,dimension,0);
}


/* calculates score for a data row in given forest
 * returns score
 */
double calculate_score(int forest_idx,double *dimension)
{
    int i;
    struct forest *f = &forest[forest_idx];
    double path_length = 0;

    if(f->t != NULL)
    {
        for(i = 0;i < tree_count;i++)
        {
            path_length += calculate_path_length(&f->t[i],dimension);
        }
    }
    
    path_length = path_length / tree_count;  // turn to average 
 
    return (1/pow(2,path_length/f->c));
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
    static char l[10240];

    l[0] = '\000';

    for(i = 0;i < label_idx_count;i++)
    {
        if(label_idx[i] < value_count)
        {
            strcat(l,values[label_idx[i]]);
            if(i < label_idx_count - 1) strcat(l,LABEL_SEPARATOR);
        }
    }
    return l;
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
               case 's':
                   fprintf(outs,"%f",score);
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
                       *printf_format != '\000' ? fprintf(outs,printf_format,dimension[i]) : fprintf(outs,float_format,decimals,dimension[i]);
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
                   fprintf(outs,"%s",CATEGORY_SEPARATOR);
                   break;
               case '.':
                   fprintf(outs,"%s",LABEL_SEPARATOR);
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
    int i;

    dimensions = 0;

    for(i = 0;i < value_count;i++) 
    {
        if((!check_idx(i,ignore_idx_count,ignore_idx) || check_idx(i,include_idx_count,include_idx)) &&
            !check_idx(i,category_idx_count,category_idx) && 
            !check_idx(i,label_idx_count,label_idx))
        {
            dim_idx[dimensions] = i;
            dimensions++;
        }
    }
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


/* analyze data from file. 
 * All lines are analyzed against loaded forest/tree data
 * and anomalies (having score > outlier_score) using printing mask
 */
void
analyze(FILE *in_stream, FILE *outs,char *not_found_format)
{
    int i;
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
            populate_dimension(dimension,values,value_count);

            forest_idx = find_forest(value_count,values,1);

            if(forest_idx >= 0)
            {
                if(aggregate)
                {
                    aggregate_values(forest_idx,dimension);
                } else
                {
                    score = calculate_score(forest_idx,dimension);
                    if(score >= outlier_score) 
                    {
                        print_(outs,score,lines,forest_idx,value_count,values,dimension,print_string,"rscldavxCt");
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
        for(i = 0;i < forest_count;i++)
        {
            if(forest[i].summary != NULL && !forest[i].filter)
            {
                score = calculate_score(i,forest[i].summary);
                if(score >= outlier_score) 
                {
                    print_(outs,score,0,i,0,NULL,forest[i].summary,print_string,"rsdaxCt");
                }
            }
        }
    }

    if(dimension != NULL) free(dimension);
}


/* Categorize dimensions
 * All lines are analyzed against loaded forest/tree data
 * All forests are analyzed and a forest having lowest anomaly score is selected as category forest
 */
void
categorize(FILE *in_stream, FILE *outs)
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
            populate_dimension(dimension,values,value_count);

            if(aggregate)
            {
                forest_idx = find_forest(value_count,values,0);
                if(forest_idx > -1) aggregate_values(forest_idx,dimension);
            } else
            {
                best_forest_idx = -1;

                for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
                {
                    if(!forest[forest_idx].filter)
                    {
                            score = calculate_score(forest_idx,dimension);
                            if(best_forest_idx == -1 || score < min_score)
                            {
                                min_score = score;
                                best_forest_idx = forest_idx;
                            }
                    }
                }

                if(best_forest_idx >= 0) print_(outs,min_score,lines,best_forest_idx,value_count,values,dimension,print_string,"rscldavxCt");
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
            if(best_forest_idx >= 0) print_(outs,min_score,0,best_forest_idx,0,NULL,forest[i].summary,print_string,"sdaxCt");
        }
    }

    if(dimension != NULL) free(dimension);
}


/* print category label of thos forests which are not
 * filtered and which ar enot used in analysis
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
