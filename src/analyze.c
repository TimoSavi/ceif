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

#define INPUT_LEN_MAX 1048576

static char input_line[INPUT_LEN_MAX];
static int first = 1;


/* find a forest for a data row
 * 
 * returns index to forest table, -1 if not found
 */
static 
int find_forest(int value_count,char **values)
{
    int i;
    char *category_string;


    category_string = make_category_string(value_count,values);

    for(i = 0;i < forest_count;i++)
    {
        if (!forest[i].filter && strcmp(category_string,forest[i].category) == 0) return i;
    }

    return -1;
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
static 
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
 */
static 
unsigned int score_to_rgb(double score)
{
    unsigned int red,blue,green;

    green = score < 0.5 ? (unsigned int) (2 * 0xff * score) : (unsigned int) (0xff * 2 * (1 - score)); 
    red = score > 0.5 ? (unsigned int) (2 * 0xff * (score - 0.5)) : 0; 
    blue = score < 0.5 ? (unsigned int) (2 * 0xff * (0.5 - score)) : 0; 
    
    return (red << 16) + (green << 8) + blue;
}                      
        
/* Print outlier info
 * printing is done using printf style string and %-directives
 */
static 
void print_outlier(FILE *outs, double score, int lines,int forest_idx,int value_count,char **values,double *dimension)
{
    int i;
    char *c = print_string;

    while(*c != '\000')
    {
       if(*c == '%' && c[1] != '\000')
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
                   if(label_dim >= 0 && label_dim < value_count) fprintf(outs,"%s",values[label_dim]);
                   break;
               case 'd':
                   for(i = 0;i < dimensions;i++) i < dimensions - 1 ? fprintf(outs,"%f,",dimension[i]) : fprintf(outs,"%f",dimension[i]);
                   break;
               case 'v':
                   for(i = 0;i < value_count;i++) i < value_count - 1 ? fprintf(outs,"%s%c",values[i],input_separator) : fprintf(outs,"%s",values[i]);
                   break;
               case 'x':
                   fprintf(outs,"%06X",score_to_rgb(score));
                   break;
               case 'C':
                   fprintf(outs,"%s",forest[forest_idx].category);
                   break;
               case '%':
                   fprintf(outs,"%c",'%');
                   break;
           }
       } else
       {
           fprintf(outs,"%c",*c);
       }
       c++;
    }
    if(*print_string) fprintf(outs,"%c",'\n');
}


/*  mark category and label dims as non dimensions dims  and populate dim_idx and
 *  category_idx tables
 *  */
void
init_dims(int value_count)
{
    int i;

    dimensions = 0;
    categories = 0;

    if(label_dim >= 0) dims_ignore[label_dim] = 1;            // do not read as normal dimemnsion

    for(i = 0;i < value_count;i++) 
    {
        if(!dims_ignore[i] && !dims_category[i])
        {
            dim_idx[dimensions] = i;
            dimensions++;
        }

        if(dims_category[i]) 
        {
            category_idx[categories] = i;
            categories++;
        }
    }
}

/* analyze data from file. 
 * All lines are analyzed against loaded forest/tree data
 * and anomalies (having score > outlier_score) using printing mask
 */
void
analyze(FILE *in_stream, FILE *outs)
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
             forest_idx = find_forest(value_count,values);
             if(forest_idx >= 0)
             {
                  populate_dimension(dimension,values,value_count);
                  score = calculate_score(forest_idx,dimension);
                  if(score >= outlier_score) 
                  {
                      print_outlier(outs,score,lines,forest_idx,value_count,values,dimension);
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
    int value_count;
    int lines = 0;
    int forest_idx;
    int best_forest_idx;
    char *values[DIM_MAX];
    double *dimension;
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

        best_forest_idx = -1;

        if(value_count)
        { 
             populate_dimension(dimension,values,value_count);

             for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
             {
                 if(!forest[forest_idx].filter)
                 {
                     min_score = calculate_score(forest_idx,dimension);
                     best_forest_idx = forest_idx;
                     forest_idx = forest_count;
                 }
             }

             if(best_forest_idx >= 0)
             {
                 forest_idx = best_forest_idx + 1;
                 for(;forest_idx < forest_count;forest_idx++)
                 {
                     if(!forest[forest_idx].filter)
                     {
                         score = calculate_score(forest_idx,dimension);
                         if(score < min_score)
                         {
                             min_score = score;
                             best_forest_idx = forest_idx;
                         }
                     }
                 }
             }

             if(best_forest_idx >= 0) print_outlier(outs,min_score,lines,best_forest_idx,value_count,values,dimension);
        }
    }
    if(dimension != NULL) free(dimension);
}


