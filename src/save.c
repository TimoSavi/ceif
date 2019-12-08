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
 *    along with ffe; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *    F607480034
 *    HJ9004-2
 *
 */
#include "ceif.h"

/* formats for write and reading data
 */

static char *W_global = "G;%d;%d;\"%s\";%d;%d;\"%s\";\"%c\";%d;%f;%f;\"%s\";\"%s\";%d\n";
static char *W_forest = "F;\"%s\";%f;%d\n";
static char *W_sample = "S;%s\n";

static char input_line[INPUT_LEN_MAX];

static 
void write_error()
{
     panic("Error while saving forest data to file",NULL,NULL);
}

static 
int valid_forest(struct forest *f)
{
    return (f->X == NULL || f->X_count >= SAMPLES_MIN);  /*  dont save forests with no decent amount of samples, f->X is null if forest is read from file */
}

/*
 * Save all global data to one line
 */
static
void write_global_data(FILE *w,int f_count)
{
    if(fprintf(w,W_global,dimensions,label_dim,print_string ? print_string : "",tree_count,samples_max,category_dims ? category_dims : "",\
                input_separator,header,outlier_score,prange_extension_factor,ignore_dims ? ignore_dims : "",include_dims ? include_dims : "",f_count) < 0)
    {
        write_error();
    }
}

/* write dimension data to csv string
   */
static 
char *dim_to_csv(int size,double *dim)
{
    static char csv[DIM_MAX*21];
    char f[200];
    int i;

    csv[0] = '\000';

    for(i = 0;i < size;i++) 
    {
        sprintf(f,"%f|",dim[i]);
        strcat(csv,f);
    }
    if(size) csv[strlen(csv) - 1] = '\000';
    return csv;
}


/*
 * save data for a forest
 */
static 
void
save_forest(int forest_idx,FILE *w)
{
    int i;
    struct forest *f = &forest[forest_idx];

    if(fprintf(w,W_forest,f->category ? f->category : "",f->c,f->heigth_limit) < 0) write_error();

    for(i = 0;i < f->X_count;i++)
    {
        if(fprintf(w,W_sample,dim_to_csv(dimensions,f->X[i].dimension)) < 0) write_error();
    }
}


/*
 * Save forest data to csv file, fields are sperated by semicolon, dimensions are seprated by pipe (|)
 * data file format:
 * ID;data1;data2;...
 * where ID:
 * G = global data
 * F = forest data 
 * T = tree data for last detected forest
 * N = node data for laste detected tree
 */
void
write_forest_file(FILE *data_file)
{
    int i;
    int forests_to_save = forest_count;
    
    for(i = 0;i < forest_count;i++) if(!valid_forest(&forest[i])) forests_to_save--;

    write_global_data(data_file,forests_to_save);

    for(i = 0;i < forest_count;i++)
    {
        if(valid_forest(&forest[i])) save_forest(i,data_file);
    }
}

/*
   Parse global parametes from forest file
   */
static 
void parse_G(char *l)
{
    int value_count;
    char *v[100];

    value_count = parse_csv_line(v,100,l,';');

    if(value_count == 14) // change this too if parameter count changes
    {
        dimensions = atoi(v[1]);
        label_dim = atoi(v[2]);
        print_string = xstrdup(v[3]);
        tree_count = atoi(v[4]);
        samples_max = atoi(v[5]);
        category_dims = xstrdup(v[6]);
        parse_dims(v[6],1,dims_category);
        input_separator = v[7][0];
        header = atoi(v[8]);
        outlier_score = atof(v[9]);
        prange_extension_factor = atof(v[10]);
        ignore_dims = xstrdup(v[11]);
        parse_dims(v[11],1,dims_ignore);
        include_dims = xstrdup(v[12]);
        parse_dims(v[12],0,dims_ignore);
        forest_count = atoi(v[13]);

        samples_total = tree_count * samples_max;
    }
}

int parse_F(int forest_idx,char *l)
{
    int value_count;
    char *v[100];

    struct forest *f = &forest[forest_idx];

    value_count = parse_csv_line(v,100,l,';');

    if(value_count == 4)
    {
        f->category = xstrdup(v[1]);
        f->c = atof(v[2]);
        f->heigth_limit = atof(v[3]);
        f->X = NULL;
        f->X_count = 0;
        f->X_cap = 0;
        f->t = xmalloc(tree_count * sizeof(struct tree));
        f->min = xmalloc(dimensions * sizeof(double));
        f->max = xmalloc(dimensions * sizeof(double));
        f->dim_density = xmalloc(dimensions * sizeof(double));
        return 1;
    }
    return 0;
}


/*
 * read saved forest structure to  memory
 * returns 1 in case read was ok, 0 other wise
 */
int 
read_forest_file(FILE *data_file)
{
    int f_count,line,value_count;
    char *values[DIM_MAX];
    

    do
    {
        if(fgets(input_line,INPUT_LEN_MAX,data_file) != NULL)
        {
            if(input_line[0] == 'G')
            {
                parse_G(input_line);
            }
        }
    } while(input_line[0] != 'G');

    if(!dimensions || !tree_count) return 0;

    forest_cap = forest_count;
    forest = xmalloc(forest_cap * sizeof(struct forest));
    f_count = 0;
        
    if(fgets(input_line,INPUT_LEN_MAX,data_file) == NULL) return 0;

    do
    {
        if(input_line[0] == 'F')
        {
            if(input_line[0] != 'F' || !parse_F(f_count,input_line)) return 0;
            line = 0;
            do
            {
                if(fgets(input_line,INPUT_LEN_MAX,data_file) == NULL) return 1;
                if(input_line[0] == 'S')
                {
                    line++;
                    value_count = parse_csv_line(values,dimensions,&input_line[2],'|');
                    add_to_X(&forest[f_count],values,value_count,line,1);
                }
            } while(input_line[0] == 'S');
        }
        f_count++;
    } while(f_count < forest_count);

    return 1;
}

