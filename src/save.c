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

static char *W_global = "G;%d;\"%s\";\"%s\";%d;%d;\"%s\";\"%c\";%d;%f;%f;\"%s\";\"%s\";%d;\"%s\";%d\n";
static char *W_forest = "F;\"%s\";%f;%d;%d\n";
static char *W_sample = "S;%s\n";

static char input_line[INPUT_LEN_MAX];

static 
void write_error()
{
     panic("Error while saving forest data to file",NULL,NULL);
}

/*
 * Save all global data to one line
 */
static
void write_global_data(FILE *w,int f_count)
{
    char *filter_string;

    filter_string = make_csv_line(cat_filter,cat_filter_count,';');

    if(fprintf(w,W_global,dimensions,label_dims ? label_dims : "",print_string ? print_string : "",tree_count,samples_max,category_dims ? category_dims : "",\
                input_separator,header,outlier_score,prange_extension_factor,ignore_dims ? ignore_dims : "",include_dims ? include_dims : "",f_count,filter_string,decimals) < 0)
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
        sprintf(f,"%.*f|",decimals,dim[i]);
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

    if(fprintf(w,W_forest,f->category ? f->category : "",f->c,f->heigth_limit,f->X_count) < 0) write_error();

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
 * S = sample
 */
void
write_forest_file(FILE *data_file)
{
    int i;
    
    write_global_data(data_file,forest_count);

    for(i = 0;i < forest_count;i++)
    {
        save_forest(i,data_file);
    }
}

/*
   Parse global parametes from forest file
   */
static 
void parse_G(char *l)
{
    int value_count,i,c;
    char *v[100];
    char *f[FILTER_MAX];

    value_count = parse_csv_line(v,100,l,';');

    if(value_count == 16) // change this too if parameter count changes
    {
        dimensions = atoi(v[1]);
        label_dims = xstrdup(v[2]);
        label_idx_count = parse_dims(v[2],label_idx);
        print_string = xstrdup(v[3]);
        tree_count = atoi(v[4]);
        samples_max = atoi(v[5]);
        category_dims = xstrdup(v[6]);
        category_idx_count = parse_dims(v[6],category_idx);
        input_separator = v[7][0];
        header = atoi(v[8]);
        outlier_score = atof(v[9]);
        prange_extension_factor = atof(v[10]);
        ignore_dims = xstrdup(v[11]);
        ignore_idx_count = parse_dims(v[11],ignore_idx);
        include_dims = xstrdup(v[12]);
        include_idx_count = parse_dims(v[12],include_idx);
        forest_count = atoi(v[13]);

        c = parse_csv_line(f,FILTER_MAX,v[14],';');
        for(i = 0;i < c;i++) add_category_filter(f[i]);

        decimals = atoi(v[15]);

        samples_total = tree_count * samples_max;
    }
}

int parse_F(int forest_idx,char *l)
{
    int value_count;
    char *v[100];

    struct forest *f = &forest[forest_idx];

    value_count = parse_csv_line(v,100,l,';');

    if(value_count == 5)
    {
        f->category = xstrdup(v[1]);
        f->c = atof(v[2]);
        f->heigth_limit = atof(v[3]);
        f->X = NULL;
        f->X_count = 0;
        f->X_cap = 0;
        f->t = xmalloc(tree_count * sizeof(struct tree));
        f->min = NULL;
        f->max = NULL;
        f->dim_density = xmalloc(dimensions * sizeof(double));

        f->X_cap = atoi(v[4]);
        f->X = xmalloc(f->X_cap * sizeof(struct sample));

        add_forest_hash(forest_idx,f->category);

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
        } else
        {
            return 0;
        }
    } while(input_line[0] != 'G');

    if(!dimensions || !tree_count) return 0;

    forest_cap = forest_count;

    if(!forest_count) return 1;

    forest = xmalloc(forest_cap * sizeof(struct forest));

    f_count = 0;
        
    if(fgets(input_line,INPUT_LEN_MAX,data_file) == NULL) return 0;

    do
    {
        if(input_line[0] == 'F')
        {
            if(!parse_F(f_count,input_line)) return 0;
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
            f_count++;
        } else
        {
            return 0;
        }
    } while(f_count < forest_count);

    return 1;
}

