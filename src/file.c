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
#include <locale.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#define VALUES_MAX (2*DIM_MAX)             // Max values for parsing a csv line

/* make a csv line from strings 
   return pointer to string
 */
char *
make_csv_line(char *values[],int value_count,char separator)
{
    static char line[INPUT_LEN_MAX];
    char s[2];
    int i;

    line[0] = '\000';
    s[0] = separator;
    s[1] = '\000';

    for(i = 0;i < value_count;i++)
    {
        strcat(line,values[i]);
        if(i < value_count - 1) strcat(line,s);
    }

    return line;
}

/* parse a csv line.
   return number of values
   pointers to each value will be written to values
   read only max values, rest are ignored
   write NULL for each separator, line will be modified
    
   values can be enclosed to double quotes (")
 */
int
parse_csv_line(char *values[],int max_values,char *line,char separator)
{
    int n = 0;
    int in_quote;
    char *new_line;
    char *start = line;

    if(line == NULL || max_values < 1) return 0;

    while(*line != '\000' && n < max_values)
    {
        if(*line == '"') 
        {
            in_quote = 1;
            line++;
        } else
        {
            in_quote = 0;
        }

        values[n++] = line;
        while(*line != '\000' && ((!in_quote && (*line != separator || (*line == separator && line > start && line[-1] == '\\'))) || (in_quote && *line != '"'))) line++;
        if (in_quote && *line == '"') *line++ = '\000';
        while(*line != '\000' && *line != separator) line++;
        if (*line == separator) *line++ = '\000';
    }

    if(n) {
        new_line = strchr(values[n - 1],'\n');
        if(new_line != NULL) *new_line = '\000';
    }

    return n;
}

static void 
check_dim_range(int index)
{
    char n[100];
    
    if(index < 1 || index > DIM_MAX)
    {
        sprintf(n,"Valid dimension numbers are 1 - %i",DIM_MAX); 
        panic(n,NULL,NULL);
    }
}

/* parse  dims list 
   list is csv list of dim numbers (starts from 1)
   dims will be samed to table as table indexes (minus 1)
   range can be expressed as first-last.

   e.g.

   1,2,10-15

   returns the number if dims

   */
int
parse_dims(char *optarg, int *array)
{
    char *value[DIM_MAX];
    char *range[2];
    int i,j;
    int ranges = 0;
    int values = 0;
    int dims = 0;

    values = parse_csv_line(value,DIM_MAX,optarg,',');

    for(i = 0;i < values;i++)
    {
        ranges = parse_csv_line(range,2,value[i],'-');
        if(ranges > 1)
        {
            for(j = atoi(range[0]);j <= atoi(range[1]);j++)
            {
                check_dim_range(j);
                array[dims++] = j - 1;
            }
        } else 
        {
                j = atoi(value[i]);
                check_dim_range(j);
                array[dims++] = j - 1;
        }
    }
    return dims;
}

#define _I fprintf(outs,"  ")
#define _P(...) fprintf(outs, __VA_ARGS__)
#define _2P(...) _I; _P( __VA_ARGS__)
#define _3P(...) _I; _2P( __VA_ARGS__)
#define _O(b) (b ? "On" : "Off")
#define _S(i,max) if(i < max - 1) _P("%c",',')
#define _D(name,count,array) _2P(name); for(i = 0;i < count;i++) {_P("%d",array[i]+1);_S(i,count);} _P("\n")

/* Print forest contents in user readbaly form
 */
void
print_forest_info(FILE *outs)
{
    int forest_idx,i;
    struct forest *f;
    char outstr[100];
    struct tm *tmp;
        
    _P("Global setting:\n");
    _2P("Number of forests: %d\n",forest_count);
    _2P("Number of analyzed dimensions: %d\n",dimensions);
    _2P("Number of samples/tree: %d\n",samples_max);
    _2P("Number of trees: %d\n",tree_count);
    _2P("Number of decimals: %d\n",decimals);
    _2P("P-range extension factor: %f\n",prange_extension_factor);
    _2P("Outlier score: %f\n",outlier_score);
    _2P("Input separator: %c\n",input_separator);
    _2P("Output separator: %c\n",list_separator);
    _2P("Header is %s\n",_O(header));
    _2P("Auto weigth is %s\n",_O(auto_weigth));
    _2P("Aggregate is %s\n",_O(aggregate));
    _2P("Unigue samples is %s\n",_O(unique_samples));
    _2P("Print string: \"%s\"\n",print_string);


    _D("Dimensions used in analysis: ",dimensions,dim_idx);
    _D("User ignored dimensions: ",ignore_idx_count,ignore_idx);
    _D("User included dimensions: ",include_idx_count,include_idx);
    _D("Category dimensions: ",category_idx_count,category_idx);
    _D("Label dimensions: ",label_idx_count,label_idx);
    _D("Dimensions treated as text: ",text_idx_count,text_idx);


    if(forest_count) _P("\nForest data:\n");
    for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
    {
        f = &forest[forest_idx];
        _2P("\nForest category string: %s\n",f->category);
        _3P("Filter is %s\n",_O(f->filter));
        _3P("Number of samples: %d\n",f->X_count);
        _3P("Average path length (c): %f\n",f->c);
        _3P("Max. tree heigth: %d\n",f->heigth_limit);

        tmp = localtime(&f->last_updated);
        if(tmp != NULL)
        {
            outstr[0] = '\000';
            strftime(outstr, sizeof(outstr), "%c", tmp);
            _3P("Last updated: %s\n",outstr);
        }

        _P("\n");
        _3P("%15s%*s\n","Dimension sample value summary:",12*dimensions,"Dimension");

        _3P("%15s","");
        for(i = 0;i < dimensions;i++) _P("%25d",i+1);
        _P("\n");

        _3P("%15s","Maximum value");
        for(i = 0;i < dimensions;i++) _P("%25.*f",decimals,f->max[i]);
        _P("\n");
        
        _3P("%15s","Minimum value");
        for(i = 0;i < dimensions;i++) _P("%25.*f",decimals,f->min[i]);
        _P("\n");
        
        _3P("%15s","Average value");
        for(i = 0;i < dimensions;i++) _P("%25.*f",decimals,f->avg[i] / (f->filter ? f->X_count : 1));
        _P("\n");
    }
}

#define DENSITY_MAX 100

/* Print sample ascii density map
 * For each forest the density of each dimensions is preinted as ascii char map,
 * There are 100 buckets between smallest and largest value of all dimension attributes
 * For each buck the ascii char is print according the number of samples in that buck (75%, 50%, 50%, 25% and 0%) (#, =, -, . and space)
 */
void
print_sample_density(FILE *outs,int common_scale)
{
    int forest_idx,i,j,first;
    static char *digit="0123456789#";
    double min = 0,max = 0,bucket_size = 0,percentage;
    struct forest *f;
    int density[DENSITY_MAX];
    size_t bucket;

    _P("Sample value density map\n");
    _P("Each dimensions is divided into %d buckets, the digit under a bucket means number of 1/10 of samples in that bucket, # means all samples belong to one bucket\n\n",DENSITY_MAX);
    _P("Empty means no samples\n\n");

    if(!forest_count) return;

    if(common_scale)
    {
        first = 1;
        for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
        {
            f = &forest[forest_idx];

            if(!f->X_count) continue;

            // find all dimensions min and max values
            if(first)
            {
                min = f->X[0].dimension[0];
                max = f->X[0].dimension[0];
                first = 0;
            }

            for(i = 0;i < f->X_count;i++)
            {
                for(j = 0;j < dimensions;j++)
                {
                    min = fmin(min,f->X[i].dimension[j]);
                    max = fmax(max,f->X[i].dimension[j]);
                }
            }
        }
        if(first) return;
        bucket_size = (max - min) / (double) DENSITY_MAX;
    }

    _P("\n### Density by forest ###\n");
 
    for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
    {
        f = &forest[forest_idx];

        if(!f->X_count) continue;

        if(!common_scale)
        {
            // find all dimensions min and max values
            min = f->X[0].dimension[0];
            max = f->X[0].dimension[0];

            for(i = 0;i < f->X_count;i++)
            {
                for(j = 0;j < dimensions;j++)
                {
                    min = fmin(min,f->X[i].dimension[j]);
                    max = fmax(max,f->X[i].dimension[j]);
                }
            }
            bucket_size = (max - min) / (double) DENSITY_MAX;
        }

        if(bucket_size == 0.0) continue;

        _P("\nForest category string: %s, bucket size %.*f\n\n",f->category,decimals,bucket_size);

        _2P("%10s","Min...Max");
        _P("%15.*f ",decimals,min);
        for(i = 0;i < DENSITY_MAX;i++) _P("-");
        _P("%15.*f\n\n",decimals,max);

        for(j = 0;j < dimensions;j++)
        {
            for(i = 0;i < DENSITY_MAX;i++) density[i] = 0;

            _2P("%10s","Dimension");
            _P("%15d ",j+1);

            for(i = 0;i < f->X_count;i++)
            {
                bucket = (size_t) ((f->X[i].dimension[j] - min) / bucket_size);
                if(bucket >= DENSITY_MAX) bucket = DENSITY_MAX - 1;
                density[bucket]++;
            }

            for(i = 0;i < DENSITY_MAX;i++) 
            {
                percentage = ((double) density[i] / f->X_count) * 10.0;
                percentage > 0.0 ? _P("%c",digit[(size_t) percentage]) : _P("%c",' ');
            }
            _P("\n");
        }
    }

    _P("\n### Density by dimension ###\n");

    for(j = 0;j < dimensions;j++)
    {
        if(!common_scale)
        {
           f = &forest[0];
            min = f->X[0].dimension[j];
            max = f->X[0].dimension[j];

            for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
            {
                f = &forest[forest_idx];

                for(i = 0;i < f->X_count;i++)
                {
                    min = fmin(min,f->X[i].dimension[j]);
                    max = fmax(max,f->X[i].dimension[j]);
                }
                bucket_size = (max - min) / (double) DENSITY_MAX;
            }
        }
        
        _P("\nDimension %d, bucket size %.*f\n\n",j + 1,decimals,bucket_size);

        _2P("%10s","Min...Max");
        _P("%35.*f ",decimals,min);
        for(i = 0;i < DENSITY_MAX;i++) _P("-");
        _P("%15.*f\n\n",decimals,max);

        for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
        {
            f = &forest[forest_idx];

            _2P("%10s%35s ","Category:",f->category);

            for(i = 0;i < DENSITY_MAX;i++) density[i] = 0;

            for(i = 0;i < f->X_count;i++)
            {
                bucket = (size_t) ((f->X[i].dimension[j] - min) / bucket_size);
                if(bucket >= DENSITY_MAX) bucket = DENSITY_MAX - 1;
                density[bucket]++;
            }

            for(i = 0;i < DENSITY_MAX;i++) 
            {
                percentage = ((double) density[i] / f->X_count) * 10.0;
                percentage > 0.0 ? _P("%c",digit[(size_t) percentage]) : _P("%c",' ');
            }
            _P("\n");
        }
    }
}
