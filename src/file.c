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
    _2P("N-vector adjust is %s\n",_O(n_vector_adjust));
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

