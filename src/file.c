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
        while(*line != '\000' && ((!in_quote && *line != separator) || (in_quote && *line != '"'))) line++;
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

   */
void
parse_dims(char *optarg, int val, int *array)
{
    char *value[DIM_MAX];
    char *range[2];
    int i,j;
    int ranges = 0;
    int values = 0;

    values = parse_csv_line(value,DIM_MAX,optarg,',');

    for(i = 0;i < values;i++)
    {
        ranges = parse_csv_line(range,2,value[i],'-');
        if(ranges > 1)
        {
            for(j = atoi(range[0]);j <= atoi(range[1]);j++)
            {
                check_dim_range(j);
                array[j - 1] = val;
            }
        } else 
        {
                j = atoi(value[i]);
                check_dim_range(j);
                array[j - 1] = val;
        }
    }
}
