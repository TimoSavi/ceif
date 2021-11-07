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
#include <ctype.h>

#define VALUES_MAX (2*DIM_MAX)             // Max values for parsing a csv line

/* Make separated string. String is initialized with NULL item
 * Every call adds single item to string
 * string is returned for every call
 */
char *
make_separated_string(char *item, char separator)
{
    static char string[10240];
    static size_t write_pos = 0;
    char *p;

    if(item == NULL)
    {
        write_pos = 0;
        string[0] = '\000';
    } else
    {
        p = item;

        while(*p) string[write_pos++] = *p++;

        if(separator)
        {
            string[write_pos++] = separator;
        } else
        {
            string[write_pos] = '\000';
        }
    }

    return string;
}



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

/*  check if config line is empty or has a comment */
static
int skip_line(char *input_line)
{
    char *c = input_line;
    
    while(isspace(*c)) c++;

    if(*c == '#' || *c == '\000') return 1;

    return 0;
}

/* parse config line. 
 * Line can have heading and trailing whitespaces
 * Comments start with #
 * Format for config is:
 * NAME VALUE
 *
 * Function checks if NAME is present and returns pointer to start of the VALUE
 * returns NULL if NAME not found
 */
static 
char *parse_config_line(char *input_line,char *name)
{
    char *c = input_line,*e;
    int quoted = 0;

    while(isspace(*c)) c++;

    if(strncasecmp(c,name,strlen(name)) == 0)
    {
        c += strlen(name);
        if(isspace(*c))
        {
            while(isspace(*c)) c++;
            if(*c == '\"')
            {
                c++;
                quoted = 1;
            }

            e = c;

            while(*e != '\000')
            {
                if(quoted && *e == '\"')
                    *e = '\000';
                else
                {
                    e++;
                    if(quoted && *e == '\000') quoted = 0;
                }
            }

            e--;

            if(!quoted) while(e > c && isspace(*e)) e--;   // remove trailing spaces
            e++;
            *e = '\000';
            if(*c != '\000') return c;
        }
    }
    return NULL;
}


/* read config file ~/.ceifrc for global default values 
 * for certain variables
 */
#define CONFIG_LEN_MAX 100
void read_config_file(char *config_file)
{
    FILE *f;
    char *value;
    char *home;
    char input_line[CONFIG_LEN_MAX];
    char input_file[1024];

    input_file[0] = '\000';

    home = getenv("HOME");

    if(config_file[0] == '~' && home != NULL)
    {
        strcat(input_file,home);
        strcat(input_file,&config_file[1]);
    } else
    {
        strcat(input_file,config_file);
    }

    f = xfopen_test(input_file,"r",'a');

    if(f == NULL)
    {
        if(strcmp(config_file,CEIF_CONFIG) == 0) return;
        panic("Cannot read rc-file",config_file,NULL);
    }


    while(fgets(input_line,CONFIG_LEN_MAX,f) != NULL)
    {
        if(skip_line(input_line)) continue;

        if((value = parse_config_line(input_line,"SAMPLES")) != NULL)
        {
            samples_max = atoi(value);
        } else if((value = parse_config_line(input_line,"TREES")) != NULL)
        {
            tree_count = atoi(value);
        } else if((value = parse_config_line(input_line,"DECIMALS")) != NULL)
        {
            decimals = atoi(value);
        } else if((value = parse_config_line(input_line,"AUTO_WEIGTH")) != NULL)
        {
            auto_weigth = atoi(value);
        } else if((value = parse_config_line(input_line,"AUTO_SCALE")) != NULL)
        {
            auto_weigth = atoi(value);
        } else if((value = parse_config_line(input_line,"CATEGORY_SEPARATOR")) != NULL)
        {
            category_separator = value[0];
            if(category_separator == '"' && value[1]) category_separator = value[1];
        } else if((value = parse_config_line(input_line,"LABEL_SEPARATOR")) != NULL)
        {
            label_separator = value[0];
            if(label_separator == '"' && value[1]) label_separator = value[1];
        } else  if((value = parse_config_line(input_line,"MAX_SAMPLES")) != NULL)
        {
             max_total_samples = atoi(value);
        } else if((value = parse_config_line(input_line,"OUTLIER_SCORE")) != NULL)
        {
             parse_user_score(value);
        } else if((value = parse_config_line(input_line,"NEAREST")) != NULL)
        {
             if(atoi(value)) nearest = 1;
             else nearest = 0;
        } else if((value = parse_config_line(input_line,"ANALYZE_SAMPLING")) != NULL)
        {
             analyze_sampling_count = atoi(value);
             if(analyze_sampling_count < 0) analyze_sampling_count = 0;
        } else if((value = parse_config_line(input_line,"DEBUG")) != NULL)
        {
             debug = atoi(value);
             if(debug < 0) debug = 0;
        } else if((value = parse_config_line(input_line,"CLUSTER_SIZE")) != NULL)
        {
            cluster_relative_size = atof(value);
            if(isnan(cluster_relative_size) || cluster_relative_size < 0.0 || cluster_relative_size > 1.0) cluster_relative_size = 0.25;
        } else if((value = parse_config_line(input_line,"PRINT_DIMENSION")) != NULL)
        {
            if(print_dimension != NULL) free(print_dimension);
            print_dimension = xstrdup(value);
        } else if((value = parse_config_line(input_line,"DIM_PRINT_WIDTH")) != NULL)
        {
            dimension_print_width = atoi(value);
            if(dimension_print_width <= 0) dimension_print_width = 25;
        } else
        {
             panic("Unknown option in config file",input_line,NULL);
        }
    }
}


#define _I fprintf(outs,"  ")
#define _P(...) fprintf(outs, __VA_ARGS__)
#define _2P(...) _I; _P( __VA_ARGS__)
#define _3P(...) _I; _2P( __VA_ARGS__)
#define _O(b) (b ? "On" : "Off")
#define _S(i,max) if(i < max - 1) _P("%c",',')
#define _D(name,count,array) _2P(name); for(i = 0;i < count;i++) {_P("%d",array[i]+1);_S(i,count);} _P("\n")

/* Print forest contents in user readable form
 */
void
print_forest_info(FILE *outs)
{
    int forest_idx,i,j;
    struct forest *f;
    char outstr[100];
    struct tm *tmp;
        
    _P("Global setting:\n");
    _2P("Number of forests: %d\n",forest_count);
    _2P("Number of analyzed dimensions: %d\n",dimensions);
    _2P("Number of samples/tree: %d\n",samples_max);
    _2P("Number of trees: %d\n",tree_count);
    _2P("Number of decimals: %d\n",decimals);

    if(percentage_score)
    {
        _2P("Outlier score is the score under which there are %.2f percent of sample scores\n",outlier_score); 
    }
    else
    {
        _2P("Outlier score: %f",outlier_score);

        if(scale_score)
        {
            _P(", scores are scaled to 0..1 using forest sample minimum and maximum score\n");
        } else
        {
            _2P("\n");
        }
    }

    _2P("Relative cluster size: %f\n",cluster_relative_size);
    _2P("Input separator: %c\n",input_separator);
    _2P("Output separator: %c\n",list_separator);
    _2P("Header is %s\n",_O(header));
    _2P("Automatic data value scaling is %s\n",_O(auto_weigth));
    _2P("Aggregate is %s\n",_O(aggregate));
    _2P("Unigue samples is %s\n",_O(unique_samples));
    _2P("Nearest distance analysis is %s\n",_O(nearest));
    _2P("Print string: \"%s\"\n",print_string);

    _2P("\n");

    _D("Dimensions used in analysis: ",dimensions,dim_idx);
    _D("User ignored dimensions: ",ignore_idx_count,ignore_idx);
    _D("User included dimensions: ",include_idx_count,include_idx);
    _D("Category dimensions: ",category_idx_count,category_idx);
    _D("Label dimensions: ",label_idx_count,label_idx);
    _D("Dimensions treated as text: ",text_idx_count,text_idx);
    _P("\n");
    _2P("Density is sample max - min range divided by sample count");
    _P("\n");
    _2P("Cluster coverage is ratio between 0 - 1, where 1 = clusters cover all samples");
    _P("\n");


    if(forest_count) _P("\nForest data:\n");
    for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
    {
        f = &forest[forest_idx];

        calculate_forest_score(forest_idx);

        _2P("\nForest category string: \'%s\'\n",f->category);
        _3P("Filter is %s\n",_O(f->filter));
        _3P("Number of samples: %d\n",f->X_count);

        if(!f->filter)
        {
            _3P("Average path length (c): %f\n",f->c);
            _3P("Max. tree heigth: %d\n",f->heigth_limit);

            if(scale_score)
            { 
                _3P("Forest score range is between %f and %f, this is used to scale data scores to 0..1 range\n",f->min_score,f->max_score);
            } else if(percentage_score)
            {
                _3P("Percentage based score: %f, %.2f%% of samples have lower score\n",f->percentage_score,outlier_score);
            }
        }

        if(nearest)
        {
            _3P("Average%ssample point distance for a single tree: %f\n",auto_weigth ? " scaled " : " ",f->avg_sample_dist);
        }

        tmp = localtime(&f->last_updated);
        if(tmp != NULL)
        {
            outstr[0] = '\000';
            strftime(outstr, sizeof(outstr), "%c", tmp);
            _3P("Last updated: %s\n",outstr);
        }

        if(!f->X_count) continue;

        _P("\n");
        _3P("%15s%*s\n","Dimension sample value summary:",dimensions * dimension_print_width / 2,"Dimension");

        _3P("%15s","");
        for(i = 0;i < dimensions;i++) _P("%*d",dimension_print_width,i + 1);
        _P("\n");

        _3P("%15s","Maximum value");
        for(i = 0;i < dimensions;i++) _P("%*.*f",dimension_print_width,decimals,f->max[i]);
        _P("\n");
        
        _3P("%15s","Minimum value");
        for(i = 0;i < dimensions;i++) _P("%*.*f",dimension_print_width,decimals,f->min[i]);
        _P("\n");
        
        _3P("%15s","Average value");
        for(i = 0;i < dimensions;i++) _P("%*f",dimension_print_width,f->avg[i]);
        _P("\n");
        
        _3P("%15s","Density");
        for(i = 0;i < dimensions;i++) _P("%*f",dimension_print_width,f->dim_density[i]);

        _P("\n");

        if(!f->filter && cluster_relative_size > 0.0)
        {
            _P("\n");
            _3P("Forest cluster centers, cluster radius: %f, cluster coverage: %f\n",f->cluster_radius,f->cluster_coverage);
            _3P("%15s\n","Cluster number");
            for(i = 0;i < f->cluster_count;i++)
            {
                _3P("%15d",i + 1);
                for(j = 0;j < dimensions;j++) _P("%*.*f",dimension_print_width,decimals,f->X[f->cluster_center[i]].dimension[j]);
                _P("\n");
            }
        }

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
            _P("%15d ",j + 1);

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

/* print samples scores. 
 * All samples of all non filterd forest are printed with sample score
 */

void
print_sample_scores(FILE *outs)
{
    int forest_idx,i,j;
    double score;
    struct forest *f;

    _P("Sample score list\n");

    for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
    {
       f = &forest[forest_idx];

       if(f->filter) continue;

       calculate_forest_score(forest_idx);
        
       _P("\nForest category string: %s\n",f->category);
       _2P("%10s%*s\n","Score",dimensions * dimension_print_width / 2,"Dimension values");
       _2P("%10s","");
       for(j = 0;j < dimensions;j++) _P("%*d",dimension_print_width,j + 1);
       _P("\n");

       for(i = 0;i < f->X_count;i++)
       {
           score = sample_score_scale(forest_idx,&f->X[i]);
           _2P("%10f",score);
           for(j = 0;j < dimensions;j++) _P("%*.*f",dimension_print_width,decimals,f->X[i].dimension[j]);
           _P("\n");
       }
    }
}

/* calculate standard deviation for all dimensions for a forest
 */
static
void calc_stddev(struct forest *f, double *stddev)
{
    int i,j;

    for(i = 0;i < dimensions;i++) stddev[i] = 0.0;

    for(i = 0;i < f->X_count;i++)
    {
        for(j = 0;j < dimensions;j++) stddev[j] += POW2(f->X[i].dimension[j] - f->avg[j]); 
    }

    for(j = 0;j < dimensions;j++) stddev[j] = sqrt(stddev[j] / (double) (f->X_count - 1));
}
        

/* print correlation coefficents for all dimension attribute pairs using Pearson method. 
 * Correlation coefficent has a value between +1 and −1. A value of +1 is total positive linear correlation, 0 is no linear correlation, and −1 is total negative linear correlation
 * */
void print_correlation_coefficent(FILE *outs)
{
    int forest_idx,i;
    int a,b;
    double psum,cc,slope;
    double *stddev;
    struct forest *f;

    if(dimensions < 2) return;  // This makes sence only for multi-dim cases

    stddev = xmalloc(dimensions * sizeof(double));

    _P("Correlation coefficent with regression line slope and y-intercept for every dimension attribute pair.\n");
    _P("Correlation coefficent has a value between +1 and −1. A value of +1 is total positive linear correlation, 0 is no linear correlation, and −1 is total negative linear correlation.\n");
    _P("Value 0 is also returned in case the correlation coefficent is undefined.\n");

    for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
    {
        f = &forest[forest_idx];

        if(f->filter) continue;
        if(f->X_count < 2) continue;

        _P("\nForest category string: %s, number of samples: %d\n",f->category,f->X_count);
        _2P("%12s %15s %15s%15s%15s\n","Coefficent","Slope","y-intercept","Dimension x","Dimension y");

        calc_stddev(f,stddev);

        for(a = 0;a < dimensions;a++)
        {
            for(b = a + 1;b < dimensions;b++)
            {
                if(stddev[a] > 0.0 && stddev[b] > 0.0)
                {
                    psum = 0.0;

                    for(i = 0;i < f->X_count;i++)
                    {
                        psum += f->X[i].dimension[a] * f->X[i].dimension[b];
                    }

                    cc = (psum - (double) f->X_count * f->avg[a] * f->avg[b]) / ((double) (f->X_count - 1) * stddev[a] * stddev[b]);
                    slope = cc * (stddev[b] / stddev[a]);
                } else
                {
                    cc = 0.0;
                    slope = 0.0;
                }
                
                _2P("%12f %15f %15f%15d%15d\n",cc,slope,f->avg[b] - slope * f->avg[a],a + 1,b + 1);
            }
        }
    }
    free(stddev);
}

