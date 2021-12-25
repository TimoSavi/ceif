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
#include <locale.h>
#include <errno.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>

/* Global data */
int dim_idx[DIM_MAX];          // final table of dimension indices to be used. Index refers to input line field index
int dimensions = 0;            // dimensions in current setup

int text_idx[DIM_MAX];         // table of dimension indices having text based input values. Texts are mapped to hash values using hash().
int text_idx_count = 0;        // number of text based input values

int ignore_idx[DIM_MAX];
int ignore_idx_count = 0;

int include_idx[DIM_MAX];
int include_idx_count = 0;

int category_idx[DIM_MAX];      // final table of category indices to be used. Index refers to input line field index
int category_idx_count  = 0;             // number of category fields

int label_idx[DIM_MAX];         // final table of label indices to be used. Index refers to input line field index
int label_idx_count = 0;                 // number of labels fields

int score_idx[DIM_MAX];        // table of dimensions indices which must have high score along total score. If none if these dims. dooes not have high score then the data is not outlier
int score_idx_count = 0;

char *cat_filter[FILTER_MAX];  // Category filters
int cat_filter_count = 0;      // Category filter count

int auto_weigth = 1;           // Should weigths be calculated automatically

char *print_string = NULL;     // How to print outlier data
char *print_dimension = NULL;  // How to print directive %m (in print_string). This combines different dimension values together
int tree_count = 100;             // trees / forest
int samples_max = 256;            // max samples / tree
int max_total_samples = 0;        // limit for samples_total, if zero use samples_max * trees
int samples_total;                // max samples / forest
char input_separator = ',';       // input separator for csv data
char category_separator = ';';    // separator for category values
char label_separator = '-';       // separator for category values
int header = 0;                   // input data has a header row to skip
double outlier_score = 0.5;      // outlier score
int decimals = 6;                 // Number of decimals when printing and saving dimension data
int unique_samples = 0;           // accept only unique samples, in some cases this yields better results
char *printf_format = "";       // User given printf format for dimension and average values
char list_separator = ',';         // seprator for dimension and average values in output
int n_vector_adjust = 0;        // should n vector to be adjust among data set
int aggregate = 0;              // should data values to be aggregated when adding new data to forest
int scale_score = 1;               // should outlier scores be scaled between foretsts, scaled score is between 0..1
int percentage_score = 0;          // outlier score is based on training data distribution, score is the largest score of the x% set of samples having the smallest score
int nearest = 1;                // the shortest distance of analyzed point to nearest sample is calulated in leaf nodes. 
int analyze_sampling_count = 0;    // number of lines / forest after sampling of analyzed lines is started, 0 = sampling disabled
int debug = 0;                     // If set print processing related info
double cluster_relative_size = 0.125; // relative distance for samples in the same cluster, must be between 0 and 1
int dimension_print_width = 25;   // dimension value printing width, used when printing forest info (option -q)

/* User given strings for dim ranges */
char *ignore_dims = "";           // which input values are ignored, user given string
char *include_dims = "";           // which input values are included, user given string
char *category_dims = "";           // list of dimensions to be used as category label, user given string
char *label_dims = "";           // list of dimensions to be used as category label, user given string
char *text_dims = "";           // list of dimensions to be used as text based input values, user given string
char *score_dims = "";           // list of dimensions which should together have high outlier score among total_score, user given string

int forest_count = 0;            // total number of forests
int forest_cap = 0;              // forest capasity in terms of items in forest table
struct forest *forest = NULL;    // forest table

struct forest_hash fhash[HASH_MAX];  // hash table for forest data, speeds search when number of forests is high

static char short_opts[] = "o:hVd:I:t:s:f:l:a:p:w:O:r:C:HSL:U:c:F:T::i:u::m:e:M::D:N::AX:qy::Ekg:Pv:R:z:=j:G:";

#ifdef HAVE_GETOPT_LONG
static struct option long_opts[] =
{
  {"output", 1, 0, 'o'},
  {"help", 0, 0, 'h'},
  {"version", 0, 0, 'V'},
  {"decimals", 1, 0, 'd'},
  {"ignore-dims", 1, 0, 'I'},
  {"include-dims", 1, 0, 'U'},
  {"trees", 1, 0, 't'},
  {"samples", 1, 0, 's'},
  {"input-separator", 1, 0, 'f'},
  {"learn", 1, 0, 'l'},
  {"analyze", 1, 0, 'a'},
  {"print", 1, 0, 'p'},
  {"write-forest", 1, 0, 'w'},
  {"outlier-score", 1, 0, 'O'},
  {"category-dim", 1, 0, 'C'},
  {"header", 0, 0, 'H'},
  {"set-locale", 0, 0, 'S'},
  {"label-dim", 1, 0, 'L'},
  {"read-forest", 1, 0, 'r'},
  {"categorize", 1, 0, 'c'},
  {"category-filter", 1, 0, 'F'},
  {"test", 2, 0, 'T'},
  {"test-interval", 1, 0, 'i'},
  {"unique-samples", 2, 0, 'u'},
  {"printf-format", 1, 0, 'm'},
  {"list-separator", 1, 0, 'e'},
  {"missing", 2, 0, 'M'},
  {"delete", 1, 0, 'D'},
  {"new", 2, 0, 'N'},
  {"aggregate", 0, 0, 'A'},
  {"text-dims", 1, 0, 'X'},
  {"query", 0, 0, 'q'},
  {"sample-density", 2, 0, 'y'},
  {"sample-scores", 0, 0, 'E'},
  {"remove-outlier", 0, 0, 'k'},
  {"rc-file", 1, 0, 'g'},
  {"correlation-coe", 0, 0, 'P'},
  {"average", 1, 0, 'v'},
  {"reset-forest", 1, 0, 'R'},
  {"inplace-forest", 1, 0, 'z'},
  {"sizeof", 0, 0, '='},
  {"print-dimension", 1, 0, 'j'},
  {"score-dims", 1, 0, 'G'},
  {NULL, 0, NULL, 0}
};
#endif

static void
help (int status)
{
  printf ("%s - Categorized extended isolation forest tool\n", PACKAGE_NAME);
  printf ("Usage: %s [OPTION]... \n", PACKAGE_NAME);
  printf ("\
Options:\n\
  -d, --decimals INTEGER      number of decimals when printing and saving dimension values\n\
  -h, --help                  display this help and exit\n\
  -V, --version               output version information and exit\n\
  -I, --ignore-dims LIST      comma separated list of dimensions not to be used, first is number 1. Ranges can be given using dash\n\
  -U, --use-dims LIST         comma separated list of dimensions to be used, first is number 1. Ranges can be given using dash. Overwrites entries from -I option\n\
  -t, --trees INTEGER         number of trees. default is 100\n\
  -s, --samples INTEGER       number of samples/tree. Default is 256\n\
  -f, --input-separator CHAR  input file field separator. Default is comma\n\
  -l, --learn FILE            file to used for training \n\
  -a, --analyze FILE          file to analyze\n\
  -c, --categorize FILE       file to categorize\n\
  -p, --print STRING          outlier printing format\n\
  -j, --print-dimension STRING print format format for printing directive %m, prints joined dimension values\n\
  -o, --output FILE           outlier data is printed to FILE. Default is stdout\n\
  -w, --write-forest FILE     write forest data to FILE\n\
  -O, --outlier-score FLOAT   outlier data is printed if score is bigger that FLOAT (0.0 - 1.0)\n\
  -O, --outlier-score FLOATs  outlier data is printed if score is bigger that FLOAT (0.0 - 1.0), actual scores are scaled to range 0..1\n\
  -O, --outlier-score FLOAT%%  outlier score is the score which covers FLOAT percent of samples, give value between 0 - 100\n\
  -r, --read-forest FILE      read forest data from FILE\n\
  -z, --inplace-forest FILE   read forest data from FILE and after any processing write forest back to FILE\n\
  -C, --category-dim LIST     comma separated list of dimensions to form a category string\n\
  -L, --label-dim LIST        comma separated list of dimensions to form a label string\n\
  -H, --header                input data file has a header\n\
  -S, --set-locale            locale information is read from environment\n\
  -T, --test FLOAT            generate test data with adjustment factor FLOAT\n\
  -i, --test-interval INTEGER number of test points for each dimension, default is 256. Used with option -T\n\
  -F, --category-filter REGEXP regular expression to filter categories. Several option can be given. If REGEXP starts with \"-v \" the matching  is inverted\n\
  -u, --unique-samples INTEGER accept INTEGER percent of samples as duplicates, default is take all samples.\n\
  -m, --printf-format STRING  printf format string for dimension and average value printing\n\
  -e, --list-separator CHAR   value separator for dimension and average value printing\n\
  -M, --missing STRING        print category value of forests which have not used in analysis. Optional printf format STRING is used for printing\n\
  -D, --delete INTEGER        before saving the forest data to file delete those forests which have not been updated INTEGER (seconds) ago\n\
  -N, --new STRING            print values which do not match any known category. Optional printf format STRING is used for printing\n\
  -A, --aggregate             instead taking samples as they are, aggregate new samples by adding values for each forest. Only one new aggregated sample for each forest is added for each usage of -l option\n\
  -X, --text-dims STRING      comma separated list of dimensions to be used as text based input values, first is number 1. Ranges can be given using dash\n\
  -G, --score-dims STRING     comma separated list of dimensions attribute indices. Combination of these dimension attributes must have outlier score along total score. Ranges can be given using dash. These are dimensions attribute indices, not input line indices (first is always number 1)\n\
  -q, --query                 print forest info and exit\n\
  -y, --sample-density        print ascii map of all forest sample value densities and exit\n\
  -yy, --sample-densityy      print ascii map of all forest sample value densities using common scale for all forests and exit\n\
  -E, --sample-scores         print samples values with sample score and exit\n\
  -k, --remove-outlier        remove the sample having largest outlier score. For each invocation of this option one sample is removed\n\
  -g, --rc-file FILE          read global settings from FILE (default is ~/.ceifrc)\n\
  -P, --correlation_coe       print list of correlation coefficents with regression line slopes and y-intercepts for every dimension attribute pair and exit. Correlation coefficent is a value between -1.0 - 1.0\n\
  -v, --average STRING        print average info for each forest after analysis using STRING as print format\n\
  -R, --reset-forest STRING   remove all samples for a forest read using option -r and having forest string STRING\n\
");
  printf ("\nSend bug reports to %s\n", PACKAGE_BUGREPORT);
  exit (status);
}

static
time_t parse_delete_interval(char *s)
{
    int len = strlen(s);
    time_t value = (time_t) atol(s);

    if(value == (time_t) 0) panic("Invalid time format for old forest data deletion",s,NULL);

    switch(s[len - 1])
    {
        case 'Y':
        case 'y':
            value *= 31556926;
            break;
        case 'M':
            value *= 2629743;
            break;
        case 'D':
        case 'd':
            value *= 86400;
            break;
        case 'm':
            value *= 60;
            break;
        case 's':
            break;
        default:
            if(s[len - 1] < '0' || s[len - 1] > '9') panic("Invalid time format for old forest data deletion",s,NULL);
            break;
    }
    return value;
}


static
void usage(int opt)
{
        fprintf(stderr,"Unknown option '-%c'\n",(char) opt);
        help(1);
}

static 
void message(char *msg,char *info,char *syserror)
{
    if(msg != NULL)
    {
        if (info == NULL && syserror == NULL)
        {
            fprintf(stderr,"%s: %s\n",PACKAGE_NAME,msg);
        } else if(info != NULL && syserror == NULL)
        {
            fprintf(stderr,"%s: %s: %s\n",PACKAGE_NAME,msg,info);
        } else if(info != NULL && syserror != NULL)
        {
            fprintf(stderr,"%s: %s: %s; %s\n",PACKAGE_NAME,msg,info,syserror);
        } else if(info == NULL && syserror != NULL)
        {
            fprintf(stderr,"%s: %s; %s\n",PACKAGE_NAME,msg,syserror);
        }
    }
}

void panic(char *msg,char *info,char *syserror)
{
    message(msg,info,syserror);
    exit(1);
}

void info(char *msg,char *info,char *syserror)
{
    message(msg,info,syserror);
}


void
print_version()
{
    printf("%s version %s\n",PACKAGE_NAME,PACKAGE_VERSION);
#if defined HAVE_JSON_JSON_H || defined HAVE_LIBFASTJSON_JSON_H || defined HAVE_JSON_C_JSON_H
    printf("Forest data is stored in JSON format\n");
#endif
    printf("Copyright (c) 2020 Timo Savinen\n\n");
    printf("This is free software; see the source for copying conditions.\n");
    printf("There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
}


/* Init forest hash table
 */
static
void init_forest_hash()
{
    int i;

    for(i = 0;i < HASH_MAX;i++)
    {
        fhash[i].idx_count = 0;
        fhash[i].idx_cap = 0;
        fhash[i].idx = NULL;
    }
}

/* parse outlier score
 */
void parse_user_score(char *score_str)
{
    char *endp;

    scale_score = 0;
    percentage_score = 0;

    outlier_score = strtod(score_str,&endp);

    if(*endp == 's' && endp[1] == '\000') 
    {
        scale_score = 1;
        *endp = '\000';
    } else if(*endp == '%' && endp[1] == '\000')
    {
        percentage_score = 1;
        *endp = '\000';
    }

    if(percentage_score)
    {
        if(outlier_score < 0.0 || outlier_score > 100.0) 
        {
            panic("Give percentage based score between 0 and 100",NULL,NULL);
        }
    } else if(outlier_score < 0.0 || outlier_score > 1.0 || *endp != '\000') 
        panic("Give outlier score between 0 and 1 (with suffix \'s\' if scaling is required) or between 0 and 100 with suffix \'%\' for percentage score",NULL,NULL);
}

int
main (int argc, char **argv)
{
    int opt;
    int set_locale = 0;
    int run_test = 0;
    int make_tree = 0;
    int make_query = 0;
    int print_density = 0;
    int print_sample_s = 0;
    int print_correlation = 0;
    int print_average = 0;
    int kill_outlier = 0;
    int common_scale = 0;
    int test_range_interval = 256;
    int print_missing = 0;
    time_t delete_interval = (time_t) 0;
    char *missing_format = "%C";
    char *average_format = NULL;
    char *not_found_format = NULL;
    double test_extension_factor = 0.0;    // extents the area from where test sample points are selected
    int score_option_given = 0;
    char *learn_file = NULL;
    char *analyze_file = NULL;
    char *categorize_file = NULL;
    char *save_file = NULL;
    char *load_file = NULL;
    char *output_file = NULL;
    FILE *learns = NULL;            // file to read learn data
    FILE *analyzes = NULL;          // file to analyze
    FILE *categorizes = NULL;          // file to categorize
    FILE *loads = NULL;           // file to read saved forest data
    FILE *outs = NULL;           // file to print results

    atexit(print_alloc_debug);

    setlocale(LC_ALL,"C");

    init_forest_hash();
    read_config_file(CEIF_CONFIG);          // config file parameters are read before options

#ifdef HAVE_GETOPT_LONG
    while ((opt = getopt_long(argc,argv,short_opts,long_opts,NULL)) != -1)
#else
        while ((opt = getopt(argc,argv,short_opts)) != -1)
#endif
        {
            switch(opt)
            {
                case 'd':
                    decimals = atoi(optarg);
                    break;
                case 'I':
                    ignore_dims = xstrdup(optarg);
                    ignore_idx_count = parse_dims(optarg,ignore_idx);
                    break;
                case 'U':
                    include_dims = xstrdup(optarg);
                    include_idx_count = parse_dims(optarg,include_idx);
                    break;
                case 't':
                    tree_count = atoi(optarg);
                    if(tree_count < 2) panic("Tree count less than two makes no sense",NULL,NULL);
                    break;
                case 's':
                    samples_max = atoi(optarg);
                    if(samples_max < SAMPLES_MIN) panic("Low sample count makes no sense",NULL,NULL);
                    break;
                case 'f':
                    input_separator = optarg[0];
                    break;
                case 'l':
                    learn_file = xstrdup(optarg);
                    break;
                case 'a':
                    analyze_file = xstrdup(optarg);
                    break;
                case 'c':
                    categorize_file = xstrdup(optarg);
                    break;
                case 'p':
                    if(print_string != NULL) free(print_string);
                    print_string = xstrdup(optarg);
                    break;
                case 'j':
                    if(print_dimension != NULL) free(print_dimension);
                    print_dimension = xstrdup(optarg);
                    break;
                case 'O':
                    parse_user_score(optarg);
                    score_option_given = 1;
                    break;
                case 'w':
                    save_file = xstrdup(optarg);
                    break;
                case 'z':   /* no break here */
                    if(save_file == NULL) save_file = xstrdup(optarg);
                case 'r':
                    load_file = xstrdup(optarg);
                    DEBUG("*** Loading forest data from %s\n",load_file);

                    /* load now, parameters after this take higher presence */
                    if(forest_count == 0)
                    {
                        if (opt == 'z')    // test file readbility
                        {
                            loads = xfopen_test(load_file,"r",'a');
                        } else
                        {
                            loads = xfopen(load_file,"r",'a');
                        }

                        if(loads != NULL)
                        {
                            fclose(loads);
                            if(!read_forest_file(load_file)) panic("Cannot load forest data from file",load_file,NULL);
                        }
                    } 
                    break;
                case 'C':
                    category_dims = xstrdup(optarg);
                    category_idx_count = parse_dims(optarg,category_idx);
                    break;
                case 'L':
                    label_dims = xstrdup(optarg);
                    label_idx_count = parse_dims(optarg,label_idx);
                    break;
                case 'H':
                    header = 1;
                    break;
                case 'S':
                    set_locale = 1;
                    break;
                case 'o':
                    output_file = xstrdup(optarg);
                    break;
                case '?':
                    help(0);
                    break;
                case 'h':
                    help(0);
                    break;
                case 'F':
                    add_category_filter(optarg);
                    break;
                case 'T':
                    if(optarg != NULL) test_extension_factor = atof(optarg);
                    run_test = 1;
                    break;
                case 'i':
                    test_range_interval = atoi(optarg);
                    break;
                case 'u':
                    unique_samples = 10;
                    if(optarg != NULL) unique_samples = atol(optarg);
                    if(unique_samples < 0 || unique_samples > 100) panic("Give unique sample percent bweteen 0 and 100",NULL,NULL);
                    break;
                case 'm':
                    printf_format = xstrdup(optarg);
                    break;
                case 'e':
                    list_separator = optarg[0];
                    break;
                case 'V':
                    print_version();
                    exit(0);
                    break;
                case 'M':
                    print_missing = 1;
                    if(optarg != NULL) missing_format = xstrdup(optarg);
                    break;
                case 'D':
                    delete_interval = parse_delete_interval(optarg);
                    break;
                case 'N':
                    if(optarg != NULL) 
                    {
                        not_found_format = xstrdup(optarg);
                    } else
                    {
                        not_found_format = "%v";
                    }
                    break;
                case 'A':
                    aggregate = 1;
                    break;
                case 'X':
                    text_dims = xstrdup(optarg);
                    text_idx_count = parse_dims(optarg,text_idx);
                    break;
                case 'G':
                    score_dims = xstrdup(optarg);
                    score_idx_count = parse_dims(optarg,score_idx);
                    break;
                case 'q':
                    make_query = 1;
                    break;
                case 'y':
                    print_density = 1;
                    if(optarg != NULL && optarg[0] == 'y') common_scale = 1;
                    break;
                case 'E':
                    print_sample_s = 1;
                    break;
                case 'k':
                    kill_outlier++;
                    break;
                case 'g':
                    read_config_file(optarg);
                    break;
                case 'P':
                    print_correlation = 1;
                    break;
                case 'v':
                    print_average = 1;
                    if(optarg != NULL)
                    {
                        average_format = xstrdup(optarg);
                    } else
                    {
                        average_format = "%C %r %h";
                    }
                    break;
                case 'R':
                    remove_samples(optarg);
                    break;
                case '=':
                    printf("sizeof double: %u\n",(unsigned int) sizeof(double));
                    printf("sizeof int: %u\n",(unsigned int) sizeof(int));
                    exit(0);
                    break;
                default:
                    usage(opt);
                    break;
            }
        }

    srand(time(NULL) + getpid());

    init_fast_n_cache();
    init_fast_c_cache();

    if(set_locale) setlocale(LC_ALL,"");

    samples_total = max_total_samples ?  max_total_samples : tree_count * samples_max;   // total samples count is trees * samples/tree, this can be limited using config MAX_SAMPLES

    if(print_string == NULL) print_string = "%s %v";
        
    if(analyze_file !=  NULL || categorize_file !=  NULL || run_test || make_query || print_sample_s || 
       kill_outlier || print_correlation || print_average) make_tree = 1;  // we need tree info

    if(output_file != NULL)
    {
        outs = xfopen(output_file,"w",'a');

    } else
    {
        outs = stdout;
    }

    if(forest_count) 
    {
         train_forest(NULL,1,make_tree); // samples read allready from saved file, run training based on that
    } else
    {
        if(learn_file != NULL) 
        {
            DEBUG("\n***read training data from file %s\n",learn_file);
            learns = xfopen(learn_file,"r",'a');
            train_forest(learns,1,make_tree); 
            fclose(learns);
            free(learn_file);
            learn_file = NULL;
        } 
    }

    while(kill_outlier--) remove_outlier();

    if(print_density)
    {
        print_sample_density(outs,common_scale);
        fclose(outs);
        exit(0);
    }

    if(analyze_file !=  NULL)
    {
        analyzes = xfopen(analyze_file,"r",'a');
        analyze(analyzes,outs,not_found_format,average_format);
        fclose(analyzes);
        if(print_missing) print_missing_categories(outs,missing_format);
    }

    if(categorize_file !=  NULL)
    {
        categorizes = xfopen(categorize_file,"r",'a');
        categorize(categorizes,score_option_given,outs);
        fclose(categorizes);
    }

    if(learn_file != NULL) 
    {
        learns = xfopen(learn_file,"r",'a');
        train_forest(learns,forest_count ? 0 : 1,0); 
        fclose(learns);
    } 

    if(make_query)
    {
        print_forest_info(outs);
        fclose(outs);
        exit(0);
    }
    
    if(print_sample_s)
    {
        print_sample_scores(outs);
        fclose(outs);
        exit(0);
    }

    if(print_correlation)
    {
        print_correlation_coefficent(outs);
        fclose(outs);
        exit(0);
    }

    if(save_file != NULL)
    {
        if(set_locale) setlocale(LC_ALL,"C");
        write_forest_file(save_file,delete_interval);
        if(set_locale) setlocale(LC_ALL,"");
    }

    if(run_test)
    {
        test2(outs,test_extension_factor,test_range_interval);
    }

    fclose(outs);

    exit(0) ;
}



