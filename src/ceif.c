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
int dim_idx[DIM_MAX];           // final table of dimension indices to be used. Index refers to input line field index
int dimensions = 0;             // dimensions in current setup

int ignore_idx[DIM_MAX];
int ignore_idx_count = 0;

int include_idx[DIM_MAX];
int include_idx_count = 0;

int category_idx[DIM_MAX];      // final table of category indices to be used. Index refers to input line field index
int category_idx_count  = 0;             // number of category fields

int label_idx[DIM_MAX];         // final table of label indices to be used. Index refers to input line field index
int label_idx_count = 0;                 // number of labels fields

char *cat_filter[FILTER_MAX];  // Category filters
int cat_filter_count = 0;      // Category filter count

char *print_string = NULL;     // How to print outlier data
int tree_count = 100;             // trees / forest
int samples_max = 256;            // max samples / tree
int samples_total;                // max samples / forest
char input_separator = ',';       // input separator for csv data
int header = 0;                         // input data has a header row to skip
double outlier_score = 0.75;            // outlier score
double prange_extension_factor = 1.0;    // extents the area from where p is selected
int decimals = 6;                 // Number of decimals when printing saving dimension data
int unique_samples = 0;           // accpet only unique samples, in some cases this yields better results
char *printf_format = "";       // User given printf format for dimension and average values
char list_separator = ',';         // seprator for dimension and average values in output
int n_vector_adjust = 0;        // should n vector to be adjust among data set
int aggregate = 0;              // should data values to be aggregated when adding new data to forest

/* User given strings for dim ranges */
char *ignore_dims = NULL;           // which input values are ignored, user given string
char *include_dims = NULL;           // which input values are included, user given string
char *category_dims = NULL;           // list of dimensions to be used as category label, user given string
char *label_dims = NULL;           // list of dimensions to be used as category label, user given string

int forest_count = 0;            // total number of forests
int forest_cap = 0;              // forest capasity in terms of items in forest table
struct forest *forest = NULL;    // forest table

struct forest_hash fhash[HASH_MAX];  // hash table for forest data, speeds search when number of forests is high

static char short_opts[] = "o:hVd:I:t:s:f:l:a:p:w:O:r:C:HSL:R:U:c:F:T::i:u::m:e:nM::D:N::A";

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
  {"save-forest", 1, 0, 'd'},
  {"read-forest", 1, 0, 'r'},
  {"prange-factor", 1, 0, 'R'},
  {"categorize", 1, 0, 'c'},
  {"category-filter", 1, 0, 'F'},
  {"test", 2, 0, 'T'},
  {"test-interval", 1, 0, 'i'},
  {"unique-samples", 2, 0, 'u'},
  {"printf-format", 1, 0, 'm'},
  {"list-separator", 1, 0, 'e'},
  {"n-adjust", 0, 0, 'n'},
  {"missing", 2, 0, 'M'},
  {"delete", 1, 0, 'D'},
  {"new", 2, 0, 'N'},
  {"aggregate", 0, 0, 'A'},
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
  -s, --samples INTEGER       number of samples/tree. default is 256\n\
  -f, --input-separator CHAR  input file field separator. Default is comma\n\
  -l, --learn FILE            file to used for training \n\
  -a, --analyze FILE          file to analyze\n\
  -c, --categorize FILE       file to categorize\n\
  -p, --print STRING          outlier printing format\n\
  -o, --output FILE           outlier data is printed to FILE. Default is stdout\n\
  -w, --write-forest FILE     write forest data to FILE\n\
  -O, --outlier-score FLOAT   outlier data is printed if score is bigger that FLOAT (0-1.0)\n\
  -r, --read-forest FILE      read forest data from FILE\n\
  -C, --category-dim INTEGER  categore dimensions number\n\
  -L, --label-dim INTEGER     label dimensions number\n\
  -H, --header                input data files have header\n\
  -S, --set-locale            locale information is read from environment\n\
  -R, --prange-factor FLOAT   prange selection adjustment factor\n\
  -T, --test FLOAT            generate test data with adjustment factor FLOAT\n\
  -i, --test-interval INTEGER number of test points for each dimension, default is 256\n\
  -F, --category-filter REGEXP regular expression to filter categories\n\
  -u, --unique-samples INTEGER accept INTEGER percent of samples as duplicates, default is take all samples.\n\
  -m, --printf-format STRING  printf format string for dimension and average value printing\n\
  -e, --list-separator CHAR   value separator for dimension and average value printing\n\
  -n, --n-adjust              adjust n-vector to be perpendicular to dimension attribute having largest value range\n\
  -M, --missing STRING        print category value of forests which have not used in analysis. Optional printf format STRING is used for printing\n\
  -D, --delete INTEGER        before saving the forest data to file delete those forests which have not been updated INTEGER (seconds) ago\n\
  -N, --new STRING            print values which do not match any known category. Optional printf format STRING is used for printing\n\
  -A, --aggregate             instead taking samples as they are, aggregate new samples by adding values for each forest. Only one new aggregated sample for each forest is added for each usage of -l option\n\
");
  printf ("\nSend bug reports to %s\n", PACKAGE_BUGREPORT);
  exit (status);
}

static
time_t parse_delete_interval(char *s)
{
    int len = strlen(s);
    time_t value = (time_t) atol(s);

    if(value == (time_t) 0) panic("Invalid time format for option -M",s,NULL);

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

void panic(char *msg,char *info,char *syserror)
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
    exit(1);
}



void
print_version()
{
    printf("%s version %s\n",PACKAGE_NAME,PACKAGE_VERSION);
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

int
main (int argc, char **argv)
{
    int opt;
    int set_locale = 0;
    int run_test = 0;
    int make_tree = 0;
    int test_range_interval = 256;
    int print_missing = 0;
    time_t delete_interval = (time_t) 0;
    char *missing_format = "%C";
    char *not_found_format = NULL;
    double test_extension_factor = 0.0;    // extents the area from where test sample points are selected
    char *learn_file = NULL;
    char *analyze_file = NULL;
    char *categorize_file = NULL;
    char *save_file = NULL;
    char *load_file = NULL;
    char *output_file = NULL;
    FILE *learns = NULL;            // file to read learn data
    FILE *analyzes = NULL;          // file to analyze
    FILE *categorizes = NULL;          // file to categorize
    FILE *saves = NULL;             // file to save forest for reuse
    FILE *loads = NULL;           // file to read saved forest data
    FILE *outs = NULL;           // file to print results

    init_forest_hash();

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
                if(tree_count < 2) panic("Tree count less than two makes no sence",NULL,NULL);
                break;
            case 's':
                samples_max = atoi(optarg);
                if(samples_max < SAMPLES_MIN) panic("Low sample count makes no sence",NULL,NULL);
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
            case 'w':
                save_file = xstrdup(optarg);
                break;
            case 'O':
                outlier_score = atof(optarg);
                if(outlier_score < 0 || outlier_score > 1) panic("Give outlier score between 0 and 1",NULL,NULL);
                break;
            case 'r':
                load_file = xstrdup(optarg);
                /* load now, parameters after this take higher presence */
                if(load_file !=  NULL && forest_count == 0)
                {
                    setlocale(LC_ALL,"C");
                    loads = xfopen(load_file,"r",'a');
                    if(!read_forest_file(loads)) panic("Cannot load forest data",load_file,NULL);
                    fclose(loads);
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
                setlocale(LC_ALL,"");
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
            case 'R':
                prange_extension_factor = atof(optarg);
                if(prange_extension_factor < 0) panic("Give P-range extension factor equal or larger than zero",NULL,NULL);
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
            case 'n':
                n_vector_adjust = 1;
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
            default:
                usage(opt);
                break;
        }
    }

    srand(time(NULL));

    init_fast_n_cache();
    init_fast_c_cache();

    set_locale ? setlocale(LC_ALL,"") : setlocale(LC_ALL,"C");

    samples_total = tree_count * samples_max;

    if(print_string == NULL) print_string = "%s %v";
        
    if(analyze_file !=  NULL || categorize_file !=  NULL || run_test) make_tree = 1;  // we need tree info

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
            learns = xfopen(learn_file,"r",'a');
            train_forest(learns,1,make_tree); 
            fclose(learns);
            free(learn_file);
            learn_file = NULL;
        } 
    }
    
    if(analyze_file !=  NULL)
    {
        analyzes = xfopen(analyze_file,"r",'a');
        analyze(analyzes,outs,not_found_format);
        fclose(analyzes);
        if(print_missing) print_missing_categories(outs,missing_format);
    }

    if(categorize_file !=  NULL)
    {
        categorizes = xfopen(categorize_file,"r",'a');
        categorize(categorizes,outs);
        fclose(categorizes);
    }

    if(learn_file != NULL) 
    {
        learns = xfopen(learn_file,"r",'a');
        train_forest(learns,forest_count ? 0 : 1,0); 
        fclose(learns);
    } 

    if(save_file != NULL)
    {
        if(set_locale) setlocale(LC_ALL,"C");
        saves = xfopen(save_file,"w",'a');
        write_forest_file(saves,delete_interval);
        fclose(saves);
        if(set_locale) setlocale(LC_ALL,"");
    }

    if(run_test) test2(outs,test_extension_factor,test_range_interval);

    fclose(outs);

    exit(0) ;
}



