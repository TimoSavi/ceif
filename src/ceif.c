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
int category_idx[DIM_MAX];           // final table of category indices to be used. Index refers to input line field index
int categories = 0;             // number of category fields

int dims_ignore[DIM_MAX];       // Which dimensions are not processed, bit map type of table
int dims_category[DIM_MAX];     // Which dimensions are used as category label, bit map type of tabe
int label_dim = -1;             // which dimension is the label dimension

char *cat_filter[FILTER_MAX];  // Category filters
int cat_filter_count = 0;      // Category filter count

char *print_string="%s %v";     // How to print outlier data
int tree_count = 100;             // trees / forest
int samples_max = 256;            // max samples / tree
int samples_total;                // max samples / forest
char input_separator = ',';       // input separator for csv data
int header = 0;                         // input data has a header row to skip
double outlier_score = 0.75;            // outlier score
double prange_extension_factor = 1.0;    // extents the area from where p is selected
char *ignore_dims = NULL;           // which input values are ignored, user given string
char *include_dims = NULL;           // which input values are included, user given string
char *category_dims = NULL;           // list of dimensions to be used as category label, user given string

int forest_count = 0;            // total number of forests
int forest_cap = 0;              // forest capasity in terms of items in forest table
struct forest *forest = NULL;    // forest table




static char short_opts[] = "o:hVDI:t:s:f:l:a:p:w:O:r:C:HSL:R:U:c:F:";

#ifdef HAVE_GETOPT_LONG
static struct option long_opts[] =
{
  {"output", 1, 0, 'o'},
  {"help", 0, 0, 'h'},
  {"version", 0, 0, 'V'},
  {"debug", 0, 0, 'D'},
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
  {"load-forest", 1, 0, 'L'},
  {"category-dim", 1, 0, 'C'},
  {"header", 0, 0, 'H'},
  {"set-locale", 0, 0, 'S'},
  {"label-dim", 1, 0, 'L'},
  {"save-forest", 1, 0, 'd'},
  {"read-forest", 1, 0, 'r'},
  {"prange-factor", 1, 0, 'R'},
  {"categorize", 1, 0, 'c'},
  {"category-filter", 1, 0, 'F'},
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
  -d, --debug                 print debug info\n\
  -h, --help                  display this help and exit\n\
  -V, --version               output version information and exit\n\
  -I, --ignore-dims LIST      comma separated list of dimensions not to be used, first is number 1. Ranges can be given using dash\n\
  -U, --use-dims LIST         comma separated list of dimensions to be used, first is number 1. Ranges can be given using dash. Overwrites entries from -I option\n\
  -t, --trees INTEGER         number of trees. default is 100\n\
  -s, --samples INTEGER       number of samples. default is 256\n\
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
  -F, --category-filter REGEXP Regular expression to filter categories\n\
");
  printf ("\nSend bug reports to %s\n", PACKAGE_BUGREPORT);
  exit (status);
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
    printf("Copyright (c) 2019 Timo Savinen\n\n");
    printf("This is free software; see the source for copying conditions.\n");
    printf("There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
}


int
main (int argc, char **argv)
{
    int opt,i;
    int set_locale = 0;
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

    for(i = 0;i < DIM_MAX;i++)   // reset dim "bit map" tables 
    {
        dims_ignore[i] = 0;
        dims_category[i] = 0;
    }

    #ifdef HAVE_GETOPT_LONG
    while ((opt = getopt_long(argc,argv,short_opts,long_opts,NULL)) != -1)
    #else
    while ((opt = getopt(argc,argv,short_opts)) != -1)
    #endif
    {
        switch(opt)
        {
            case 'I':
                ignore_dims = xstrdup(optarg);
                parse_dims(optarg,1,dims_ignore);
                break;
            case 'U':
                include_dims = xstrdup(optarg);
                parse_dims(optarg,0,dims_ignore);
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
                    loads = xfopen_test(load_file,"r",'a');
                    if(loads) 
                    {
                        if(!read_forest_file(loads)) panic("Cannot load forest data",load_file,NULL);
                        fclose(loads);
                    }
                }
                break;
            case 'C':
                category_dims = xstrdup(optarg);
                parse_dims(optarg,1,dims_category);
                break;
            case 'L':
                label_dim = atoi(optarg) - 1;
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
            case 'V':
                print_version();
                exit(0);
                break;
            default:
                usage(opt);
                break;
        }
    }

    srand(time(NULL));

    set_locale ? setlocale(LC_ALL,"") : setlocale(LC_ALL,"C");

    samples_total = tree_count * samples_max;

    if(learn_file != NULL) 
    {
        learns = xfopen(learn_file,"r",'a');
        train_forest(learns,forest_count ? 0 : 1); 
        fclose(learns);
    } else
    {
        if(forest_count) train_forest(NULL,1); // samples read allready from saved file, run training based on that
    }

    if(save_file != NULL)
    {
        if(set_locale) setlocale(LC_ALL,"C");
        saves = xfopen(save_file,"w",'a');
        write_forest_file(saves);
        fclose(saves);
        if(set_locale) setlocale(LC_ALL,"");
    }

    if(output_file != NULL)
    {
        outs = xfopen(output_file,"w",'a');

    } else
    {
        outs = stdout;
    }

    if(analyze_file !=  NULL)
    {
        analyzes = xfopen(analyze_file,"r",'a');
        analyze(analyzes,outs);
        fclose(analyzes);
    }

    if(categorize_file !=  NULL)
    {
        categorizes = xfopen(categorize_file,"r",'a');
        categorize(categorizes,outs);
        fclose(categorizes);
    }

    fclose(outs);

    exit(0) ;
}



