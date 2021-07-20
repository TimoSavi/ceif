#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if __STDC__
# define VOID void
#else
# define VOID char
#endif


#if STDC_HEADERS

#include <sys/types.h>
#include <string.h>             /* for strlen etc. */
#include <stdlib.h>

#endif  /* !STDC_HEADERS */

#include <stdio.h>
#include <errno.h>

/* Global constants */

#define DIM_MAX 1024
#define INPUT_LEN_MAX 1048576
#define SAMPLES_MIN 24            // minimun number of samples for a forest
#define FILTER_MAX 100            // maximun number of category filters
#define WEIGTH_MAX 100            // maximun number of weigth options
#define HASH_MAX 32771            // max hash value
#define TEST_SAMPLES 10240        // number of samples when making analysis test
#define N_ADJUST_COUNT 24         // How many N vectors are tested for the best n adjust result
#define NODE_MIN_SAMPLE 3         // Minimum number of samples in node 

/* Special outlier_score values for automatic and average based scores */
#define AVERAGE_SCORE -1
#define AUTO_SCORE -2

/* should normal distributed values be written to cache for faster execution */
#define FAST_N_SAMPLES 32771

/* cache size for c values (the average depth in an unsuccessful search in a Binary Search Tree) */
#define FAST_C_SAMPLES 2048

/* Power of 2 */
#define POW2(a) ((a)*(a))

/* Data structures */
/* All forest related structures will be managed by dynamic tables
 */

struct sample
{
    double *dimension;      // dimension table, dynamically reserved
};

struct node
{
    int sample_count;       // number of samples
    struct sample *samples; // samples array for this node, samples are scaled if auto scale is on
    double *n;               // random normal vector having dimensions count of coordinates
    double pdotn;           // calculate p dot n for performance issues
    int left;               // first left node, -1 if not existing
    int rigth;              // first rigth node, -1 if not existing
};

struct tree
{
    int node_count;         // number of nodes in n
    int node_cap;           // space reserved for node table
    int sample_count;       // number of samples used
    struct node *n;         // table of nodes
    int first;              // index to first node in table.
};

struct forest
{
    char *category;         // Category string, all data having this category string in input data will behandled by this
    int filter;             // true if this forest should not be used in analysis or categorize
    int X_count;	    // Number of samples
    int X_current;          // whis smaple should be taken next to tree 
    int X_cap;              // Memory allocated for X, in terms of units of struct sample
    int X_summary;          // Which sample is used for data summary
    double c;               // Average path length for the forest
    int heigth_limit;       // max tree depth
    struct sample *X;       // samples for this forest
    double *min;            // learn data min values dimension
    double *max;            // learn data max values dimension
    double avg_sample_dist; // Average sample distance in hypercube 
    int scale_range_idx;    // dimension index to range to be used in scaling (-W). Points tomin and max arrays, -1 if no ranges (all attributes have the same value)
    double *avg;            // dimension averages
    double *dim_density;    // learn data average attribute distance dimension
    double *summary;        // aggregated values when analysing or categorizing
    int analyzed;           // Is this forest used in analysis 
    time_t last_updated;    // time when the forest data was last updated in save file. Can be used clean up old forests
    double auto_score;      // automatic socre calculated based on sample value having max score value
    double average_score;   // average score of all saved samples
    double percentage_score;   // percentage based score of saved samples
    double min_score;       // Minimum score of all saved samples
    double max_score;       // Maximum score of all saved samples
    double test_average_score; // average score of tested data, can be compared with average_score
    int total_rows;         // Number of rows read from input
    int analyzed_rows;      // Number of rows used in analysis
    int high_analyzed_rows; // Number of rows having score higher than avaerage score
    int extra_rows;         // Number of rows read after train file after max number of samples reached
    struct tree *t;         // Tree table, NULL if not initialized
};

struct forest_hash
{
    int idx_count;          // number of entries in idx table
    int idx_cap;            // space reserverd for idx
    size_t *idx;               // indices to forest table 
};


/* Global data */
extern int dim_idx[];
extern int ignore_idx[];
extern int include_idx[];
extern int category_idx[];
extern int label_idx[];
extern int text_idx[];

extern int dimensions;            // dimensions in current setup
extern int ignore_idx_count;
extern int include_idx_count;
extern int category_idx_count;
extern int label_idx_count;
extern int text_idx_count;

extern char *cat_filter[];
extern int cat_filter_count;

extern int auto_weigth;      

extern int tree_count;                // trees / forest
extern int samples_max;              // max samples / tree
extern int max_total_samples;
extern int samples_total;              // max samples / forest
extern char *print_string;
extern char input_separator;
extern int header;
extern double outlier_score;
extern double auto_score_factor;
extern double test_extension_factor;
extern int decimals;
extern int unique_samples;
extern char *printf_format;
extern char list_separator;
extern int n_vector_adjust;
extern int aggregate;
extern double average_score_factor;
extern int scale_score;
extern int nearest;
extern int percentage_score;
extern int analyze_sampling_count;

extern char category_separator;       // separator for category values
extern char label_separator;       // separator for category values

extern char *include_dims;
extern char *ignore_dims;
extern char *category_dims;
extern char *label_dims;
extern char *text_dims;

extern int forest_count;
extern int forest_cap;
extern struct forest *forest;           // forest table

extern struct forest_hash fhash[];


/* ceif.c prototypes */
void panic(char *,char *,char *);
void info(char *,char *,char *);
int parse_dims(char *,int *);
void parse_user_score(char *,int);



/* xmalloc.c prototypes */
VOID *xmalloc (size_t);
VOID *xcalloc (size_t, size_t);
VOID *xrealloc (VOID *, size_t);
char *xstrdup (char *);
FILE * xfopen(char *, char *, char);
FILE * xfopen_test(char *, char *, char);

/* file.c prototypes */
char *make_csv_line(char **,int,char);
int parse_csv_line(char **,int,char *,char);
void print_forest_info(FILE *);
void print_sample_density(FILE *,int);
void print_sample_scores(FILE *);
void print_correlation_coefficent(FILE *);
void read_config_file(char *);
char * make_separated_string(char *, char);


/* learn.c prototypes */
void train_forest(FILE *,int,int);
double parse_dim_attribute(char *);
double parse_dim_hash_attribute(char *);
double dot(double *, double *);
double wdot(double *, double *,int, double *,double *);
double c(int);
int dim_ok(int,int);
void add_to_X(struct forest *,char **, int , int);
void add_category_filter(char *);
int search_forest_hash(char *); 
void add_forest_hash(int, char *);
void test2(FILE *,double,int);
void init_fast_n_cache();
void init_fast_c_cache();
double *v_expand(double *,double *,double *,int);
double scale_double(double,double,double,double,double);
void v_copy(double *,double *);
double v_dist_nosqrt(double *,double *);
double v_dist(double *,double *);
double *scale_dimension(double *,struct forest *);
void parse_values(double *,char **, int, int);
int ri(int, int);



/* analyze.c prototypes */
void analyze(FILE *, FILE *,char *,char *);
void categorize(FILE *, int, FILE *);
void init_dims(int);
char *make_category_string(int,char **);
double calculate_score(int ,double *);
void print_missing_categories(FILE *,char *);
void print_(FILE *, double, int,int,int,char **,double *,char *,char *);
int check_idx(int ,int , int *);
void calculate_forest_auto_score(int);
void init_auto_scores();
void remove_outlier();
void calculate_average_sample_score(int);
void calculate_forest_score(int);
double get_forest_score(int);
void remove_samples(char *);

/* save.c prototypes */
void write_forest_file(FILE *,time_t);
int read_forest_file(FILE *);


