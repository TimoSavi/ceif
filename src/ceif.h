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
#define SAMPLES_MIN 32            // minimun number of samples for a forest
#define TREE_UPDATE_PERCENT 50    // Percentage of tree to rebuild when updating a forest

/* Data structures */
/* All forest related structures will be managed by dynamic tables
 */

struct sample
{
    double *dimension;      // dimension table, dynamically reserved
};

struct node
{
    int level;              // level count;
    int sample_count;       // number of samples
    double *n;               // random normal vector having dimensions count of coordinates
    double pdotn;           // calculate p dot n for performance issues
    int left;               // first left node, -1 if not existing
    int rigth;              // first rigth node, -1 if not existing
};

struct tree
{
    int node_count;         // number of nodes in n
    int node_cap;           // space reserved for node table
    struct node *n;         // table of nodes
    int first;              // index to first node in table.
};

struct forest
{
    char *category;         // Category string, all data having this category string in input data will behandled by this
    int X_count;	        // Number of dimensions
    int X_current;          // whis smaple should be taken next to tree 
    int X_cap;              // Memory allocated for X, in terms of units of struct sample
    double c;               // Average path length for the forest
    int heigth_limit;       // max tree depth
    struct sample *X;       // dimensions for this forest
    double *min;            // learn data min values dimension
    double *max;            // learn data max values dimension
    double *center;         // learn data center values dimension
    double *dim_density;    // learn data average attribute distance dimension
    struct tree *t;         // Tree table, NULL if not initialized
};


/* Global data */
extern int dim_idx[];
extern int category_idx[];
extern int dims_ignore[];       // Which dimensions are not processed
extern int dims_category[];       // dimensions to be used as category

extern int tree_count;                // trees / forest
extern int samples_max;              // max samples / tree
extern int samples_total;              // max samples / forest
extern int dimensions;            // dimensions in current setup
extern int categories;
extern int label_dim;             // which dimension is the lable dimension
extern char *print_string;
extern char input_separator;
extern int header;
extern double outlier_score;
extern double prange_extension_factor;
extern char *include_dims;
extern char *ignore_dims;
extern char *category_dims;

extern int forest_count;
extern int forest_cap;
extern struct forest *forest;           // forest table


/* ceif.c prototypes */
void panic(char *,char *,char *);
void parse_dims(char *,int,int *);


/* xmalloc.c prototypes */
VOID *xmalloc (size_t);
VOID *xcalloc (size_t, size_t);
VOID *xrealloc (VOID *, size_t);
char *xstrdup (char *);
FILE * xfopen(char *, char *, char);
FILE * xfopen_test(char *, char *, char);

/* file.c prototypes */
int parse_csv_line(char **,int,char *,char);

/* learn.c prototypes */
void train_forest(FILE *,int);
double parse_dim_attribute(char *);
double dot(double *, double *);
double c(int);
int dim_ok(int,int);
void add_to_X(struct forest *,char **, int ,int ,int);


/* analyze.c prototypes */
void analyze(FILE *, FILE *);
void categorize(FILE *, FILE *);
void init_dims(int);
char *make_category_string(int,char **);

/* save.c prototypes */
void write_forest_file(FILE *);
int read_forest_file(FILE *);


