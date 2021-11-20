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
#include <time.h>

#ifdef HAVE_LIBFASTJSON_JSON_H

#include <libfastjson/json.h>

// json_object_to_file and json_object_from_file json-c compability  mappings might be missing for libfastjson
#ifndef json_object_to_file
#define json_object_to_file fjson_object_to_file
#endif

#ifndef json_object_from_file
#define json_object_from_file fjson_object_from_file
#endif

#else

#ifdef HAVE_JSON_JSON_H
#include <json/json.h>
#else 
#ifdef HAVE_JSON_C_JSON_H
#include <json-c/json.h>
#endif
#endif

#endif

#if defined HAVE_JSON_JSON_H || defined HAVE_LIBFASTJSON_JSON_H || HAVE_JSON_C_JSON_H   // headers are defined only if devel library is present
#define NVL(a) ((a) ? (a) : "")

#define GLOBALS "globals"
#define FORESTS "forests"
#define CATEGORY "category"
#define SAMPLE_COUNT "sampleCount"
#define LAST_UPDATED "lastUpdated"
#define SAMPLES "samples"

#define DIMENSIONS "dimensions"
#define FOREST_COUNT "forestCount"
#define PRINT_STRING "printString"
#define PRINTF_FORMAT "printfFormat"
#define TREE_COUNT "treeCount"
#define SAMPLES_MAX "samplesMax"
#define INPUT_SEPARATOR "inputSeparator"
#define LIST_SEPARATOR "listSeparator"
#define HEADER "header"
#define OUTLIER_SCORE "outlierScore"
#define CATEGORY_DIMS "categoryDims"
#define LABEL_DIMS "labelDims"
#define INCLUDE_DIMS "includeDims"
#define IGNORE_DIMS "ignoreDims"
#define SCORE_DIMS "scoreDims"
#define TEXT_DIMS "textDims"
#define FILTER_STR "filter"
#define DECIMALS "decimals"
#define UNIQUE_SAMPLES "uniqueSamples"
#define AGGREGATE "aggregate"

#if defined HAVE_JSON_C_SET_SERIALIZATION_DOUBLE_FORMAT || defined HAVE_JSON_OBJECT_NEW_DOUBLE_S   // in versions 0.15 and 0.12 or fastjson
    static char double_print[10];
#endif

/* makes global variable json object 
 */
static json_object *
write_globals(void)
{
    char str[] = "X";
    char scorestr[20];
    char *filter_str;

    filter_str = xstrdup(make_csv_line(cat_filter,cat_filter_count,';'));

    json_object *globals = json_object_new_object();
    json_object *jdims = json_object_new_int(dimensions);
    json_object *jforest_count = json_object_new_int(forest_count);
    json_object *jprint_string = json_object_new_string(print_string);
    json_object *jprintf_format = json_object_new_string(printf_format);
    json_object *jtree_count = json_object_new_int(tree_count);
    json_object *jsamples_max = json_object_new_int(samples_max);
    str[0] = input_separator;
    json_object *jinput_separator = json_object_new_string(str);
    str[0] = list_separator;
    json_object *jlist_separator = json_object_new_string(str);
    json_object *jheader = json_object_new_int(header);
    sprintf(scorestr,"%f%s",outlier_score,scale_score ? "s" : (percentage_score ? "%" : ""));
    json_object *joutlier_score = json_object_new_string(scorestr);
    json_object *jcategory_dims = json_object_new_string(category_dims);
    json_object *jlabel_dims = json_object_new_string(label_dims);
    json_object *jscore_dims = json_object_new_string(score_dims);
    json_object *jignore_dims = json_object_new_string(ignore_dims);
    json_object *jinclude_dims = json_object_new_string(include_dims);
    json_object *jtext_dims = json_object_new_string(text_dims);
    json_object *jfilter_str = json_object_new_string(filter_str);
    json_object *jdecimals = json_object_new_int(decimals);
    json_object *junique_samples = json_object_new_int(unique_samples);
    json_object *jaggregate = json_object_new_int(aggregate);

    json_object_object_add(globals,DIMENSIONS,jdims);
    json_object_object_add(globals,FOREST_COUNT,jforest_count);
    json_object_object_add(globals,PRINT_STRING,jprint_string);
    json_object_object_add(globals,PRINTF_FORMAT,jprintf_format);
    json_object_object_add(globals,TREE_COUNT,jtree_count);
    json_object_object_add(globals,SAMPLES_MAX,jsamples_max);
    json_object_object_add(globals,INPUT_SEPARATOR,jinput_separator);
    json_object_object_add(globals,LIST_SEPARATOR,jlist_separator);
    json_object_object_add(globals,HEADER,jheader);
    json_object_object_add(globals,OUTLIER_SCORE,joutlier_score);
    json_object_object_add(globals,CATEGORY_DIMS,jcategory_dims);
    json_object_object_add(globals,LABEL_DIMS,jlabel_dims);
    json_object_object_add(globals,INCLUDE_DIMS,jinclude_dims);
    json_object_object_add(globals,IGNORE_DIMS,jignore_dims);
    json_object_object_add(globals,TEXT_DIMS,jtext_dims);
    json_object_object_add(globals,SCORE_DIMS,jscore_dims);
    json_object_object_add(globals,FILTER_STR,jfilter_str);
    json_object_object_add(globals,DECIMALS,jdecimals);
    json_object_object_add(globals,UNIQUE_SAMPLES,junique_samples);
    json_object_object_add(globals,AGGREGATE,jaggregate);

    return globals;
}

/* make a json object for a forest, includes samples too
 */
static json_object *
write_forest(int forest_idx)
{
    int i,j;
    struct forest *f = &forest[forest_idx];

    json_object *jforest = json_object_new_object();
    json_object *jsamples = json_object_new_array();
    json_object *jsample;
    json_object *jdouble;
    json_object *jcategory = json_object_new_string(NVL(f->category));
    json_object *jsample_count = json_object_new_int(f->X_count);
    json_object *jlast_updated = json_object_new_int64(f->last_updated);   // Might work after 19 January 2038...
    
    json_object_object_add(jforest,CATEGORY,jcategory);
    json_object_object_add(jforest,SAMPLE_COUNT,jsample_count);
    json_object_object_add(jforest,LAST_UPDATED,jlast_updated);

    for(i = 0;i < f->X_count;i++)
    {
        jsample = json_object_new_array();
        for(j = 0;j < dimensions;j++)
        {
#if !defined HAVE_JSON_C_SET_SERIALIZATION_DOUBLE_FORMAT && defined HAVE_JSON_OBJECT_NEW_DOUBLE_S 
            char ds[1024];
            sprintf(ds,double_print,f->X[i].dimension[j]);
            jdouble = json_object_new_double_s(f->X[i].dimension[j],ds);
#else
            jdouble = json_object_new_double(f->X[i].dimension[j]);
#endif
            json_object_array_add(jsample,jdouble);
        }
        json_object_array_add(jsamples,jsample);
    }

    json_object_object_add(jforest,SAMPLES,jsamples);

    return jforest;
}



/* write forest data to json file
 * returns true if json is supported
 */

int write_forest_file_json(char *file_name,time_t delete_interval)
{
    int i;
    time_t now = time(NULL);

#if defined HAVE_JSON_C_SET_SERIALIZATION_DOUBLE_FORMAT || defined HAVE_JSON_OBJECT_NEW_DOUBLE_S   // in versions 0.15 and 0.12 or libfastjson
    sprintf(double_print,"%%.%df",decimals);
#ifdef HAVE_JSON_C_SET_SERIALIZATION_DOUBLE_FORMAT
    json_c_set_serialization_double_format(double_print,JSON_C_OPTION_GLOBAL);
#endif
#endif

    json_object *root = json_object_new_object();
    json_object *jglobals;
    json_object *jforests = json_object_new_array();
    json_object *jforest;

    jglobals = write_globals();
    json_object_object_add(root,GLOBALS,jglobals);

    for(i = 0;i < forest_count;i++)
    {
        if(delete_interval == (time_t) 0 || (delete_interval > (time_t) 0 && forest[i].last_updated >= now - delete_interval)) 
        {
            jforest = write_forest(i);
            json_object_array_add(jforests,jforest);
        }
    }

    json_object_object_add(root,FORESTS,jforests);

    if(json_object_to_file(file_name,root) == -1) panic("Cannot write to file",file_name,"");
        
    json_object_put(root);
    return 1;
}

/* read global variables
*/
static
void read_globals(json_object *globals)
{
    char str[10240];
    char *f[FILTER_MAX];
    int i,c;

    json_object *jdims  ;
    json_object *jforest_count  ;
    json_object *jprint_string  ;
    json_object *jprintf_format  ;
    json_object *jtree_count  ;
    json_object *jsamples_max  ;
    json_object *jinput_separator  ;
    json_object *jlist_separator  ;
    json_object *jheader  ;
    json_object *joutlier_score  ;
    json_object *jcategory_dims  ;
    json_object *jlabel_dims  ;
    json_object *jscore_dims  ;
    json_object *jignore_dims  ;
    json_object *jinclude_dims  ;
    json_object *jtext_dims  ;
    json_object *jfilter_str  ;
    json_object *jdecimals  ;
    json_object *junique_samples  ;
    json_object *jaggregate  ;

    if(!json_object_object_get_ex(globals,DIMENSIONS,&jdims)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,FOREST_COUNT,&jforest_count)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,PRINT_STRING,&jprint_string)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,PRINTF_FORMAT,&jprintf_format)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,TREE_COUNT,&jtree_count)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,SAMPLES_MAX,&jsamples_max)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,INPUT_SEPARATOR,&jinput_separator)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,LIST_SEPARATOR,&jlist_separator)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,HEADER,&jheader)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,OUTLIER_SCORE,&joutlier_score)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,CATEGORY_DIMS,&jcategory_dims)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,LABEL_DIMS,&jlabel_dims)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,SCORE_DIMS,&jscore_dims)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,IGNORE_DIMS,&jignore_dims)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,INCLUDE_DIMS,&jinclude_dims)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,TEXT_DIMS,&jtext_dims)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,FILTER_STR,&jfilter_str)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,DECIMALS,&jdecimals)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,UNIQUE_SAMPLES,&junique_samples)) panic("Error in globals object","","");
    if(!json_object_object_get_ex(globals,AGGREGATE,&jaggregate)) panic("Error in globals object","","");

    dimensions = json_object_get_int(jdims);
    if(dimensions > DIM_MAX) dimensions = DIM_MAX;

    forest_count = json_object_get_int(jforest_count);
    print_string = xstrdup(json_object_get_string(jprint_string));
    printf_format = xstrdup(json_object_get_string(jprintf_format));
    tree_count = json_object_get_int(jtree_count);
    samples_max = json_object_get_int(jsamples_max);
    strcpy(str,json_object_get_string(jinput_separator));
    input_separator = str[0];
    strcpy(str,json_object_get_string(jlist_separator));
    list_separator = str[0];
    header = json_object_get_int(jheader);

    strcpy(str,json_object_get_string(joutlier_score));
    parse_user_score(str);

    strcpy(str,json_object_get_string(jcategory_dims));
    category_dims = xstrdup(str);
    category_idx_count = parse_dims(str,category_idx);

    strcpy(str,json_object_get_string(jlabel_dims));
    label_dims = xstrdup(str);
    label_idx_count = parse_dims(str,label_idx);

    strcpy(str,json_object_get_string(jscore_dims));
    score_dims = xstrdup(str);
    score_idx_count = parse_dims(str,score_idx);

    strcpy(str,json_object_get_string(jignore_dims));
    ignore_dims = xstrdup(str);
    ignore_idx_count = parse_dims(str,ignore_idx);
    
    strcpy(str,json_object_get_string(jinclude_dims));
    include_dims = xstrdup(str);
    include_idx_count = parse_dims(str,include_idx);

    strcpy(str,json_object_get_string(jtext_dims));
    text_dims = xstrdup(str);
    text_idx_count = parse_dims(str,text_idx);

    strcpy(str,json_object_get_string(jfilter_str));
    c = parse_csv_line(f,FILTER_MAX,str,';');
    for(i = 0;i < c;i++) add_category_filter(f[i]);

    decimals = json_object_get_int(jdecimals);
    unique_samples = json_object_get_int(junique_samples);
    aggregate = json_object_get_int(jaggregate);

    samples_total = max_total_samples ?  max_total_samples : tree_count * samples_max;  
}

/* read forest data including samples
 */
static
void init_forest(int forest_idx,json_object *jforest)
{
    struct forest *f;
    int i,j,sample_count,attr_count;
    static double new[DIM_MAX];

    json_object *jsamples;
    json_object *jsample;
    json_object *jcategory;
    json_object *jlast_updated;
    
    f = &forest[forest_idx];
    
    if(!json_object_object_get_ex(jforest,CATEGORY,&jcategory)) panic("Missing category string","","");
    f->category = xstrdup(json_object_get_string(jcategory));

    if(!json_object_object_get_ex(jforest,LAST_UPDATED,&jlast_updated)) panic("Missing last updated date","","");;
    f->last_updated = (time_t) json_object_get_int64(jlast_updated);

    f->c = 0;
    f->heigth_limit = 0;
    f->X = NULL;
    f->X_count = 0;
    f->X_cap = 0;
    f->X_summary = -1;
    f->t = NULL;
    f->min = NULL;
    f->max = NULL;
    f->scale_range_idx = -1;
    f->avg = NULL;
    f->summary = NULL;
    f->dim_density = NULL;
    f->analyzed = 0;
    f->filter = 0;
    f->total_rows = 0;
    f->analyzed_rows = 0;
    f->high_analyzed_rows = 0;
    f->extra_rows = 0;
    f->percentage_score = 0.0;
    f->min_score = 1.0;
    f->test_average_score = 0.0;
    
    add_forest_hash(forest_idx,f->category);

    f->X_count = 0;

    if(json_object_object_get_ex(jforest,SAMPLES,&jsamples))
    {
        sample_count = json_object_array_length(jsamples);   // get number of samples from array
        f->X_cap = sample_count + 1;
        f->X = xmalloc(f->X_cap * sizeof(struct sample));

        for(i = 0;i < sample_count;i++)
        {
            jsample = json_object_array_get_idx(jsamples,i);
            attr_count = json_object_array_length(jsample);

            for(j = 0;j < attr_count;j++) new[j] = json_object_get_double(json_object_array_get_idx(jsample,j));
            for(j = attr_count;j < dimensions;j++)  new[j] = 0.0;     // reset rest, if array is shorter than expected, should not happen

            add_to_X(f,new,attr_count,1);
        }
    } else
    {
        f->X_cap = f->X_count + 1;
        f->X = xmalloc(f->X_cap * sizeof(struct sample));
    }
}


/* read json formatted forest data into memory
 */
int read_forest_file_json(char *file_name)
{
    int forest_idx;

    json_object *root = json_object_from_file(file_name);
    json_object *globals;
    json_object *forests;

    if(!root) panic("Cannot read file: ",file_name,"");

    if(json_object_object_get_ex(root,GLOBALS,&globals))
    {
        read_globals(globals);
    } else
    {
        panic("No globals in JSON file","","");
    }
    
    if(json_object_object_get_ex(root,FORESTS,&forests))
    {
        forest_count = json_object_array_length(forests);         // get actual forest count from array

        forest_cap = forest_count + 1;
        forest = xmalloc(forest_cap * sizeof(struct forest));     // allocate forest table

        for(forest_idx = 0;forest_idx < forest_count;forest_idx++)
        {
            init_forest(forest_idx,json_object_array_get_idx(forests,forest_idx));
        }
    } else
    {
        forest_count = 0;
        forest_cap = forest_count + 1;
        forest = xmalloc(forest_cap * sizeof(struct forest));     // allocate forest table
    }

    json_object_put(root);

    return 1;
}

#else

/* no json-c support cases
 *  */
int
write_forest_file_json(char *file_name,time_t delete_interval)
{
    return 0;
}

int read_forest_file_json(char *file_name)
{
    panic("JSON not implemented","","");
    return 0;
}

#endif
