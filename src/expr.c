/*
 *    ceif - categorized extended isolation forest
 *
 *    Copyright (C) 2024 Timo Savinen
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
#include "tinyexpr.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#define FORMULA_MAX 100            // Maximum nmber of formulas
#define DATA_REFERENCE '$'         // data references have notation $n, where n is the data index number (starts from 1)
#define DECIMAL_SEPARATOR ':'      // expression can have trailer :n, where n is the number of decimals to be used when converting double to string

struct data_value_formula formula[FORMULA_MAX];
int formulas = 0;        // Number of formulas


/* parse input data reference 
 * reference has syntax $n where n is an integer 1-MAX_DIM
 * writes it to struct expression_data_reference
 *
 * s should start with $
 *
 * if ref->data_idx == -1 then an error has occured and reference is not valid
 */
static
void parse_data_reference(char *s, struct expression_data_reference *ref)
{
    char *num,save,*position;

    ref->data_idx = -1;

    if(*s != DATA_REFERENCE) return;

    position = s;
    
    s++;

    num = s;        // Number start

    while(isdigit(*s)) s++;

    if(num == s) return;  // No number after $

    save = *s;
    *s = '\000';
    ref->data_idx = atoi(num);
    *s = save;

    if(ref->data_idx < 1) 
    {
        ref->data_idx = -1;
        return;
    }

    ref->data_idx--;     // use as array index
    ref->length = s - position;
}




/* parse left and right parts of a formula
 * write left value (datavalue index) to left
 * write expression after = to right
 * if left is -1 or rigth null, the formula syntax is wrong
 */

static
void parse_left_right(char *fstring,int *left,char **right)
{
    char *c = fstring;
    struct expression_data_reference l;

    *left = -1;
    *right = NULL;

    do
    {
        switch(*c)
        {
            case DATA_REFERENCE:
                if(*left == -1)
                {
                    parse_data_reference(c,&l);
                    if(l.data_idx == -1) return;
                    *left = l.data_idx;
                    c+=l.length;
                } else
                {
                    return;
                }
                break;
            case '=':
                if(*left != -1)
                {
                    c++;
                    while(isspace(*c)) c++;
                    *right = c;
                    return;
                } else
                {
                    return;
                }
                break;
            default:
                if(!isspace(*c)) return;
                while(isspace(*c)) c++;
                break;
        }
    } while(*c);
}

/* parse decimals, expression can have trailer :n, where n is the number of 
 * decimals to be used when converting double to string
 *
 * : and everything after that is removed
 *
 * returns the number of decimals, default is six
 */
static 
int parse_decimals(char *expr)
{
    char *colon_pos;
    int decimals = 6;

    colon_pos = strchr(expr,DECIMAL_SEPARATOR);

    if(colon_pos != NULL)
    {
        *colon_pos = '\000';
        colon_pos++;
        decimals = atoi(colon_pos);
    }
    return decimals;
}


/* Remove one formula from the formula array
 */
static 
void remove_formula(char *remove)
{
    int i,j;

    for(i = 0;i < formulas;i++)
    {
        if(strcmp(formula[i].formula,remove) == 0)
        {
            formulas--;
            free(formula[i].formula);
            free(formula[i].expression);
            for(j = i;j < formulas;j++) formula[j] = formula[j + 1];
        }
    }
}





/* Parse one formula which is used to replace input data value
 * Formula is assumed to have following syntax:
 * $<ti> = <expression>
 * ti = target data value index (starts from 1)
 * expression = tinyexpr expression having $n references to input data. $n is the nth input data value 
 *
 * If formula starts with -, the formula is removed.
 */
void
parse_expression(char *new)
{
    struct data_value_formula *dvf;
    int value_idx;
    char *expr,*s;

    if(new[0] == '-')
    {
        remove_formula(&new[1]);
        return;
    }

    if(formulas == FORMULA_MAX) return;    // Silently ignore rest of formulas

    parse_left_right(new,&value_idx,&expr);

    if(value_idx == -1 || expr == NULL || *expr == '\000') panic("Error in data value expression ",new,NULL);

    dvf = &formula[formulas];
    
    formulas++;

    dvf->formula = xstrdup(new);
    dvf->expression = xstrdup(expr);
    dvf->decimals = parse_decimals(dvf->expression);
    dvf->target_data_idx = value_idx;
    dvf->dref_count = 0;

    s = dvf->expression;

    /* Add all data references to dvf->dref, even if it is not a valid reference */
    while(*s)
    {
        switch(*s)
        {
            case DATA_REFERENCE:
                if(dvf->dref_count < EXPRESSION_DATA_REFERENCE_MAX)
                {
                    parse_data_reference(s,&dvf->dref[dvf->dref_count]);
                    dvf->dref_count++;
                }
                break;
        }
        s++;
    }
}

/* check if a string is a valid floating point number in ascii representation
 * strtod is used for this.
 * If end points to \000 it is assumed that string is valid float
 * returns true if string is a valid float
 */
static
int check_float_string(char *string)
{
    char *end;

    if(*string == '\000') return 0;

    (void) strtod(string,&end);

    return *end == '\000';
}

/* Check if a input field (value index in data_idx) should be replaced using an expression
 * Expressions are only evaluated for numeric fields, text fields remain intact
 * If field is to be modified then it is evaluated using te_interp
 * If te_interp return NaN an error is raised
 */
char *
evaluate_data_expression(int data_idx, int value_count,char **values)
{
    int f,r;
    struct data_value_formula *dvf;
    char *s,*t;
    double newval;
    static char expr[2048];
    static char retval[50];

    // likely case, check this first
    if(!formulas) return values[data_idx];

    // No formulas to change data or data is not a float
    if(data_idx >= value_count) return "";
    if(!check_float_string(values[data_idx])) return values[data_idx];


    for(f = 0;f < formulas;f++)
    {
        if(formula[f].target_data_idx == data_idx)
        {
            dvf = &formula[f];
            s = dvf->expression;
            r = 0;
            t = expr;
            expr[0] = '\000';

            // copy expression to expr, replace all $-references with actual data values
            while(*s)
            {
                switch(*s)
                {
                    case DATA_REFERENCE:
                        if(dvf->dref[r].data_idx > -1 && dvf->dref[r].data_idx < value_count && r < dvf->dref_count)
                        {
                            strcpy(t,values[dvf->dref[r].data_idx]);
                            while(*t) t++;
                            s += dvf->dref[r].length;
                        } else
                        {
                            s++;
                        }
                        r++;
                        break;
                    default:
                        *t = *s;
                        t++;
                        s++;
                        break;
                }
            }
            *t = '\000';
           
            newval = te_interp(expr,0);
            if(isnormal(newval) || newval == 0.0)
            {
                sprintf(retval,"%.*f",dvf->decimals,newval);
            } else
            {

                info("Expression cannot be interpreted or evaluated ",dvf->expression,NULL);

                if(ignore_expression_errors)
                {
                    info("Expression with parameters expanded, this will be replace by zero",expr,NULL);
                    sprintf(retval,"%.*f",dvf->decimals,0.0);
                } else
                {
                    panic("Expression with parameters expanded",expr,NULL);
                }
            }
            return retval;
        }
    }
    return values[data_idx];   // default is to return the original data
}

