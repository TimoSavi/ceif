/* xmalloc.c -- malloc with out of memory checking 
   Copyright (C) 1990, 91, 92, 93, 94, 95, 96, 99 Free Software Foundation, Inc.
   Copyright (C) 2009 Timo Savinen

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  */

#include "ceif.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if STDC_HEADERS

#else  /* !STDC_HEADERS */

extern size_t strlen ();
extern char *strcpy ();

VOID *calloc ();
VOID *malloc ();
VOID *realloc ();
void free ();
#endif

static size_t total_allocation = 0;


static VOID *
fixup_null_alloc (n)
     size_t n;
{
  VOID *p;

  p = 0;
  if (n == 0)
    p = malloc ((size_t) 1);
  if (p == 0)
    {
      /* possible revisions: release some memory and re-try, print
	 more information (e.g. line number of input file) */
      panic("Memory exhausted",NULL,NULL);
      exit(1);
    }
  return p;
}

/* Allocate N bytes of memory dynamically, with error checking.  */

VOID *
xmalloc (n)
     size_t n;
{
  VOID *p;

  p = malloc (n);
  if (p == 0)
    p = fixup_null_alloc (n);
  total_allocation += n;
  return p;
}

/* Allocate memory for N elements of S bytes, with error checking.  */

VOID *
xcalloc (n, s)
     size_t n, s;
{
  VOID *p;

  p = calloc (n, s);
  if (p == 0)
    p = fixup_null_alloc (n);
  total_allocation += n * s;
  return p;
}

/* Change the size of an allocated block of memory P to N bytes,
   with error checking.
   If P is NULL, run xmalloc.  */

VOID *
xrealloc (VOID *p, size_t n)
{
  size_t old;

  if (p == 0)
    return xmalloc (n);
  old = sizeof(p);
  p = realloc (p, n);
  if (p == 0)
    p = fixup_null_alloc (n);
  total_allocation += n - old;
  return p;
}

/* Make a copy of a string in a newly allocated block of memory. */

char *
xstrdup (str)
     char *str;
{
  VOID *p;
  size_t len = strlen (str);

  p = xmalloc (len + 1);
  strcpy (p, str);
  total_allocation += len;
  return p;
}

FILE *
xfopen(char *name, char *mode, char bin_asc)
{
    register FILE *ret;
    static int stdin_opened = 0;
    static int stdout_opened = 0;

    if(name[0] == '-' && name[1] == '\000')
    {
        if( mode[0] == 'r' &&  mode[1] == '\000') 
        {
            if(stdin_opened) panic("stdin allready open",NULL,NULL);
            stdin_opened = 1;
            return stdin;
        } else 
        {
            if(stdout_opened) panic("stdout allready open",NULL,NULL);
            stdout_opened = 1;
            return stdout;
        }

   }


   ret = fopen(name,mode);
   if(ret == NULL) panic("Error in opening file",name,strerror(errno));

#if defined(HAVE_SETMODE) && defined(WIN32)
   if(bin_asc == 'a') setmode(fileno(ret),O_TEXT);
   if(bin_asc == 'b') setmode(fileno(ret),O_BINARY);
#endif
   return ret;
}

FILE *
xfopen_test(char *name, char *mode, char bin_asc)
{
   register FILE *ret;

   ret = fopen(name,mode);
   if(ret == NULL) return ret;

#if defined(HAVE_SETMODE) && defined(WIN32)
   if(bin_asc == 'a') setmode(fileno(ret),O_TEXT);
   if(bin_asc == 'b') setmode(fileno(ret),O_BINARY);
#endif
   return ret;
}

void
print_alloc_debug()
{
    DEBUG("Total dynamic memory allocation = %u\n",(unsigned int) total_allocation);
}
