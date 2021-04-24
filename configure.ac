AC_INIT([ceif], [2.1], [tjsa@iki.fi])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC
AC_HEADER_STDC
AC_SEARCH_LIBS([cos],[m])
AC_CONFIG_HEADERS([config.h])
AC_CHECK_FUNCS([getopt_long])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_DEFINE([CEIF_CONFIG],["~/.ceifrc"],[Some variables can be set here])
AC_OUTPUT
