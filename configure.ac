AC_INIT([ceif],[2.5],[tjsa@iki.fi])
AC_USE_SYSTEM_EXTENSIONS
AM_INIT_AUTOMAKE([-Wall -Werror foreign tar-pax])
AC_PROG_CC
AC_HEADER_STDC
AC_SEARCH_LIBS([cos],[m])

have_fastjson=yes

AC_SEARCH_LIBS([fjson_object_new_object],[fastjson],[],[have_fastjson=no
  AC_MSG_WARN([Is fastjson-c development library installed? JSON may not be supported])])

if test "x${have_fastjson}" = xno; then
     AC_SEARCH_LIBS([json_object_new_object],[json-c],[],[AC_MSG_WARN([Is json-c development library installed? JSON may not be supported])])
fi

AC_CHECK_HEADERS([libfastjson/json.h json/json.h json-c/json.h])

AC_CONFIG_HEADERS([config.h])
AC_CHECK_FUNCS([getopt_long json_c_set_serialization_double_format json_object_new_double_s fjson_object_new_double_s])
AC_CHECK_DECLS([HAVE_FJSON_OBJECT_NEW_DOUBLE_S],[AC_DEFINE([HAVE_JSON_OBJECT_NEW_DOUBLE_S],[1],[Make a single symbol for json_object_new_double_s])])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_DEFINE([CEIF_CONFIG],["~/.ceifrc"],[Some variables can be set here])
AC_OUTPUT
