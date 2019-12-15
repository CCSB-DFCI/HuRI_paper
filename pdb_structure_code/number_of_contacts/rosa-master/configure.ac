#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT([rosa], [1.0], [bug-report@librosa.org])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR([src/rosa/atom.cpp])
AC_CONFIG_HEADERS([config.h])
m4_pattern_allow([AM_INIT_AUTOMAKE])
AM_INIT_AUTOMAKE([-Wall -Werror -Wno-extra-portability foreign])

# Marenostrum optimization
AC_ARG_ENABLE( marenostrum,
               [AS_HELP_STRING([--enable-marenostrum], [customize for compilation on Marenostrum])],
               [AS_IF( [test "x$enable_marenostrum" != xno], [mymarenostrum=true] ) ],
               [mymarenostrum=false] )

# Sets default language to C++
AC_LANG(C++)

# Checks for programs
AC_PROG_CXX
AC_PROG_LEX
AC_PROG_YACC
AC_PROG_RANLIB
AC_PROG_INSTALL

# Checks for program paths
AC_PATH_PROGS(BASH, bash sh)


# Checks for conditionals

# Debug informations
AC_ARG_ENABLE( debug,
               [AS_HELP_STRING([--enable-debug], [include debug information])],
               [AS_IF( [test "x$enable_debug" != xno], [mydebug=true] ) ] )

AS_IF(
  [test "x$mydebug" != xtrue],
    [AC_SUBST( CXXFLAGS, ["-O3"] )] ,
    [AC_SUBST( CXXFLAGS, ["-Wall -ansi -g -O2"] )
     AC_DEFINE([ACTIVATE_DEBUG], [1],
               [Activates inclusion of debug information.])] )

# Static executables
AC_ARG_ENABLE( static,
               [AS_HELP_STRING([--enable-static], [compile static executables])],
               [AS_IF( [test "x$enable_static" != xno], [mystatic=true] ) ] )

AS_IF(
  [test "x$mystatic" == xtrue],
    [AC_SUBST( LDFLAGS, ["$LDFLAGS -static"] )] )

# Exceptions
AC_ARG_ENABLE( exceptions,
               [AS_HELP_STRING([--enable-exceptions], [include exceptions handling code])],
               [AS_IF( [test "x$enable_exceptions" != xno], [myexceptions=true] ) ],
               [myexceptions=true] )

AS_IF(
  [test "x$myexceptions" != xtrue],
    [AC_SUBST( CXXFLAGS, ["$CXXFLAGS -fno-exceptions"] )] )

AM_CONDITIONAL([INCLUDE_DEBUG], [test x$mydebug = xtrue])

# 32 bits
#AC_ARG_ENABLE( 32bits,
#               [AS_HELP_STRING([--enable-32bits], [compile for 32 bits machines])],
#               [AS_IF( [test "x$enable_32bits" != xno], [my32bits=true] ) ] )

#AS_IF(
#  [test "x$my32bits" != xtrue],
#    [AC_SUBST( CXXFLAGS, ["$CXXFLAGS -m32"] )] )

#AS_IF(
#  [test "x$my32bits" == xtrue],
#    [AC_SUBST( LDFLAGS, ["$LDFLAGS -m32"] )] )

# Aggressive optimizations
AC_ARG_ENABLE( optimization,
               [AS_HELP_STRING([--enable-optimization], [include extreme optimization options])],
               [AS_IF( [test "x$enable_optimization" != xno], [myoptimization=true] ) ] )

AS_IF(
  [test "x$myoptimization" = xtrue],
    [AS_IF( [test "x$CXX" == xg++],
      [ACX_CHECK_CXX_FLAGS( [-ffast-math], fast_math_opt, [AC_SUBST( CXXFLAGS, ["$CXXFLAGS -ffast-math"] )] )
       ACX_CHECK_CXX_FLAGS( [-funroll-all-loops], unroll_all_loops, [AC_SUBST( CXXFLAGS, ["$CXXFLAGS -funroll-all-loops"] )] )
       ACX_CHECK_CXX_FLAGS( [-ftree-vectorize], tree_vectorize, [AC_SUBST( CXXFLAGS, ["$CXXFLAGS -ftree-vectorize"] )] )
       AS_IF( [test "x$mymarenostrum" == xtrue], [AC_SUBST( CXXFLAGS, ["$CXXFLAGS -mcpu=970 -mtune=970"] )] )],
       AS_IF( [test "x$CXX$mymarenostrum" == xxlCtrue], [AC_SUBST( CXXFLAGS, ["$CXXFLAGS -qstrict -qtune=ppc970 -qcache=auto"])]))])

AS_IF(
  [test "x$myoptimization" = xtrue],
    [AS_IF( [test "x$CXX" == xllvm-g++],
      [AC_SUBST( CXXFLAGS, ["$CXXFLAGS -ffast-math -funroll-all-loops -ftree-vectorize"] )] ) ])

# Boost libraries
AC_ARG_WITH([boost],
            [AS_HELP_STRING([--with-boost],
              [enable use of Boost libraries instead of tr1 extensions to the standard library])],
            [],
            [with_boost=no])

LIBBOOST=
AS_IF([test "x$with_boost" != xno],
  [AC_CHECK_HEADER([boost/shared_ptr.hpp],
    [AC_DEFINE([HAVE_LIBBOOST], [1],
               [Define if you have Boost libraries])
    ],
    [AC_MSG_FAILURE(
       [--with-boost was given, but test for boost failed])])])

# Additional flags

# Checks for libraries

# Checks for header files
AC_CHECK_HEADER(getopt.h)
AC_CHECK_HEADER(memory)


# Checks for typedefs, structures, and compiler characteristics
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions
AC_CHECK_FUNCS(strncmp)
AC_CHECK_FUNCS(getopt_long)
AC_CHECK_FUNCS(exit)

# In case I want to build the script using configure directly. In this case
# there is no need of the EXTRA_DIST variable in test/Makefile.am. A drawback
# is the fact that variables as $(datadir) are not fully expanded and may
# refer to other directory variables.
# AC_CONFIG_FILES([test/test_atomsel.sh], [chmod +x test/test_atomsel.sh])

AC_OUTPUT( Makefile             \
           src/Makefile         \
           src/rosa/Makefile    \
           src/tools/Makefile    \
           test/Makefile        \
         )

AC_OUTPUT
