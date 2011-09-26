dnl Process this file with autoconf to produce configure
dnl ---------------------------------------------------------------------------|
dnl Configuration of Open Optimization Library
dnl
dnl Luis Alberto D'Afonseca
dnl since: Feb, 27, 2004
dnl ---------------------------------------------------------------------------|

AC_INIT(conmin/conmin.c)

OOL_LT_VERSION="1:0:0" 
AC_SUBST(OOL_LT_VERSION)

AM_INIT_AUTOMAKE(ool,0.2.0)
AM_CONFIG_HEADER(config.h)

AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET

AC_LANG_C
AC_HEADER_STDC

dnl Check for inline functions
dnl ---------------------------------------------------------------------------|
AC_C_INLINE

dnl Cheking for libraries
dnl ---------------------------------------------------------------------------|
AC_CHECK_LIB(m,main,,       [AC_MSG_ERROR(ool requires math library)])
AC_CHECK_LIB(gslcblas,main,,[AC_MSG_ERROR(ool requires cblas library)])
AC_CHECK_LIB(gsl,main,,     [AC_MSG_ERROR(ool requires gsl)])

dnl ---------------------------------------------------------------------------|
AC_OUTPUT( ool_version.h    \
           Makefile         \
           ool/Makefile     \
           conmin/Makefile  \
           pgrad/Makefile   \
           spg/Makefile     \
           gencan/Makefile  \
           tools/Makefile   \
           doc/Makefile     )

dnl ---------------------------------------------------------------------------|
