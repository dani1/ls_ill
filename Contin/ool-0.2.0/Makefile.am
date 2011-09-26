## Process this file with automake to produce Makefile.in
#------------------------------------------------------------------------------#
# Top Level
#
# Luis Alberto D'Afonseca
# since: Feb, 28, 2004
#------------------------------------------------------------------------------#


SUBDIRS = ool conmin pgrad spg gencan tools doc

SUBLIBS = conmin/liboolconmin.la  \
          pgrad/liboolpgrad.la \
          spg/liboolspg.la \
          gencan/liboolgencan.la  \
          tools/libooltools.la

pkginclude_HEADERS = ool_version.h 

lib_LTLIBRARIES   = libool.la
libool_la_SOURCES = version.c
libool_la_LIBADD  = $(SUBLIBS)
libool_la_LDFLAGS = -version-info $(OOL_LT_VERSION)

CLEANFILES = *~

#------------------------------------------------------------------------------#
