## Process this file with automake to produce Makefile.in
#------------------------------------------------------------------------------#
# gencan subdirectory
#
# Luis Alberto D'Afonseca
# since: Apr, 04, 2004
#------------------------------------------------------------------------------#

noinst_LTLIBRARIES = liboolgencan.la

liboolgencan_la_SOURCES = ool_conmin_gencan.h \
                          gencan.h            \
                          gencan.c            \
                          gencan_tools.c      \
                          gencan_iterate.c    \
                          gencan_spgls.c      \
                          gencan_tnls.c       \
                          gencan_cg.c

pkginclude_HEADERS = ool_conmin_gencan.h

noinst_HEADERS = gencan.h	

#------------------------------------------------------------------------------#
check_PROGRAMS = test_gencan

TESTS = test_gencan

test_gencan_SOURCES = test_gencan.c
test_gencan_LDADD   = ../libool.la

CLEANFILES = *~

#------------------------------------------------------------------------------#