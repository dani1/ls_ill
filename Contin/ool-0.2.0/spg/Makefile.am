## Process this file with automake to produce Makefile.in
#------------------------------------------------------------------------------#
# spg subdirectory
#
# Ricardo Biloti
# since: May 3rd, 2005
#------------------------------------------------------------------------------#

noinst_LTLIBRARIES = liboolspg.la

liboolspg_la_SOURCES = ool_conmin_spg.h \
                          spg.h            \
                          spg.c

pkginclude_HEADERS = ool_conmin_spg.h

noinst_HEADERS = spg.h	

#------------------------------------------------------------------------------#
check_PROGRAMS = test_spg

TESTS = test_spg

test_spg_SOURCES = test_spg.c
test_spg_LDADD   = ../libool.la

CLEANFILES = *~

#------------------------------------------------------------------------------#
