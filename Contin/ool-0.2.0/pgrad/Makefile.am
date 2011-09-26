## Process this file with automake to produce Makefile.in
#------------------------------------------------------------------------------#
# pgrad subdirectory
#
# Iara da Cunha
# since: July, 07, 2004
#------------------------------------------------------------------------------#

noinst_LTLIBRARIES = liboolpgrad.la

liboolpgrad_la_SOURCES = ool_conmin_pgrad.h \
                          pgrad.h            \
                          pgrad.c

pkginclude_HEADERS = ool_conmin_pgrad.h

noinst_HEADERS = pgrad.h	

#------------------------------------------------------------------------------#
check_PROGRAMS = test_pgrad

TESTS = test_pgrad

test_pgrad_SOURCES = test_pgrad.c
test_pgrad_LDADD   = ../libool.la

CLEANFILES = *~

#------------------------------------------------------------------------------#