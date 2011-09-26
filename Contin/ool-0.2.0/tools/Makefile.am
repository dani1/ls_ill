## Process this file with automake to produce Makefile.in
#------------------------------------------------------------------------------#
# tools subdirectory
#
# Luis Alberto D'Afonseca
# since: May, 23, 2004
#------------------------------------------------------------------------------#

noinst_LTLIBRARIES = libooltools.la

libooltools_la_SOURCES = ool_tools_diff.h tools_diff.c

pkginclude_HEADERS = ool_tools_diff.h

#------------------------------------------------------------------------------#
check_PROGRAMS = test_tools

TESTS = test_tools

test_tools_SOURCES = test_tools.c
test_tools_LDADD   = ../libool.la

CLEANFILES = *~

#------------------------------------------------------------------------------#
