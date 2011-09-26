## Process this file with automake to produce Makefile.in
#------------------------------------------------------------------------------#
# ool subdirectory
#
# Luis Alberto D'Afonseca
# since: Mar, 25, 2004
#------------------------------------------------------------------------------#

header-links: remove-links
	HEADERLIST="../*.h $(top_srcdir)/*/*.h"; \
	for h in $$HEADERLIST; do \
	  BASENAME=`basename $$h`; \
	  test -e $$BASENAME || $(LN_S) $$h $$BASENAME; \
	done

remove-links: 
	rm -f *.h

all: all-am header-links

clean:     clean-am     remove-links
distclean: distclean-am remove-links

CLEANFILES = *~

#------------------------------------------------------------------------------#