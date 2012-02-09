%change -I_folder to include folders in which have been installed ool and
%gsl
mex -I/usr/local/include -lool -lgsl -lgslcblas -lm contin.c
