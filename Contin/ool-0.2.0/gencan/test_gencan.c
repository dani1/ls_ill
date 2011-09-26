/*----------------------------------------------------------------------------*
 * Open Optimization Library - Example Programs
 * 
 * gencan/test_gencan.c
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * Sergio Drumond Ventura
 * Luis Alberto D'Afonseca
 * since: Feb, 16, 2004
 *
 * $Id: test_gencan.c,v 1.8 2005/05/11 01:09:09 biloti Exp $
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <ool/ool_conmin.h>

/*----------------------------------------------------------------------------*/
double quadratic    ( const gsl_vector *X, void *p  );
void   quadratic_df ( const gsl_vector *X, void *p, gsl_vector *G );
void   quadratic_fdf( const gsl_vector *X, void *p, double     *f,           
		      gsl_vector *G );
void   quadratic_Hv ( const gsl_vector *X, void *p, const gsl_vector *V, 
		      gsl_vector *hv );

void iteration_echo( ool_conmin_minimizer *M );

/*----------------------------------------------------------------------------*/
int main( void )
{
  size_t nn   = 100; 
  size_t nmax = 1000;

  size_t ii;
  int status;

  const ool_conmin_minimizer_type *T = ool_conmin_minimizer_gencan;
  ool_conmin_gencan_parameters P; 

  ool_conmin_function   F;
  ool_conmin_constraint C;
  ool_conmin_minimizer *M;

  gsl_vector *X;

  /* Starting objective function */
  F.n   = nn; 
  F.f   = &quadratic;
  F.df  = &quadratic_df;
  F.fdf = &quadratic_fdf;
  F.Hv  = &quadratic_Hv;

  /* Contraints */
  C.n = F.n;
  C.L = gsl_vector_alloc( C.n );
  C.U = gsl_vector_alloc( C.n );

  gsl_vector_set_all( C.L, -3.0 );
  gsl_vector_set_all( C.U,  3.0 );

  /* Starting Point */
  X = gsl_vector_alloc( F.n );
  for(ii = 0; ii < nn; ii++ )
     gsl_vector_set( X, ii, 1.0+ii );
  
  /* Allocation of minimizer */
  M = ool_conmin_minimizer_alloc( T, nn ); 
  
  /* Starting minimizer parameters */
  ool_conmin_parameters_default( T, (void*)(&P) );

  /* Starting minimizer */
  ool_conmin_minimizer_set( M, &F, &C, X, (void*)(&P) );
  
  /* Iterating */
  ii = 0;
  status = OOL_CONTINUE;
  
  printf( "%4i : ", ii );
  iteration_echo ( M );

  while( ii < nmax && status == OOL_CONTINUE ){
     ii++;
     
     /* New iteration */
     status = ool_conmin_minimizer_iterate( M );
     status = ool_conmin_is_optimal( M );
     
     printf( "%4i : ", ii );
     iteration_echo( M );
     
  }

  if(status == OOL_SUCCESS){
     printf("Minimum found: \n");
     iteration_echo( M );
  }
     
  printf("\nfunction evalutaions = %i\ngradient evaluations = %i\n",
	 ool_conmin_minimizer_fcount( M ),
	 ool_conmin_minimizer_gcount( M ));

  /* Freeing memory */
  gsl_vector_free( C.L );
  gsl_vector_free( C.U );
  gsl_vector_free( X );
  
  ool_conmin_minimizer_free( M );
  
  return OOL_SUCCESS;
}

/*----------------------------------------------------------------------------*/
double 
quadratic( const gsl_vector *X, void* params )
{
  size_t ii, nn;
  double f, x;

  nn = X->size;

  f = 0;

  for( ii = 0; ii < nn; ii++ )
    {
      x  = gsl_vector_get( X, ii ) - (ii+1)/10.0;
      f += (ii+1) * gsl_pow_2( x );
    }

  return f;
}

/*----------------------------------------------------------------------------*/
void 
quadratic_df( const gsl_vector *X, void* params, gsl_vector *G )
{
  size_t ii, nn;
  double xx, gg;

  nn = X->size;

  for( ii = 0; ii < nn; ii++ )
    {
      xx = gsl_vector_get( X, ii ) - (ii+1)/10.0;
      gg = 2.0 * (ii+1) * xx;
      gsl_vector_set( G, ii, gg );
    }
}

/*----------------------------------------------------------------------------*/
void 
quadratic_fdf( const gsl_vector *X, void* params,
	       double *f, gsl_vector *G )
{
  size_t ii, ip1, nn;
  double xx;

  nn = X->size;
  (*f) = 0;

  for( ii = 0; ii < nn; ii++ )
    {
      ip1 = ii + 1;

      xx = gsl_vector_get( X, ii ) - (ii+1)/10.0;

      (*f) += ip1 * gsl_pow_2( xx );

      gsl_vector_set( G, ii, 2.0 * ip1 * xx );
    }
}
 
/*----------------------------------------------------------------------------*/
void 
quadratic_Hv( const gsl_vector *X, void *params,
	      const gsl_vector *V, gsl_vector *hv )
{
  size_t ii, nn;
  double aux;

  nn = X->size;

  for( ii = 0; ii < nn; ii++ )
    {
      aux = 2.0 * (ii+1) * gsl_vector_get( V, ii );
      gsl_vector_set( hv, ii, aux );
    }
}

/*----------------------------------------------------------------------------*/
void iteration_echo( ool_conmin_minimizer *M )
{
  gsl_vector *X = M->x;

  size_t ii, nn;

  nn = X->size;

  printf( "f( " );

  for( ii = 0; ii < 4; ii++ )
    printf( "%+6.3e, ", gsl_vector_get( X, ii ) );

  printf( "... ) = %+7.4e ||g|| = %6.3e\n",
	  ool_conmin_minimizer_minimum( M ),
	  ool_conmin_minimizer_size( M ));
}

/*----------------------------------------------------------------------------*/
