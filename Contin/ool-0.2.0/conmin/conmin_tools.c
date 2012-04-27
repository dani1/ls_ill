/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * conmin/conmin_tools.c
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
 * since: Mar, 28, 2004
 *
 * $Id: conmin_tools.c,v 1.5 2005/04/15 21:18:52 biloti Exp $
 *----------------------------------------------------------------------------*/

#include <conmin_tools.h>
#include <conmin_vector.h>

/* Shrink vector V from the full space (n) to the reduced space (nind).
 *----------------------------------------------------------------------------*/
void conmin_shrink( const size_t nind, gsl_vector_uint *Ind, gsl_vector *V )
{
  size_t ii, indii;
 
  for( ii = 0; ii < nind; ii++ )
    {
      indii = gsl_vector_uint_get( Ind, ii );

      gsl_vector_swap_elements( V, ii, indii );
    }
}

/* Expands vector v from the reduced space (nind) to the full space (n).
 *----------------------------------------------------------------------------*/
void conmin_expand( const size_t nind, gsl_vector_uint *Ind, gsl_vector *V )
{
  size_t jj, ii, indii;

  ii = nind;
  for( jj = 0; jj < nind; jj++ ) 
    {
      ii--;
      indii = gsl_vector_uint_get( Ind, ii );
      gsl_vector_swap_elements( V, ii, indii );
    }
} 

/* Evaluate objective function from the reduced space
 *----------------------------------------------------------------------------*/
double conmin_calc_f( ool_conmin_minimizer *M, 
		      const size_t          nind, 
		      gsl_vector_uint      *Ind,
		      gsl_vector           *X,
		      gsl_vector           *Xc )
{
  const size_t missing = X->size - nind;
  double f;

  if( missing > 0 )
    {
      double *x  = &(X ->data[nind]);
      double *xc = &(Xc->data[nind]);

      /* Complete X with values from Xc */
      conmin_vector_memcpy( missing, x, xc );

      /* Expand to full space*/
      conmin_expand( nind, Ind, X );
    }
      
  /* Compute f */
  OOL_CONMIN_EVAL_F( M, X, f );

  if( missing > 0 )
    {
      /* Shrink to reduced space */
      conmin_shrink( nind, Ind, X );
    }

  return f;
}

/* Evaluate gradient from the reduced space
 *----------------------------------------------------------------------------*/
void conmin_calc_g( ool_conmin_minimizer *M, 
		    const size_t          nind, 
		    gsl_vector_uint      *Ind,
		    gsl_vector           *X,
		    gsl_vector           *Xc,
		    gsl_vector           *G  )
{
  const size_t missing = X->size - nind;

  if( missing > 0 )
    {
      double *x  = &(X ->data[nind]);
      double *xc = &(Xc->data[nind]);

      /* Complete X with values from Xc */
      conmin_vector_memcpy( missing, x, xc );

      /* Expand to full space*/
      conmin_expand( nind, Ind, X );
    }

  /* Compute gradient */
  OOL_CONMIN_EVAL_DF( M, X, G );

  if( missing > 0 )
    {
      /* Shrink to reduced space */
      conmin_shrink( nind, Ind, X );
      conmin_shrink( nind, Ind, G );
    }
}

/* Evaluate hessian times a vector from the reduced space
 *----------------------------------------------------------------------------*/
void conmin_calc_Hv( ool_conmin_minimizer *M, 
		     const size_t          nind, 
		     gsl_vector_uint      *Ind,
		     gsl_vector           *X,
		     gsl_vector           *Xc,
		     gsl_vector           *V,
		     gsl_vector           *Hv   )
{
  const size_t missing = X->size - nind;

  if( missing > 0 )
    {
      double *x  = &(X ->data[nind]);
      double *v  = &(V ->data[nind]);
      double *xc = &(Xc->data[nind]);

      /* Complete X with values from Xc */
      conmin_vector_memcpy( missing, x, xc );

      /* Complete V with zeros */
      conmin_vector_set_zero( missing, v );

      /* Expand to full space*/
      conmin_expand( nind, Ind, X );
      conmin_expand( nind, Ind, V );
    }

  /* Compute gradient */
  OOL_CONMIN_EVAL_HV( M, X, V, Hv );

  if( missing > 0 )
    {
      /* Shrink to reduced space */
      conmin_shrink( nind, Ind, X  );
      conmin_shrink( nind, Ind, V  );
      conmin_shrink( nind, Ind, Hv );
    }
}

/*----------------------------------------------------------------------------*/
