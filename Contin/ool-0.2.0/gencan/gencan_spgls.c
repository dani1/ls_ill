/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * gencan/gencan_spgls.c
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
 * Ricardo Biloti
 * since: Mar 28th, 2004
 *
 * $Id: gencan_spgls.c,v 1.9 2005/05/17 18:10:05 biloti Exp $
 *----------------------------------------------------------------------------*/

#include<gencan.h>

/*----------------------------------------------------------------------------*/
int gencan_spgls( ool_conmin_minimizer         *M, 
		  conmin_gencan_state          *st,
		  ool_conmin_gencan_parameters *P   )
{
  gsl_vector *X        = M->x;
  gsl_vector *gradient = M->gradient;

  size_t nn = X->size;

  /* Direct access to vector data */
  double *l = st->L->data;
  double *u = st->U->data;
  double *d = st->D->data;
  double *x = X->data;

  double *xtrial = st->Xtrial->data;

  /* Internal variables */
  size_t interp;
  size_t imax;

  double alpha;
  double dinfn;
  double gtd;
  double ftrial;

  /* Compute first trial point, spectral projected gradient direction,
   * and directional derivative <g,d> */
  alpha = 1;

  /* Xtrial = min{ U, max[ L, ( X-lambda G ) ] } */
  gsl_vector_memcpy( st->Xtrial, X );
  gsl_blas_daxpy( -(st->lambda), gradient, st->Xtrial );
  conmin_vector_minofmax( st->n, xtrial, u, l, xtrial );

  /* D = Xtrial - X */
  gsl_vector_memcpy( st->D, st->Xtrial );
  gsl_vector_sub( st->D, X );

  /* Inifite norm of D and < G, D > */
  imax  = gsl_blas_idamax( st->D );
  dinfn = fabs( gsl_vector_get( st->D, imax ) );
  gsl_blas_ddot( gradient, st->D, &gtd );

  /* Evaluate objective function */
  OOL_CONMIN_EVAL_F( M, st->Xtrial, ftrial ); 

  interp = 0;

  /* While Armijo isn't satizfied repeat */
  while( ftrial > M->f + P->gamma*alpha*gtd ) 
  {
     /* Test if the trial point has a function value lower than fmin */
     if( ftrial < M->f ) 
	{
	   M->f = ftrial;
	   gsl_vector_memcpy( X, st->Xtrial );

	   return OOL_UNBOUNDEDF;
	}

     interp++;
     
     if( alpha < P->sigma1 )
	{
	   alpha /= P->nint;
	}
     else 
	{
	   /* quadratic model */
	   double atemp = ( -gtd*alpha*alpha )/( 2.0*(ftrial-M->f-alpha*gtd) );
	   
	   if( atemp < P->sigma1 || 
	       atemp > P->sigma2*alpha ) alpha /= P->nint;
	   else                          alpha  = atemp;
	}
     
     /* Compute new trial point 
      * Xtrial = X + alpha D */
     gsl_vector_memcpy( st->Xtrial, X );
     gsl_blas_daxpy( alpha, st->D, st->Xtrial );
     
     /* Evaluate objective function */
     OOL_CONMIN_EVAL_F( M, st->Xtrial, ftrial ); 

     /* Test whether at least mininterp interpolations were made 
      * and two consecutive iterates are close enough */
     if( interp > P->mininterp && 
	 gencan_are_close( nn, alpha, d, x, P->epsrel, P->epsabs )
	 ) 
      {
	 M->f = ftrial;
	 gsl_vector_memcpy( X, st->Xtrial );
	 
	 return OOL_FLSEARCH;
      }
  }

  /* Set new values of f and X */
  M->f = ftrial;
  gsl_vector_memcpy( X, st->Xtrial );

  return OOL_SUCCESS;
}

/*----------------------------------------------------------------------------*/
