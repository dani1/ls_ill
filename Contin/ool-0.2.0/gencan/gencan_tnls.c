/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * gencan/gencan_tnls.c
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
 * $Id: gencan_tnls.c,v 1.8 2005/05/17 18:10:05 biloti Exp $
 *----------------------------------------------------------------------------*/

#include <gencan.h>

/*----------------------------------------------------------------------------*/
int gencan_tnls_extrapolation( ool_conmin_minimizer         *M, 
			       conmin_gencan_state          *st,
			       ool_conmin_gencan_parameters *P,
			       double                        alpha,
			       double                        fplus  );

int gencan_tnls_interpolation( ool_conmin_minimizer         *M, 
			       conmin_gencan_state          *st,
			       ool_conmin_gencan_parameters *P,
			       double alpha,
			       double fplus,
			       double gtd );

/*----------------------------------------------------------------------------*/
int gencan_tnls( ool_conmin_minimizer         *M, 
		 conmin_gencan_state          *st,
		 ool_conmin_gencan_parameters *P   )
{
  gsl_vector *X        = M->x;
  gsl_vector *gradient = M->gradient;
  gsl_vector *Xplus    = st->Xtrial;

  /* Direct access to vector data */
  double *x = X->data;
  double *g = gradient->data;
  double *d = st->D->data;
  double *xplus = Xplus->data;

  /* Constant values */
  const size_t nind = st->nind;

  /* Internal variables */
  double fplus;
  double gtd;
  double alpha;
  double gptd;

  /* Compute directional derivative */
  gtd = cblas_ddot( (int)nind, g, 1, d, 1 );

  /* Compute first trial */  
  alpha = GSL_MIN( 1.0, st->tnls_amax );

  /* Xplus = X + alpha D */
  conmin_vector_memcpy( nind, xplus, x );
  cblas_daxpy( (int)nind, alpha, d, 1, xplus, 1 );

  /* Evaluate objective function */
  fplus = conmin_calc_f( M, nind, st->Ind, Xplus, X );

  /* Test Armijo and beta-condition and decide for:
   * 1 - accepting the trial point,
   * 2 - interpolate or
   * 3 - extrapolate. */
  if( st->tnls_amax > 1.0 ) 
    {
      /* X+D belongs to the interior of the feasible set (amax > 1) */
    
      /* Verify Armijo */
      if( fplus <= M->f + P->gamma * alpha * gtd ) 
	{
	  /* Armijo condition holds */

	  /* Evaluate the gradient of objective function */
	  conmin_calc_g( M, nind, st->Ind, Xplus, X, gradient );

	  /* Eval gptd = < g, d > */
	  gptd = cblas_ddot( (int)nind, g, 1, d, 1 );

	  /* Verify directional derivative (beta condition) */
	  if ( gptd >= P->beta * gtd ) 
	    {
	      /* Step = 1 was ok, finish the line search */

	      M->f = fplus;
	      conmin_vector_memcpy( nind, x, xplus );

	      return OOL_SUCCESS;
	    }
	  else
	    {
	      return gencan_tnls_extrapolation( M, st, P, alpha, fplus );
	    }
	}
      else
	{
	  return gencan_tnls_interpolation(M, st, P, alpha, fplus, gtd);
	}
    }
  else 
    {
      /* x + d does not belong to the feasible set (amax <= 1) */
      if( fplus < M->f )
	{
	  return gencan_tnls_extrapolation( M, st, P, alpha, fplus );
	}
      else
	{
	  return gencan_tnls_interpolation(M, st, P, alpha, fplus, gtd);
	}
    }
}

/*----------------------------------------------------------------------------*/
int gencan_tnls_extrapolation( ool_conmin_minimizer         *M, 
			       conmin_gencan_state          *st,
			       ool_conmin_gencan_parameters *P,
			       double                        alpha,
			       double                        fplus  )
{
  gsl_vector *X        = M->x;
  gsl_vector *gradient = M->gradient;
  gsl_vector *Xplus    = st->Xtrial;

  /* Direct access to vector data */
  double *x = X->data;
  double *d = st->D->data;
  double *l = st->L->data;
  double *u = st->U->data;

  double *xplus = Xplus->data;
  double *xtemp = st->tnls_Xtemp->data;

  /* Constant values */
  const size_t nind = st->nind;

  /* Internal variables */
  double atemp;
  double ftemp;

  size_t ii, extrap;
  short same;

  /* Iterations */
  extrap = 0;
  do
    {
      extrap++;

      /* Test if maximum number of extrapolation was exceeded */
      if ( extrap > P->maxextrap )
	 {
	    M->f = fplus;
	    conmin_vector_memcpy( nind, x, xplus );
		 
	    if (extrap > 0 || st->tnls_amax < 1){
	       conmin_calc_g( M, nind, st->Ind, Xplus, X, gradient );
	    }
	    return GENCAN_TNLS_MAXEXTRAP;
	 }

      /* Chose new step */
      if( alpha < st->tnls_amax && st->tnls_amax < P->next*alpha ) 
	atemp = st->tnls_amax;
      else
	atemp = P->next * alpha;
  
      /* Compute new trial point. Xtemp = X + atemp*D */
      conmin_vector_memcpy( nind, xtemp, x );
      cblas_daxpy( (int)nind, atemp, d, 1, xtemp, 1 );

      /* Project */
      if( atemp > st->tnls_amax )
	conmin_vector_maxofmin( nind, xtemp, l, xtemp, u );

      /* Test if this is not the same point as the previous one.
       * This test is performed only when alpha >= amax. */
      if( alpha > st->tnls_amax ) 
	{
	   same = 1;
	   for (ii = 0; ii<nind && same; ii++)
	      {
		 double aux;
		 
		 aux = P->epsrel * fabs( xplus[ii] );
		 
		 if ( fabs( xtemp[ii] - xplus[ii] ) >
		      GSL_MAX( aux, P->epsabs ))
		    {
		       same = 0;
		    }
	      }

	   if (same)
	      {
		 /* Finish the extrapolation with the current point */
		 M->f = fplus;
		 
		 conmin_vector_memcpy( nind, x, xplus );
		 
		 if (extrap > 0 || st->tnls_amax < 1){
		    conmin_calc_g( M, nind, st->Ind, Xplus, X, gradient );
		 }
		 return OOL_SUCCESS;
	      }
	}

      ftemp = conmin_calc_f( M, nind, st->Ind, st->tnls_Xtemp, X );

      if( ftemp < fplus ) 
	{
	  /* If the functional value decreases then set the current
	   * point and continue the extrapolation */
	  
	  alpha = atemp;
	  fplus = ftemp;
	  conmin_vector_memcpy( nind, xplus, xtemp );

	  continue;
	}
      else 
	{
	  /* If the functional value does not decrease then discard the
	   * last trial and finish the extrapolation with the previous
	   * point */

	  M->f = fplus;

	  conmin_vector_memcpy( nind, x, xplus );
	  if (extrap > 0 || st->tnls_amax < 1){
	     conmin_calc_g( M, nind, st->Ind, X, X, gradient );
	  }

	  return OOL_SUCCESS;
	}
    }
  while( 1 );

  /* Just to make gcc happy */
  return OOL_SUCCESS;
}

/*----------------------------------------------------------------------------*/
int gencan_tnls_interpolation( ool_conmin_minimizer         *M, 
			       conmin_gencan_state          *st,
			       ool_conmin_gencan_parameters *P,
			       double alpha,
			       double fplus,
			       double gtd )
{
  gsl_vector *X        = M->x;
  gsl_vector *gradient = M->gradient;
  gsl_vector *Xplus    = st->Xtrial;

  /* Direct access to vector data */
  double *x = X->data;
  double *d = st->D->data;
  double *xplus = Xplus->data;

  /* Constant values */
  const size_t nind = st->nind;

  /* Internal variables */
  size_t interp;
  double atemp;

  /* Initialization */
  interp = 0;

  /* Iterations */
  do 
    {
      interp++;

      /* Test Armijo condition */
      if( fplus <= M->f + P->gamma * alpha * gtd ) 
	{
	  /* Finish the line search */
	  M->f = fplus;

	  /* X = Xplus */
	  conmin_vector_memcpy( nind, x, xplus );

	  /* Evaluate objective function gradient */
	  conmin_calc_g( M, nind, st->Ind, X, X, gradient );

	  return OOL_SUCCESS;
      }
   
      /* Compute new step */
      if( alpha < P->sigma1 ) 
	  alpha= alpha / P->nint;
      else
	{
	  /* quadratic model */
	  atemp = -gtd * alpha*alpha /
	           (2 * (fplus - M->f - alpha * gtd));

	  if( atemp < P->sigma1       || 
	      atemp > P->sigma2*alpha  ) alpha = alpha / P->nint;
	  else                           alpha = atemp;
	}

      /* Compute new trial point: xplus = x + alpha d */
      conmin_vector_memcpy( nind, xplus, x );
      cblas_daxpy( (int)nind, alpha, d, 1, xplus, 1 );

      /* Evaluate objective function */
      fplus = conmin_calc_f( M, nind, st->Ind, Xplus, X );

      /* Test whether at least mininterp interpolations were made
       * and the steplength is soo small */
      if ( gencan_are_close( nind, alpha, d, x,  P->epsrel, P->epsabs ) &&
	   interp > P->mininterp ){
	 return OOL_FLSEARCH;
      } 
  }
  while( 1 );

  /* Just to make gcc happy */
  return OOL_SUCCESS;
}

/*----------------------------------------------------------------------------*/
