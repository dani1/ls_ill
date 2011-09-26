/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * gencan/gencan_cg.c
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
 * $Id: gencan_cg.c,v 1.13 2005/05/13 22:08:00 biloti Exp $
 *----------------------------------------------------------------------------*/

#include<gencan.h>
#include<gsl/gsl_poly.h>

/*----------------------------------------------------------------------------*/
int gencan_cg( ool_conmin_minimizer         *M, 
	       conmin_gencan_state          *st,
	       ool_conmin_gencan_parameters *P   )
{
  gsl_vector *X        = M->x;
  gsl_vector *gradient = M->gradient;

  /* Direct access to vector data */
  double *x = X->data;
  double *g = gradient->data;

  double *l = st->L->data;
  double *u = st->U->data;

  double *w = st->cg_W->data;
  double *r = st->cg_R->data;
  double *d = st->cg_D->data;

  double *sprev = st->cg_Sprev->data;

  /* To keep the same variable internal name */
  double *s = st->D->data;

  /* Internal variables */
  const size_t nind = st->nind;

  size_t ii;
  size_t itnqmp;

  double q, qprev, qamax, qamaxn;
  double currprog, bestprog;
  double gts, dtr, dtw, dts;
  double snorm2, snorm2_prev;
  double rnorm2, rnorm2_prev;
  double dnorm2, gnorm2;
  double cg_eps, cg_eps2;
  double amax, amax1, amax1n, amax2, amax2n, amaxn;
  double alpha;
  double aux;

  int    iter;
  int    cg_maxit;

  /* Pointers to avoid index arithmetics */
  double *sii, *xii, *lii, *uii, *dii;

  /* cg_maxit */
  if ( P->ucgmia < 0 || P->ucgmib < 0 )
     {
	if ( P->nearlyq ) cg_maxit = nind;
	else
	   {
	      double kappa;
	      double aux;
	      if ( P->cg_scre == 1 )
		 kappa = ( log10( st->gpeucn2 / st->gpeucn20 ) / 
			   log10( st->epsgpen2 / st->gpeucn20 ));
	      else
		 kappa = ( log10( st->gpsupn / st->gpsupn0 ) /
			   log10( P->epsgpsn / st->gpsupn0 ));

	      aux = GSL_MIN( 1, kappa );
	      kappa = GSL_MAX( 0, aux );

	      aux = GSL_MIN( nind, 10*log10((double)nind) );
	      cg_maxit = (int)( (1-kappa) * GSL_MAX( 1, aux ) + 
				kappa * nind );
	   }
     }
  else{
     int iaux = (int) ( P->ucgmia * nind + P->ucgmib );
     cg_maxit = GSL_MAX( 1, iaux );
  }
  
  
  if ( P->cg_scre == 1 ){
     cg_eps = sqrt( pow( 10,
			 st->acgeps*log10( st->gpeucn2 ) + st->bcgeps
			 ));
  }
  else{
     cg_eps = pow (10, ( st->acgeps * log10( st->gpsupn ) + st->bcgeps ));
  }

  aux = GSL_MIN( P->cg_epsi, cg_eps );
  cg_eps = GSL_MAX( P->cg_epsf, aux );
  
  cg_eps2 = gsl_pow_2( cg_eps );

  gnorm2 = cblas_ddot( (int)nind, g, 1, g, 1 );

  conmin_vector_set_zero( nind, s );
  conmin_vector_memcpy( nind, r, g );  

  q      = 0.0;
  gts    = 0.0;
  snorm2 = 0.0;
  rnorm2 = gnorm2;

  iter = 0;
  itnqmp = 0;
  bestprog = 0;

  /* repeat while |r|_2 = |H s + g|_2 > eps * |g|_2
   *--------------------------------------------------------------------------*/
  while( rnorm2 > cg_eps2 * gnorm2 ) 
    {
      /* Stopping criteria */
      if( iter > cg_maxit ) return OOL_FINNERIT;

      /* Compute direction */ 
      if( iter == 0 ) 
	{
	  conmin_vector_memcpy( nind, d, r );
	  cblas_dscal( (int)nind, -1.0, d, 1 );

	  dnorm2 =  rnorm2;
	  dtr    = -rnorm2;
	}
      else 
	{
	  double beta;

	  beta = rnorm2 / rnorm2_prev;

	  /* d = -r + beta d */
	  cblas_dscal( (int)nind, beta, d, 1 );
	  cblas_daxpy( (int)nind, -1.0, r, 1, d, 1 );

	  aux    = dtr + alpha * dtw;
	  dnorm2 =  rnorm2 + beta * ( beta*dnorm2 - 2.0*aux );
	  dtr    = -rnorm2 + beta * aux;
	}

      /* Force d to be a descent direction of q(s), i.e.,
       * <\nabla q(s), d> = <H s + g, d> = <r, d> \le 0 */
      if( dtr > 0.0 ) 
	{
	  cblas_dscal( (int)nind, -1.0, d, 1 );
	  dtr = -dtr;
	}

      /* w = H d */
      /* Why is not the macro used here ? */
      conmin_calc_Hv( M, nind, st->Ind, X, X, st->cg_D, st->cg_W );

      /* dtw = < d, w > */
      dtw = cblas_ddot( (int)nind, d, 1, w, 1 );

      /*======================================================================
	void gencan_cg_maximum_step( ... )
	======================================================================*/

      /* Compute maximum step 
       *----------------------------------------------------------------------*/
  
      /* amax1 > 0 and amax1n < 0 are the values of alpha such that 
       * ||s + alpha * d||_2 or ||s + alpha * d||_\infty = delta */

      /* dts = < d, s > */
      dts = cblas_ddot( (int)nind, d, 1, s, 1 );

      if( P->trtype == 0 ) /* 2-norm trust-region radius */
	{
	  double aa = dnorm2;
	  double bb = 2 * dts;
	  double cc = snorm2 - st->cg_delta * st->cg_delta;
	  
	  gsl_poly_solve_quadratic (aa, bb, cc, &amax1n, &amax1);

	}
      else  /* Infinite-norm trust-region radius */
	{
	   
	  sii = s;
	  dii = d;

	  amax1  =  P->infabs;
	  amax1n = -P->infabs;

	  for( ii = 0; ii < nind; ii++ ) 
	    {
	      if( *dii > 0.0 ) 
		{
		  amax1  = GSL_MIN( amax1,  ( st->cg_delta - *sii) / *dii );
		  amax1n = GSL_MAX( amax1n, (-st->cg_delta - *sii) / *dii );
		}
	      else if( *dii < 0.0 ) 
		{
		  amax1  = GSL_MIN( amax1,  (-st->cg_delta - *sii) / *dii );
		  amax1n = GSL_MAX( amax1n, ( st->cg_delta - *sii) / *dii );
		}

	      sii++;
	      dii++;
	    }
	  
	}

      /* amax2 > 0 and amax2n < 0 are the maximum and the minimum values
       * of alpha such that l - x <= s + alpha * d <= u - x, respectively */
  
      amax2  =  P->infabs;
      amax2n = -P->infabs;

      uii = u;
      lii = l;
      xii = x;
      dii = d;
      sii = s;
      for( ii = 0; ii < nind; ii++ ) 
	{
	  double aux1 = (*uii - *xii - *sii) / *dii;
	  double aux2 = (*lii - *xii - *sii) / *dii;

	  if ( *dii > 0.0 ) 
	     {
		amax2  = GSL_MIN( amax2,  aux1 );
		amax2n = GSL_MAX( amax2n, aux2 );
	     }
	  else if( *dii < 0.0 ) 
	    {
	       amax2  = GSL_MIN( amax2,  aux2 );
	       amax2n = GSL_MAX( amax2n, aux1 );
	    }
	  
	  uii++;
	  lii++;
	  dii++;
	  xii++;
	  sii++;
	}

      /* Compute amax as the minimum among amax1 and amax2, and amaxn as the
       * minimum among amax1n and amax2n. Moreover change amaxn by -amaxn to
       * have amax and amaxn as maximum steps along d direction (and not -d
       * in the case of amaxn) */
      amax  = GSL_MIN( amax1,  amax2  );
      amaxn = GSL_MAX( amax1n, amax2n );

      /*======================================================================
	======================================================================*/

      /* Compute the step, 
       * and the quadratic functional value at the new point */
      qprev = q;

      /* If d^T H d > 0 then take the conjugate gradients step */
      if( dtw > 0.0 ) 
	{
	  alpha = GSL_MIN( amax, rnorm2 / dtw );
	  q     = q + alpha * ( alpha * dtw / 2 + dtr );
	}
      else 
	{
	  /* If d^T H d <= 0 and function f is nearly quadratic then 
	     take the point with the minimum functional value (q) among 
	     the actual one and the ones which are at the boundary, i.e., 
	     the best one between q(s), q(s + amax*d) and q(s + amaxn*d) */
	
	  qamax= q + amax * ( amax * dtw / 2 + dtr );
	  
	  /* If we are at iteration zero then take the maximum positive
	     step in the minus gradient direction */
    
	  if ( iter == 0 )
	    {
	      alpha = amax;
	      q     = qamax;
	    }
	  else{

	     if ( P->nearlyq && ( qamax < q || qamaxn < q ) ) 
		{
		   /* If we are not in the first iteration then if
		      function f is nearly quadratic and q(s + amax *
		      d) or q(s + amaxn * d) is smaller than q(s), go
		      to the best point in the boundary */

		   qamaxn = q + amaxn * ( 0.5 * amaxn * dtw + dtr );

		   if( qamax < qamaxn ) 
		      {
			 alpha = amax;
			 q     = qamax;
		      }
		   else 
		      {
			 alpha = amaxn;
			 q     = qamaxn;
		      }
		}
	     else 
		{
		   /* Else, stop at the current point */
		   return OOL_FDDIR;
		}
	  }
	}

      /* Compute new s */
      conmin_vector_memcpy ( nind, sprev, s );
      cblas_daxpy( (int)nind, alpha, d, 1, s, 1 );

      snorm2_prev = snorm2;
      snorm2      = snorm2 + alpha *( alpha * dnorm2 + 2.0 * dts );

      /* Compute the residual r = H s + g  */
      rnorm2_prev = rnorm2;

      cblas_daxpy( (int)nind, alpha, w, 1, r, 1 );
      rnorm2 = cblas_ddot( (int)nind, r, 1, r, 1 );

      /* Increment number of iterations */
      iter++;

      /* Test other stopping criteria */
  
      /* Test angle condition */
      gts = cblas_ddot( (int)nind, g, 1, s, 1 );

      aux = gsl_pow_2( P->theta ) * gnorm2 * snorm2;

      if ( gts > 0.0 || gsl_pow_2( gts ) < aux ) 
	{
	  conmin_vector_memcpy( nind, s, sprev );

	  snorm2 = snorm2_prev;
	  q      = qprev;
	  
	  return GENCAN_CG_S;
	}

      /* If we are in the boundary of the trust region then stop */
      /* if ( alpha == amax1 || alpha == amax1n )  */
      /*    return GENCAN_CG_TRUST_REGION;         */

      if ( fabs( alpha - amax1 ) < P->epsabs )
	{
	  alpha = amax1;
	  return GENCAN_CG_TRUST_REGION;
	}

      if ( fabs( alpha - amax1n ) < P->epsabs )
	{
	  alpha = amax1n;
	  return GENCAN_CG_TRUST_REGION;
	}

      /* If we are in the boundary of the box also stop */  
      /* if ( alpha == amax2 || alpha == amax2n )  */
      /*    return GENCAN_CG_BOUNDARY;             */

      if ( fabs( alpha - amax2 ) < P->epsabs )
	{
	  alpha = amax2;
	  return GENCAN_CG_BOUNDARY; 
	}

      if ( fabs( alpha - amax2n ) < P->epsabs )
	{
	  alpha = amax2n;
	  return GENCAN_CG_BOUNDARY; 
	}

      if ( gencan_are_close( nind, alpha, d, s, P->epsrel, P->epsabs )){
	 return GENCAN_CG_CLOSE_ITERATE;
      }

      currprog = qprev - q;
      bestprog = GSL_MAX( currprog, bestprog );
      
      if ( currprog < P->cg_epsnqmp * bestprog ){
	 itnqmp++;
	 
	 if ( itnqmp >= P->cg_maxitnqmp ){
	    
	    return GENCAN_CG_INSUF_PROG;
	 }
      }
      else itnqmp = 0;

    } /* main loop */
  

  return OOL_SUCCESS;
}

/*----------------------------------------------------------------------------*/
