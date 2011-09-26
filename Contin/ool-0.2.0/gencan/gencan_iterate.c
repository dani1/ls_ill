/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * gencan/gencan_iterate.h
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
 * $Id: gencan_iterate.c,v 1.11 2005/05/15 14:10:56 biloti Exp $
 *----------------------------------------------------------------------------*/

#include<gencan.h>

/* Truncated Newton maximum step
 *----------------------------------------------------------------------------*/
static double 
gencan_tnls_maximum_step( ool_conmin_minimizer         *M, 
			  conmin_gencan_state          *st, 
			  ool_conmin_gencan_parameters *P  )
{
  /* Direct access to vector data */
  double *x =  M->x->data;
  double *l = st->L->data;
  double *u = st->U->data;
  double *d = st->D->data;

  double step = P->infabs;
  size_t ii;
      
  for( ii = 0; ii < st->nind; ii++ ) 
    {
      if( *d > 0 )
	{
	   double aux = ( *u - *x ) / *d;
	   step = GSL_MIN( step, aux );
	}
      else if( *d < 0 )
	{
	   double aux = ( *l - *x ) / *d;
	   step = GSL_MIN( step, aux );
	}

      x++;
      l++;
      u++;
      d++;
    }

  return step;
}

/* Spectral Step Length
 *----------------------------------------------------------------------------*/
static void
gencan_spg_steplength( conmin_gencan_state          *st, 
		       ool_conmin_gencan_parameters *P  )
{
  if( st->sty <= 0.0 ) 
    {
      st->lambda = GSL_MAX( 1.0, st->xeucn ) / sqrt( st->gpeucn2 );
    }
  else
    {
      double aux;
      double ss   = st->sts / st->sty;
		
      aux = GSL_MAX( P->lspgmi, ss );
      st->lambda = GSL_MIN( P->lspgma, aux );
    }
}

/* Iterating
 *----------------------------------------------------------------------------*/
int 
gencan_actual_iterate( ool_conmin_minimizer         *M, 
		       conmin_gencan_state          *st,
		       ool_conmin_gencan_parameters *P )
{
  /* Direct access to vector data */
  double *x =  M->x->data;
  double *l = st->L->data;
  double *u = st->U->data;
  /*  double *d = st->D->data; */

  /* Status of internal iterations */
  int lsflag;
  int cgflag;

  /* Internal variables */
  size_t ii, imax;

  /* Saving previous values */
  gsl_vector_memcpy( st->S, M->x        );
  gsl_vector_memcpy( st->Y, M->gradient );

  /* The actual iteration 
   *--------------------------------------------------------------------------*/
  if ( st->gieucn2 <= st->ometa2 * st->gpeucn2 ) 
    /* Compute the new iterate using an SPG iteration */
    {
      /* Perform a line search with spectral continuous projected gradient */
      lsflag = gencan_spgls( M, st, P );

      /* Compute the gradient for the new iterate */
      OOL_CONMIN_EVAL_DF( M, M->x, M->gradient );
    }
  else 
    /* The new iterate will belong to the closure of the actual face */
    {
      /* Shrink the point, its gradient and the bounds */
      conmin_shrink( st->nind, st->Ind, M->x        );
      conmin_shrink( st->nind, st->Ind, M->gradient );
      conmin_shrink( st->nind, st->Ind, st->L       );
      conmin_shrink( st->nind, st->Ind, st->U       );
    
      /* Compute the descent direction solving the newtonian system */
      cgflag = gencan_cg( M, st, P );

      /* Compute maximum step for truncated newton line search */
      if ( cgflag == GENCAN_CG_BOUNDARY )
	  st->tnls_amax = 1.0;
      else
	  st->tnls_amax = gencan_tnls_maximum_step( M, st, P );

      /* Perform the line search */
      lsflag = gencan_tnls( M, st, P );

      /* Expand the point, its gradient and the bounds */
      conmin_expand( st->nind, st->Ind, M->x        );
      conmin_expand( st->nind, st->Ind, M->gradient );
      conmin_expand( st->nind, st->Ind, st->L       );
      conmin_expand( st->nind, st->Ind, st->U       );

      /* If the line search in the Truncated Newton direction stopped due to
       * a very small step discard this iteration and force a SPG iteration */
      if ( lsflag == OOL_FLSEARCH )
	{
	  /* Perform a line search with spectral projected gradient */
	  lsflag = gencan_spgls( M, st, P );

	  /* Compute the gradient for the new iterate */
	  OOL_CONMIN_EVAL_DF( M, M->x, M->gradient );
	}
    }

  /* Prepare for the next iteration 
   *--------------------------------------------------------------------------*/
  /* Adjustment */
  for( ii = 0; ii < st->n; ii++ )
     {
	/* In principle, the current point could be slightly changed
	 * here, requiring a new function and gradient
	 * evaluation. But, according to the algorithms authors, this
	 * is done just to account for points that are "numerically¨
	 * at faces already. Thus, no additional evaluations are
	 * performed. (May 11th, 2005).
	 */
	if     ( x[ii] <= st->near_l[ii] ) x[ii] = l[ii];
	else if( x[ii] >= st->near_u[ii] ) x[ii] = u[ii];
     }

  /* Compute infinite and Euclidian-norm of X */
  imax      = gsl_blas_idamax( M->x );
  st->xsupn = fabs( gsl_vector_get( M->x, imax ) );
  st->xeucn = gsl_blas_dnrm2 ( M->x ); 

  /* Until now S = X_prev, now S = X - X_prev 
   * Compute s = x_{k+1} - x_k = X - S
   * and     y = g_{k+1} - g_k = G - Y  */
  gsl_vector_sub  ( st->S, M->x        ); /* S = S - X */
  gsl_vector_scale( st->S, -1.0        ); /* S = -S = X - S_prev */
  gsl_vector_sub  ( st->Y, M->gradient ); /* Y = Y - G */
  gsl_vector_scale( st->Y, -1.0        ); /* Y = -Y = G - Y_prev */

  /* Compute sts = s dot s 
   *         sty = s dot y
   * and     sinf = |s|_inf */
  gsl_blas_ddot( st->S, st->S, &(st->sts) );    
  gsl_blas_ddot( st->S, st->Y, &(st->sty) );
  imax      = gsl_blas_idamax( st->S );
  st->sinf  = fabs( gsl_vector_get( st->S, imax ) );

  /* Compute continuous project gradient */
  gencan_projected_gradient( st, M->x, M->gradient );

  /* Update spectral steplength */
  gencan_spg_steplength( st, P );

  /* Update trust-region radius */
  if ( P->trtype ) st->cg_delta = GSL_MAX( P->delmin, 10*sqrt( st->sts ) );
  else             st->cg_delta = GSL_MAX( P->delmin, 10*    ( st->sinf) );

  return lsflag;
}

/*----------------------------------------------------------------------------*/
