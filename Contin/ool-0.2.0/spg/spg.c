/*----------------------------------------------------------------------------*
 * Open Optimization Library - Spectral Projected Gradient Method
 * 
 * spg/spg.c
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
 * Ricardo Biloti
 * since: May 3rd, 2005
 *
 * $Id: spg.c,v 1.4 2005/05/10 20:24:27 biloti Exp $
 *----------------------------------------------------------------------------*/
#include <spg.h>

#define def_P ool_conmin_spg_parameters*  P = \
             (ool_conmin_spg_parameters*) M->parameters

#define def_st conmin_spg_state*  st = \
              (conmin_spg_state*) M->state

/* Set default values for spg parameters
 *----------------------------------------------------------------------------*/
static void
spg_parameters_default( void  *vpar )
{
   ool_conmin_spg_parameters *P = (ool_conmin_spg_parameters*) vpar;

   /* Setting default values of parameters */
   P->fmin     = -1.0e+99;
   P->tol      =  1.0e-4;
   P->alphamin =  1.0e-30;
   P->alphamax =  1.0e+30;
   P->gamma    =  1.0e-4; /* Armijo */
   P->sigma1   =  0.1;
   P->sigma2   =  0.9;
   P->M        =  10;
}

/* Set values for spg parameters
 *----------------------------------------------------------------------------*/
static int
spg_parameters_set( ool_conmin_minimizer *M, void *new_param )
{
   /* Copy new state values */
   memcpy( M->parameters, new_param, sizeof(ool_conmin_spg_parameters) );

   return OOL_SUCCESS;
}

/* Allocation  
 *----------------------------------------------------------------------------*/
static int
spg_alloc( void *vstate, size_t n )
{
   conmin_spg_state *state = (conmin_spg_state*)vstate;

   /* Space dimension */
   state->n = n;

   state->xx = gsl_vector_alloc( n );
   state->d  = gsl_vector_alloc( n );
   state->s  = gsl_vector_alloc( n );
   state->y  = gsl_vector_alloc( n );
   state->L  = gsl_vector_alloc( n );
   state->U  = gsl_vector_alloc( n );
  
   if( state->xx == NULL || 
       state->s  == NULL || 
       state->d  == NULL || 
       state->y  == NULL || 
       state->L  == NULL || 
       state->U  == NULL  )
      {
	 gsl_vector_free( state->xx );
	 gsl_vector_free( state->s  );
	 gsl_vector_free( state->d  );
	 gsl_vector_free( state->y  );
	 gsl_vector_free( state->L  );
	 gsl_vector_free( state->U  );

	 return OOL_ENOMEM;
      }

   /* This space can only be allocated after knowing parameter M */
   state->f = NULL;
  
   return OOL_SUCCESS;
}

/* Free memory allocated by spg_alloc
 *----------------------------------------------------------------------------*/
static void
spg_free( void *vstate )
{
   size_t ii;
   conmin_spg_state *state = (conmin_spg_state*)vstate;

   free( state->f );

   gsl_vector_free( state->xx );
   gsl_vector_free( state->s  );
   gsl_vector_free( state->d  );
   gsl_vector_free( state->y  );
   gsl_vector_free( state->L  );
   gsl_vector_free( state->U  );

}

/* Setting
 *----------------------------------------------------------------------------*/
static int
spg_set( ool_conmin_minimizer *M )
{
   def_st;
   CBLAS_INDEX ii;
   size_t jj, kk;

   /* Copy boundary vectors */
   gsl_vector_memcpy( st->L, M->con->L );
   gsl_vector_memcpy( st->U, M->con->U );

   /* Turns M->x into a feasible point */
   spg_proj( st->L, st->U, M->x );
   OOL_CONMIN_EVAL_FDF( M, M->x, &(M->f), M->gradient );
   
   /* Infinite norm of d <- g1 = [P(x - g) - x] */
   gsl_vector_memcpy( st->d, M->x );
   gsl_blas_daxpy( -1, M->gradient, st->d );
   spg_proj( st->L, st->U, st->d );
   gsl_vector_sub( st->d, M->x );

   M->size = fabs(gsl_vector_get( st->d, gsl_blas_idamax( st->d )));

   /* alpha_0 */
   ii = gsl_blas_idamax( M->gradient );
   st->alpha = 1/fabs(gsl_vector_get( M->gradient, ii ));

   /* Allocate space for g of previous iterations */
   st->M = ((ool_conmin_spg_parameters *) (M->parameters))->M;
   st->f = (double *) malloc( sizeof(double) * st->M );

   if ( st->f == NULL ){
      return OOL_ENOMEM;
   }

   st->m = 1;
   st->tail = 0;
   st->f[st->tail] = M->f;

   return OOL_SUCCESS;
}

/* Restarting
 *----------------------------------------------------------------------------*/
static int
spg_restart( ool_conmin_minimizer *M )
{
   def_st;
   size_t ii;
   
   /* Turns M->x into a feasible point */
   spg_proj( st->L, st->U, M->x );
   OOL_CONMIN_EVAL_FDF( M, M->x, &(M->f), M->gradient );

   /* alpha_0 */
   ii = gsl_blas_idamax( M->gradient );
   st->alpha = 1/fabs(gsl_vector_get( M->gradient, ii ));

   /* Infinite norm of d <- g1 = [P(x - g) - x] */
   gsl_vector_memcpy( st->d, M->x );
   gsl_blas_daxpy( -1, M->gradient, st->d );
   spg_proj( st->L, st->U, st->d );
   gsl_vector_sub( st->d, M->x );

   M->size = fabs(gsl_vector_get( st->d, gsl_blas_idamax( st->d )));

   st->m = 1;
   st->tail = 0;
   st->f[st->tail] = M->f;

   return OOL_SUCCESS;
}

/* Iterating
 *----------------------------------------------------------------------------*/
static int
spg_iterate( ool_conmin_minimizer *M )
{
   def_P;
   def_st;
   double bk, ak;
   
   spg_line_search( M );
   
   /* st->y = - (g(x_{k+1}) - g(x_k)) */
   gsl_blas_ddot( st->s, st->y, &bk );
   
   if (bk >= 0)
      {
	 st->alpha = P->alphamax;
      }
   else
      {
	 ak = gsl_blas_dnrm2( st->s );
	 ak = - ak * ak /bk;
	 st->alpha = GSL_MIN ( P->alphamax, GSL_MAX (P->alphamin, ak) );
      }

   return OOL_SUCCESS;
}

/* Optimality test
 *----------------------------------------------------------------------------*/
static int
spg_is_optimal( ool_conmin_minimizer *M )

{
   def_P;

   if( (M->size > P->tol) && (M->f > P->fmin) )
      return OOL_CONTINUE;
   else
      return OOL_SUCCESS;
}

/* Spectral line search
 *----------------------------------------------------------------------------*/
static void
spg_line_search( ool_conmin_minimizer *M )
{
   def_st;
   def_P;
  
   double lambda, lambdanew;
   double fmax, faux, fxx, dTg; 
   short armijo_failure;
   size_t ii;
  
   /* Saving the previous gradient */
   gsl_vector_memcpy( st->y, M->gradient );

   /* d = P(x - alpha * g) - x */
   gsl_vector_memcpy( st->d, M->x );
   gsl_blas_daxpy( -st->alpha, M->gradient, st->d );
   spg_proj( st->L, st->U, st->d );
   gsl_vector_sub( st->d, M->x );

   gsl_blas_ddot( st->d, M->gradient, &dTg );

   lambda = 1;
   armijo_failure = 1;
  
   while ( armijo_failure ){
      /* x trial */
      gsl_vector_memcpy( st->xx, M->x );
      gsl_blas_daxpy( lambda, st->d, st->xx );
      OOL_CONMIN_EVAL_F( M, st->xx, fxx );

      fmax = 0;
      for (ii = 0; ii < st->m; ii++)
	 {
	    faux = st->f[ii] + P->gamma * lambda * dTg;
	    fmax = GSL_MAX( fmax, faux );
	 }

      if ( fxx > fmax )
	 {
	    armijo_failure = 1;
	    /* build a quadratic model and minimize it */
	    lambdanew = - lambda * lambda * dTg / (2*(fxx - M->f -lambda*dTg));

	    lambda = GSL_MAX( P->sigma1 * lambda,
			      GSL_MIN ( P->sigma2 *lambda,
					lambdanew )
			     );
	 }
      else
	 {
	    armijo_failure = 0;
	 }
   }
 
   /* st->s = x_{k+1} - x_k */
   gsl_vector_memcpy( st->s, st->xx );
   gsl_vector_sub( st->s, M->x );

   /* M->x = x_{k+1} */
   gsl_vector_memcpy( M->x, st->xx );
   OOL_CONMIN_EVAL_DF( M, M->x, M->gradient );
   M->f = fxx;

   /* Infinite norm of g1 = d <- [P(x - g) - x] */
   gsl_vector_memcpy( st->d, M->x );
   gsl_blas_daxpy( -1, M->gradient, st->d );
   spg_proj( st->L, st->U, st->d );
   gsl_vector_sub( st->d, M->x );

   M->size = fabs(gsl_vector_get( st->d, gsl_blas_idamax( st->d )));

   /* st->y = - (g(x_{k+1}) - g(x_k)) */
   gsl_vector_sub( st->y, M->gradient );

   st->m++;
   st->m = GSL_MIN(st->M, st->m);

   st->tail++;
   st->tail = st->tail % st->M;
     st->f[ st->tail ] = M->f;

} 

/* Projection onto feasible region
 *----------------------------------------------------------------------------*/
void spg_proj( gsl_vector *L,
	       gsl_vector *U,
	       gsl_vector *X )
{

   double *l = L->data;
   double *u = U->data;
   double *x = X->data;

   size_t ii, n = X->size;
   
   for(ii=0; ii<n; ii++){
      *x = GSL_MAX (*l, GSL_MIN (*x, *u));
      x++;
      l++;
      u++;
   }
}

/* Setting spg type
 *----------------------------------------------------------------------------*/
static const ool_conmin_minimizer_type spg_type = 
   { "spg",
     sizeof( conmin_spg_state ),
     sizeof( ool_conmin_spg_parameters ),
     &spg_parameters_default,
     &spg_parameters_set,
     &spg_alloc,
     &spg_free,
     &spg_set,
     &spg_restart,
     &spg_iterate,
     &spg_is_optimal
   };

const ool_conmin_minimizer_type *ool_conmin_minimizer_spg = &spg_type;

/*----------------------------------------------------------------------------*/
