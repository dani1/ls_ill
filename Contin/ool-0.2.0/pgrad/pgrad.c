/*----------------------------------------------------------------------------*
 * Open Optimization Library - Projected Gradient Method
 * 
 * pgrad/pgrad.c
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
 * Iara da Cunha
 * since: June, 29, 2004 
 *
 * $Id: pgrad.c,v 1.14 2005/05/19 19:37:07 biloti Exp $
 *----------------------------------------------------------------------------*/
#include <pgrad.h>

#define def_P ool_conmin_pgrad_parameters*  P = \
             (ool_conmin_pgrad_parameters*) M->parameters

#define def_st conmin_pgrad_state*  st = \
              (conmin_pgrad_state*) M->state

/*----------------------------------------------------------------------------*
  Gradient Projection
  *----------------------------------------------------------------------------*/

/* Set default values for pgrad parameters
 *----------------------------------------------------------------------------*/
static void
pgrad_parameters_default( void  *vpar )
{
  ool_conmin_pgrad_parameters *P = (ool_conmin_pgrad_parameters*) vpar;

  /* Setting default values of parameters */
  P->fmin = -1.0e+99;
  P->tol  =  1.0e-4;
  P->alpha=  1.0e-4; /* Armijo */

  P->sigma1 =  0.1;
  P->sigma2 =  0.9;
}

/* Set values for pgrad parameters
 *----------------------------------------------------------------------------*/
static int
pgrad_parameters_set( ool_conmin_minimizer *M, void *new_param )
{
  /* Copy new state values */
  memcpy( M->parameters, new_param, sizeof(ool_conmin_pgrad_parameters) );

  return OOL_SUCCESS;
}

/* Allocation  
 *----------------------------------------------------------------------------*/
static int
pgrad_alloc( void *vstate, size_t n )
{
  conmin_pgrad_state *state = (conmin_pgrad_state*)vstate;

  /* Space dimension */
  state->n = n;

  state->xx = gsl_vector_alloc( n );
  state->L  = gsl_vector_alloc( n );
  state->U  = gsl_vector_alloc( n );
  
  if( state->xx == NULL || 
      state->L  == NULL || 
      state->U  == NULL  )
    {
      gsl_vector_free( state->xx );
      gsl_vector_free( state->L  );
      gsl_vector_free( state->U  );

      return OOL_ENOMEM;
    }
  
  return OOL_SUCCESS;
}

/* Free memory allocated by pgrad_alloc
 *----------------------------------------------------------------------------*/
static void
pgrad_free( void *vstate )
{
  conmin_pgrad_state *state = (conmin_pgrad_state*)vstate;

  gsl_vector_free( state->xx );
  gsl_vector_free( state->L  );
  gsl_vector_free( state->U  );  
}

/* Setting
 *----------------------------------------------------------------------------*/
static int
pgrad_set( ool_conmin_minimizer *M )
{
  def_st;

  /* Copy boundary vectors */
  gsl_vector_memcpy( st->L, M->con->L );
  gsl_vector_memcpy( st->U, M->con->U );

  /* Turns M->x into a feasible point */
  pgrad_proj( st->L, st->U, M->x );
  
  OOL_CONMIN_EVAL_FDF( M, M->x, &(M->f), M->gradient );

  return OOL_SUCCESS;
}

/* Restarting
 *----------------------------------------------------------------------------*/
static int
pgrad_restart( ool_conmin_minimizer *M )
{
  def_st;
   
  /* Turns M->x into a feasible point */
  pgrad_proj( st->L, st->U, M->x );

  OOL_CONMIN_EVAL_FDF( M, M->x, &(M->f), M->gradient );

  return OOL_SUCCESS;
}

/* Iterating
 *----------------------------------------------------------------------------*/
static int
pgrad_iterate( ool_conmin_minimizer *M )
{
  def_st;
  double t;
  
  pgrad_line_search( M );
  
  /* Infinite norm of g1 = d <- [P(x - g) - x] */
  gsl_vector_memcpy( st->xx, M->x );
  gsl_vector_sub( st->xx, M->gradient );
  pgrad_proj( st->L, st->U, st->xx );
  gsl_vector_sub( st->xx, M->x );
  
  M->size = fabs(gsl_vector_get( st->xx, gsl_blas_idamax( st->xx )));

  return OOL_SUCCESS;
 }

/* Optimality test
 *----------------------------------------------------------------------------*/
static int
pgrad_is_optimal( ool_conmin_minimizer *M )

{
  def_st;
  def_P;

  if( M->size > P->tol && M->f > P->fmin )
    return OOL_CONTINUE;
  else
    return OOL_SUCCESS;
}

/* Line search (as in Kelley's book)
 *----------------------------------------------------------------------------*/
static void
pgrad_line_search( ool_conmin_minimizer *M )
{
  def_st;
  def_P;
  
  double t, tnew, fx, fxx, dif2, gtd;
   
  fx = M->f;

  tnew = 1;

  do{

     t = tnew;
     /* xx = x - t df */
     gsl_vector_memcpy( st->xx, M->x );
     gsl_blas_daxpy( -t, M->gradient, st->xx );  
     pgrad_proj( st->L, st->U, st->xx );
     OOL_CONMIN_EVAL_F( M, st->xx, fxx );

     gsl_vector_memcpy( M->dx, st->xx );
     gsl_vector_sub( M->dx, M->x );

     dif2 = gsl_blas_dnrm2( M->dx );
     dif2 = dif2 * dif2;

     gsl_blas_ddot( M->gradient, M->dx, &gtd );
     
     tnew = - t * t * gtd / (2*(fxx - M->f - gtd));

     double arg1 = P->sigma1 * t;
     double arg2 = GSL_MIN ( P->sigma2 * t, tnew );
     tnew = GSL_MAX( arg1, arg2 );

  /* sufficient decrease condition (Armijo) */
  }while( fxx > fx - (P->alpha/t) * dif2 );

  gsl_vector_memcpy( M->x, st->xx );
  M->f = fxx;
  OOL_CONMIN_EVAL_DF( M, M->x, M->gradient );
}

/* Projection onto feasible region
   a = max(b, min(c, d))
 *----------------------------------------------------------------------------*/
void pgrad_proj( gsl_vector *L,
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

/* Setting pgrad type
 *----------------------------------------------------------------------------*/
static const ool_conmin_minimizer_type pgrad_type = 
  { "pgrad",
    sizeof( conmin_pgrad_state ),
    sizeof( ool_conmin_pgrad_parameters ),
    &pgrad_parameters_default,
    &pgrad_parameters_set,
    &pgrad_alloc,
    &pgrad_free,
    &pgrad_set,
    &pgrad_restart,
    &pgrad_iterate,
    &pgrad_is_optimal
  };

const ool_conmin_minimizer_type *ool_conmin_minimizer_pgrad = &pgrad_type;

/*----------------------------------------------------------------------------*/
