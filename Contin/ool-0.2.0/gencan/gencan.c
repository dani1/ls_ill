/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * gencan/gencan.c
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
 * since: Feb 16th, 2004
 *
 * $Id: gencan.c,v 1.16 2005/05/17 19:08:18 biloti Exp $
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <gencan.h>

#define def_P ool_conmin_gencan_parameters*  P = \
             (ool_conmin_gencan_parameters*) M->parameters

#define def_st conmin_gencan_state*  st = \
              (conmin_gencan_state*) M->state

/* Set default values for gencan parameters 
 *----------------------------------------------------------------------------*/
static void gencan_parameters_default( void  *vpar )
{
  ool_conmin_gencan_parameters *P = (ool_conmin_gencan_parameters*) vpar;

  /* Setting default values of parameters */
  P->epsgpen   =  1.0e-05;
  P->epsgpsn   =  1.0e-05;

  P->fmin      = -1.0e+99;

  P->udelta0   = -1;
  P->ucgmia    = -1;
  P->ucgmib    = -1;
  P->cg_scre   =  1;
  P->cg_gpnf   =  P->epsgpen;
  P->cg_epsi   =  1.0e-1;
  P->cg_epsf   =  1.0e-5;

  P->cg_epsnqmp= 1.0e-4;
  P->cg_maxitnqmp = 5;

  P->nearlyq   =  0;

  P->nint      = 2.0;
  P->next      = 2.0;
  P->mininterp = 4;
  P->maxextrap = 100;

  P->trtype    = 0;

  P->eta       = 0.9;
  P->delmin    = 0.1;

  P->lspgmi    = 1.0e-10;
  P->lspgma    = 1.0e+10;

  P->theta     = 1.0e-06;
  P->gamma     = 1.0e-04; 
  P->beta      = 0.5;

  P->sigma1    = 0.1;
  P->sigma2    = 0.9;

  P->epsrel    = 1.0e-07;
  P->epsabs    = 1.0e-10;

  P->infrel    = 1.0e+20;
  P->infabs    = 1.0e+99;

}

/* Set values for gencan parameters 
 *----------------------------------------------------------------------------*/
static int gencan_parameters_set( ool_conmin_minimizer *M, void *new_param )
{
   ool_conmin_gencan_parameters*  P = 
      (ool_conmin_gencan_parameters*) new_param;

   /* Check parameters */

   if ( P->epsgpsn < 0              ||
	P->epsgpen < 0              ||
	P->cg_gpnf < 0              ||
	P->cg_epsi < 0              ||
	P->cg_epsf < 0              ||
	P->cg_epsnqmp < 0           ||
	P->cg_maxitnqmp < 1         ||
	( P->nearlyq != 0 &&
	  P->nearlyq != 1   )       ||
	P->nint <= 1                ||
	P->next <= 1                ||
	P->mininterp < 1            ||
	( P->trtype != 0 &&
	  P->trtype != 1   )        ||
	( P->eta <= 0 ||
	  P->eta >= 1  )            ||
	P->delmin <= 0              ||
	P->lspgmi <= 0              ||
	P->lspgma < P->lspgmi       ||
	( P->theta <= 0 ||
	  P->theta >= 1  )          ||
	( P->gamma <= 0 ||
	  P->gamma >= 0.5 )         ||
	( P->beta <= 0 ||
	  P->beta >= 1  )           ||
	P->sigma1 <= 0              ||
	( P->sigma2 <= P->sigma1 ||
	  P->sigma2 >=1           ) ||
	P->epsrel < 0               ||
	P->epsabs < 0               ||
	P->infrel < 0               ||
	P->infabs < 0
	)
      return OOL_EINVAL;

   /* Copy new parameter values */
   memcpy( M->parameters, new_param, sizeof(ool_conmin_gencan_parameters) );
   
   return OOL_SUCCESS;
}

/* Allocation  
 *----------------------------------------------------------------------------*/
static int gencan_alloc( void *vstate, size_t n )
{
  conmin_gencan_state *st = (conmin_gencan_state *) vstate;

  /* Space dimension */
  st->n = n;

  /* Allocation of working vectors. */
  /* We should test of proper allocation after each vector allocation! */

  st->Ind = gsl_vector_uint_alloc( n );

  st->L = gsl_vector_calloc( n );
  st->U = gsl_vector_calloc( n );
  st->S = gsl_vector_alloc( n );
  st->Y = gsl_vector_alloc( n );
  st->D = gsl_vector_alloc( n );

  st->near_l = (double*) malloc( sizeof(double) * n );
  st->near_u = (double*) malloc( sizeof(double) * n );

  st->Xtrial = gsl_vector_alloc( n );

  st->tnls_Xtemp = gsl_vector_alloc( n );

  st->cg_W     = gsl_vector_alloc( n );
  st->cg_R     = gsl_vector_alloc( n );
  st->cg_D     = gsl_vector_alloc( n );
  st->cg_Sprev = gsl_vector_alloc( n );

  if( st->Ind          == NULL  ||
      st->L            == NULL  ||
      st->U            == NULL  ||
      st->S            == NULL  ||
      st->Y            == NULL  ||
      st->D            == NULL  ||
      st->near_l       == NULL  ||
      st->near_u       == NULL  ||
      st->Xtrial       == NULL  ||
      st->tnls_Xtemp   == NULL  ||
      st->cg_W         == NULL  ||
      st->cg_R         == NULL  ||
      st->cg_D         == NULL  ||
      st->cg_Sprev     == NULL    )
    {
      gsl_vector_uint_free( st->Ind );

      gsl_vector_free( st->L );
      gsl_vector_free( st->U );
      gsl_vector_free( st->S );
      gsl_vector_free( st->Y );
      gsl_vector_free( st->D );

      free(st->near_l);
      free(st->near_u);

      gsl_vector_free( st->Xtrial );

      gsl_vector_free( st->tnls_Xtemp );

      gsl_vector_free( st->cg_W );
      gsl_vector_free( st->cg_R );
      gsl_vector_free( st->cg_D );
      gsl_vector_free( st->cg_Sprev );

      GSL_ERROR_VAL( "failed to allocate space for gencan working vectors",
		     OOL_ENOMEM, 0 );
    }

  return OOL_SUCCESS;
}

/* Free memory allocated by gencan_alloc
 *----------------------------------------------------------------------------*/
static void gencan_free( void *vstate )
{
  conmin_gencan_state *st = (conmin_gencan_state *) vstate;

  gsl_vector_uint_free( st->Ind );

  gsl_vector_free( st->L );
  gsl_vector_free( st->U );

  gsl_vector_free( st->S );
  gsl_vector_free( st->Y );
  gsl_vector_free( st->D );

  free(st->near_l);
  free(st->near_u);
      
  gsl_vector_free( st->Xtrial );

  gsl_vector_free( st->tnls_Xtemp );

  gsl_vector_free( st->cg_W );
  gsl_vector_free( st->cg_R );
  gsl_vector_free( st->cg_D );
  gsl_vector_free( st->cg_Sprev );
}

/* Define state variables and prepare to initiate iterations
 *----------------------------------------------------------------------------*/
int gencan_prepare_iteration( ool_conmin_minimizer *M )
{
  def_st;
  def_P;

  /* Direct access to vector data */
  double *x =  M->x->data;
  double *l = st->L->data;
  double *u = st->U->data;

  /* Internal variables */
  size_t nn = M->x->size;
  size_t ii, imax;

  /* Impose factibility */
  conmin_vector_maxofmin( nn, x, l, u, x );

  /* Eval Euclidean and infinity norms of X */
  st->xeucn = gsl_blas_dnrm2 ( M->x );
  imax      = gsl_blas_idamax( M->x );
  st->xsupn = fabs( gsl_vector_get( M->x, imax ) );

  /* Evaluate objective function and its gradient */
  OOL_CONMIN_EVAL_FDF( M, M->x, &(M->f), M->gradient );

  /* Define near_l and near_u vector */
  for (ii=0; ii < nn; ii++){
     st->near_l[ii] = l[ii] + GSL_MAX( P->epsrel*fabs( l[ii] ), P->epsabs );
     st->near_u[ii] = u[ii] - GSL_MAX( P->epsrel*fabs( u[ii] ), P->epsabs );
  }
  
  /* Setting constant parameters */
  st->ometa2 = gsl_pow_2( 1.0 - P->eta );
  st->epsgpen2 = gsl_pow_2( P->epsgpen );

  /* Compute continuous project gradient */
  gencan_projected_gradient( st, M->x, M->gradient );

  /* Compute a linear relation between gpeucn2 and cgeps2, i.e.,
   * scalars a and b such that 
   *
   *     a * log10(||g_P(x_0)||_2^2) + b = log10(cgeps_0^2) and
   *
   *     a * log10(||g_P(x_f)||_2^2) + b = log10(cgeps_f^2),
   *
   *  where cgeps_0 and cgeps_f are provided. Note that if 
   *  cgeps_0 is equal to cgeps_f then cgeps will be always 
   *  equal to cgeps_0 and cgeps_f.
   *  
   *  We introduce now a linear relation between gpsupn and cgeps also.
   */
  if ( P->cg_scre == 1 ){
     st->acgeps = 2 *( log10( P->cg_epsf / P->cg_epsi ) /
		       log10( P->cg_gpnf * P->cg_gpnf / st->gpeucn2 ));
     
     st->bcgeps = 2 * log10( P->cg_epsi ) -
	          st->acgeps * log10( st->gpeucn2 );
  }
  else {
     st->acgeps = ( log10( P->cg_epsf / P->cg_epsi ) /
		    log10( P->cg_gpnf / st->gpsupn ) );
     st->bcgeps = log10( P->cg_epsi ) - st->acgeps * log10( st->gpsupn );
  }

  /*     And it will be used for the linear relation of cgmaxit */
  st->gpsupn0  = st->gpsupn;
  st->gpeucn20 = st->gpeucn2;

  /* Initialize the spectral steplength */
  if ( st->gpeucn2 != 0.0 ) 
    st->lambda = GSL_MAX( 1.0, st->xeucn ) / sqrt( st->gpeucn2 );

  /* Initialize the trust-region radius */
  if( P->udelta0 < 0.0 ) 
    {
       double aux;
       if ( P->trtype ) aux = 0.1 * GSL_MAX( 1.0, st->xeucn );
       else             aux = 0.1 * GSL_MAX( 1.0, st->xsupn );

       st->cg_delta = GSL_MAX( P->delmin, aux );

    }
  else st->cg_delta = GSL_MAX( P->delmin, P->udelta0 );

  /* Export size */
  M->size = st->gpsupn;

  return OOL_SUCCESS;
}

/* Setting
 *----------------------------------------------------------------------------*/
static int gencan_set( ool_conmin_minimizer *M )
{
  def_st;

  /* Internal variables */
  size_t nn = M->x->size;

  /* Checking dimensions */
  if( nn != st->n || nn != M->fdf->n || nn != M->con->n  )
    {
      return OOL_EINVAL;
    }

  /* Copy boundary vectors */
  gsl_vector_memcpy( st->L, M->con->L );
  gsl_vector_memcpy( st->U, M->con->U );

  return gencan_prepare_iteration( M );
}

/* Restarting
 *----------------------------------------------------------------------------*/
static int gencan_restart( ool_conmin_minimizer *M )
{
  /* Restarting dx */
  gsl_vector_set_zero( M->dx );

  return gencan_prepare_iteration( M );
}

/* Iterating
 * Actual code in file gencan_iterate.c
 *----------------------------------------------------------------------------*/
static int gencan_iterate( ool_conmin_minimizer *M )
{
  def_st;
  def_P;

  int status;

  status = gencan_actual_iterate( M, st, P );

  /* Export size and dx variables */
  M->size = st->gpsupn;

  /* In the next version does dx replace st->S ? */
  gsl_vector_memcpy( M->dx, st->S );

  return status;
}

/* Optimality test
 *----------------------------------------------------------------------------*/
static int gencan_is_optimal( ool_conmin_minimizer *M )
{
  def_st;
  def_P;

  return (( st->gpeucn2 <= st->epsgpen2 || 
	    st->gpsupn  <= P->epsgpsn   ||
	    M->f        <= P->fmin      )? OOL_SUCCESS : OOL_CONTINUE );
}

/* Setting gencan type
 *----------------------------------------------------------------------------*/
static const ool_conmin_minimizer_type gencan_type = 
  { "gencan",
    sizeof( conmin_gencan_state ),
    sizeof( ool_conmin_gencan_parameters ),
    &gencan_parameters_default,
    &gencan_parameters_set,
    &gencan_alloc,
    &gencan_free,
    &gencan_set,
    &gencan_restart,
    &gencan_iterate,
    &gencan_is_optimal,
};

const ool_conmin_minimizer_type *ool_conmin_minimizer_gencan = &gencan_type;

/*----------------------------------------------------------------------------*/
