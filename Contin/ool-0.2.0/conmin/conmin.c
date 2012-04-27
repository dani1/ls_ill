/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * conmin/conmin.c
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
 * $Id: conmin.c,v 1.5 2005/05/07 19:52:30 biloti Exp $
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <gsl/gsl_errno.h>
#include <ool/ool_conmin_common.h>

/* Set default values for method parameters 
 *----------------------------------------------------------------------------*/
void
ool_conmin_parameters_default( const ool_conmin_minimizer_type *T, 
			       void  *parameters )
{
  (T->parameters_default)( parameters );
}

/* Set new method parameters
 *----------------------------------------------------------------------------*/
int 
ool_conmin_parameters_set( ool_conmin_minimizer *M, void *new_param )
{
  return (M->type->parameters_set)( M, new_param );
}

/* Get parameters in use
 *----------------------------------------------------------------------------*/
void
ool_conmin_parameters_get( const ool_conmin_minimizer *M,
			   void *param_view )
{
  memcpy( param_view, M->parameters, M->type->parameters_size );
}

/* Allocate new instance of a minimizer
 *----------------------------------------------------------------------------*/
ool_conmin_minimizer*
ool_conmin_minimizer_alloc( const ool_conmin_minimizer_type *T, 
			    const size_t n )
{
  int status;

  ool_conmin_minimizer *M;

  /* Allocation of minimizer */
  M = (ool_conmin_minimizer *) malloc( sizeof(ool_conmin_minimizer) );

  if( M == 0 )
    GSL_ERROR_VAL( "failed to allocate space for minimizer", OOL_ENOMEM, 0 );

  /* Store pointer to minimizer type */
  M->type = T;

  /* Allocation of x */
  M->x = gsl_vector_calloc( n );

  if( M->x == 0 ) 
    {
      free( M );
      GSL_ERROR_VAL ("failed to allocate space for x", OOL_ENOMEM, 0);
    }

  /* Allocation of gradient */
  M->gradient = gsl_vector_calloc( n );

  if( M->gradient == 0 ) 
    {
      gsl_vector_free( M->x );
      free( M );
      GSL_ERROR_VAL ("failed to allocate space for gradient", OOL_ENOMEM, 0);
    }

  /* Allocation of dx */
  M->dx = gsl_vector_calloc( n );

  if( M->dx == 0 ) 
    {
      gsl_vector_free ( M->x );
      gsl_vector_free ( M->gradient );
      free( M );
      GSL_ERROR_VAL ("failed to allocate space for dx", OOL_ENOMEM, 0);
    }

  /* Allocation of minimizer parameters */
  M->parameters = malloc( T->parameters_size );

  if( M->parameters == 0 )
    {
      gsl_vector_free( M->x );
      gsl_vector_free( M->gradient );
      gsl_vector_free( M->dx );
      free( M );
      GSL_ERROR_VAL( "failed to allocate space for minimizer state", 
		     OOL_ENOMEM, 0);
    }

  /* Allocation of minimizer state */
  M->state = malloc( T->state_size );

  if( M->state == 0 )
    {
      gsl_vector_free( M->x );
      gsl_vector_free( M->gradient );
      gsl_vector_free( M->dx );
      free( M->parameters );
      free( M );
      GSL_ERROR_VAL( "failed to allocate space for minimizer state", 
		     OOL_ENOMEM, 0);
    }

  /* Allocation of internal state variables */
  status = (T->alloc)( M->state, n );

  if( status != OOL_SUCCESS )
    {
      gsl_vector_free( M->x );
      gsl_vector_free( M->gradient );
      gsl_vector_free( M->dx );
      free( M->parameters );
      free( M->state );
      free( M );
      GSL_ERROR_VAL ("failed to initialize minimizer state", OOL_ENOMEM, 0);
    }

  return M;
}

/* Frees al the memory associated to minimizer
 *----------------------------------------------------------------------------*/
void ool_conmin_minimizer_free( ool_conmin_minimizer *M )
{
  (M->type->free)( M->state );

  free( M->state );
  free( M->parameters );

  gsl_vector_free( M->dx       );
  gsl_vector_free( M->gradient );
  gsl_vector_free( M->x        );

  free( M );
}

/* Initializes the  minimizer
 *----------------------------------------------------------------------------*/
int 
ool_conmin_minimizer_set( ool_conmin_minimizer  *M,
			  ool_conmin_function   *fdf,
			  ool_conmin_constraint *con,
			  gsl_vector            *x,
			  void                  *parameters )
{
  int status;

  if( M->x->size != fdf->n )
    GSL_ERROR( "function incompatible with solver size", OOL_EBADLEN );
  
  if( x->size != fdf->n ) 
    GSL_ERROR( "vector length not compatible with function", OOL_EBADLEN );

  if( x->size != con->n ) 
    GSL_ERROR( "vector length not compatible with constraint", OOL_EBADLEN );

  M->fdf = fdf;
  M->con = con;

  status = (M->type->parameters_set)( M, parameters );

  if( status != OOL_SUCCESS )
    {
      GSL_ERROR( "Improper parameters", OOL_EINVAL );
    }

  gsl_vector_memcpy  ( M->x, x );
  gsl_vector_set_zero( M->dx   );
  
  /* Start evaluation conunters */
  M->fcount = 0;
  M->gcount = 0;
  M->hcount = 0;

  return (M->type->set)( M );
}

/* Restart Method Instance
 *----------------------------------------------------------------------------*/
int ool_conmin_minimizer_restart( ool_conmin_minimizer *M )
{
  /* Reset dx vector */
  gsl_vector_set_zero( M->dx );

  /* Restart evaluation counters */
  M->fcount = 0;
  M->gcount = 0;
  M->hcount = 0;

  return (M->type->restart)( M );
}

/* Iteration
 *----------------------------------------------------------------------------*/
int ool_conmin_minimizer_iterate( ool_conmin_minimizer *M )
{
  return (M->type->iterate)( M );
}

/* Optimality check of last visited point 
 *----------------------------------------------------------------------------*/
int 
ool_conmin_is_optimal( ool_conmin_minimizer *M )
{
  return (M->type->is_optimal)( M );
}

/* State information
 *----------------------------------------------------------------------------*/
gsl_vector* 
ool_conmin_minimizer_x( const ool_conmin_minimizer *M )
{
  return M->x;
}

gsl_vector* 
ool_conmin_minimizer_dx( const ool_conmin_minimizer *M )
{
  return M->dx;
}

double 
ool_conmin_minimizer_minimum( const ool_conmin_minimizer *M )
{
  return M->f;
}

gsl_vector* 
ool_conmin_minimizer_gradient( const ool_conmin_minimizer *M )
{
  return M->gradient;
}

const char* 
ool_conmin_minimizer_name( const ool_conmin_minimizer *M )
{
  return M->type->name;
}

double 
ool_conmin_minimizer_size( const ool_conmin_minimizer *M )
{
  return M->size;
}

size_t 
ool_conmin_minimizer_fcount ( const ool_conmin_minimizer *M )
{
  return M->fcount;
}

size_t 
ool_conmin_minimizer_gcount ( const ool_conmin_minimizer *M )
{
  return M->gcount;
}

size_t 
ool_conmin_minimizer_hcount ( const ool_conmin_minimizer *M )
{
  return M->hcount;
}

/*----------------------------------------------------------------------------*/
