/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * conmin/ool_conmin_common.h
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
 * since: Feb, 13, 2004
 *
 * $Id: ool_conmin_common.h,v 1.8 2005/04/20 19:38:49 biloti Exp $
 *----------------------------------------------------------------------------*/

#ifndef __OOL_CONMIN_COMMON_H__
#define __OOL_CONMIN_COMMON_H__

#include <stdlib.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_vector.h>

#include <ool/ool_conmin_defs.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS   /* empty */
#endif

__BEGIN_DECLS

/* Target function
 *----------------------------------------------------------------------------*/
typedef struct
{
  /* Dimention of the system */
  size_t n;

  /* Return the function value at point X */
  double(*f)( const gsl_vector *X, 
	      void  *params );

  /* Evaluate the gradient of the function at X */
  void( *df)( const gsl_vector *X, 
	      void  *params, 
	      gsl_vector       *gradient );

  /* Evaluate the value and gradient of the function at X */
  void(*fdf)( const gsl_vector *X, 
	      void  *params, 
	      double           *f, 
	      gsl_vector       *gradient );

  /* Evaluate the Hessian times vector V product */
  void( *Hv)( const gsl_vector *X, 
	      void  *params, 
	      const gsl_vector *V,
	      gsl_vector       *hv );

  /* Parameters of the function */
  void *params;

} ool_conmin_function;

/* Evaluation of the target function */
#define OOL_CONMIN_EVAL_F( M, x, F )       \
  F = (*(M->fdf->f))( x, M->fdf->params ); \
  (M->fcount)++

#define OOL_CONMIN_EVAL_DF( M, x, g )      \
  (*(M->fdf->df))( x, M->fdf->params, g ); \
  (M->gcount)++

#define OOL_CONMIN_EVAL_FDF( M, x, f, g )      \
  (*(M->fdf->fdf))( x, M->fdf->params, f, g ); \
  (M->fcount)++;                                \
  (M->gcount)++

#define OOL_CONMIN_EVAL_HV( M, x, v, Hxv )      \
  (*(M->fdf->Hv))( x, M->fdf->params, v, Hxv ); \
  (M->hcount)++

/* Constraints 
 *----------------------------------------------------------------------------*/
typedef struct
{
  /* Dimension of the system */
  size_t n;

  /* Box constraints */
  gsl_vector *L;
  gsl_vector *U;

  /* Number of inequality and equality constraints */
  size_t ng;
  size_t nh;

  /* Inequality constraints */
  double(*g)( const size_t      i,
	      const gsl_vector *X,
	      const void       *g_params );

  /* Equality constraints */
  double(*h)( const size_t      i, 
	      const gsl_vector *X,
	      const void       *h_params );

  /* Tolerance on equality constraints */
  double eq_tol;

  /* Parameters */
  void *g_params;
  void *h_params;

} ool_conmin_constraint;

struct ool_conmin_minimizer_type_struct;

/* Constrained Minimizer
 *----------------------------------------------------------------------------*/
typedef struct 
{
  const struct ool_conmin_minimizer_type_struct *type;

  ool_conmin_function   *fdf;
  ool_conmin_constraint *con;

  double f;
  double size;

  gsl_vector *x;
  gsl_vector *gradient;
  gsl_vector *dx;

  void *parameters;
  void *state;

  size_t fcount;
  size_t gcount;
  size_t hcount;

} ool_conmin_minimizer;

/* Constrained Minimizer Type
 *----------------------------------------------------------------------------*/
typedef struct ool_conmin_minimizer_type_struct
{
  /* Method name */
  const char *name;

  /* Struct sizes */
  size_t state_size;
  size_t parameters_size;

  /* Parameters manipulation */
  void (*parameters_default)( void *parameters );
  int  (*parameters_set    )( ool_conmin_minimizer *M, void *newparam );

  /* Manipulating memory for state struct */
  int  (*alloc)( void *state, size_t n );
  void (*free )( void *state );

  /* Method functions */
  int (*set       )( ool_conmin_minimizer *M );
  int (*restart   )( ool_conmin_minimizer *M );
  int (*iterate   )( ool_conmin_minimizer *M );
  int (*is_optimal)( ool_conmin_minimizer *M );

} ool_conmin_minimizer_type;

/* Parameters manipulation
 *----------------------------------------------------------------------------*/

/* Set default values for method parameters */
void
ool_conmin_parameters_default( const ool_conmin_minimizer_type *T, 
			       void  *parameters );

/* Set new parameters */
int 
ool_conmin_parameters_set( ool_conmin_minimizer *M, 
			   void  *new_param );

/* Get parameters in use */
void
ool_conmin_parameters_get( const ool_conmin_minimizer *M, 
			   void  *param_view );

/* Allocation, Initialization and Finalization 
 *----------------------------------------------------------------------------*/

/* Allocate new instance of a minimizer */
ool_conmin_minimizer*
ool_conmin_minimizer_alloc( const ool_conmin_minimizer_type *T, 
			    const size_t n );

/* Initializes the  minimizer */
int 
ool_conmin_minimizer_set( ool_conmin_minimizer  *M,
			  ool_conmin_function   *fdf,
			  ool_conmin_constraint *con,
			  gsl_vector            *x,
			  void                  *parameters  );

/* Frees al the memory associated to minimizer */
void ool_conmin_minimizer_free( ool_conmin_minimizer *M );

/* Iterations
 *----------------------------------------------------------------------------*/

/* Iterates Method */
int ool_conmin_minimizer_iterate( ool_conmin_minimizer *M );

/* Restart Method Instance */
int ool_conmin_minimizer_restart( ool_conmin_minimizer *M );

/* State information
 *----------------------------------------------------------------------------*/

/* Return the last point */
gsl_vector* ool_conmin_minimizer_x( const ool_conmin_minimizer *M );

/* Return the last step */
gsl_vector* ool_conmin_minimizer_dx( const ool_conmin_minimizer *M );

/* Return the function value */
double ool_conmin_minimizer_minimum( const ool_conmin_minimizer *M );

/* Return function gradient */
gsl_vector* ool_conmin_minimizer_gradient( const ool_conmin_minimizer *M );

/* Return the method name */
const char* ool_conmin_minimizer_name( const ool_conmin_minimizer *M );

/* Return a method dependent size value */
double ool_conmin_minimizer_size( const ool_conmin_minimizer *M );

/* Return the number of evaluations of objective function */
size_t ool_conmin_minimizer_fcount( const ool_conmin_minimizer *M );

/* Return the number of evaluations of gradient */
size_t ool_conmin_minimizer_gcount( const ool_conmin_minimizer *M );

/* Return the number of evaluations of hessian */
size_t ool_conmin_minimizer_hcount ( const ool_conmin_minimizer *M );

/* Optimality tests
 *----------------------------------------------------------------------------*/

/* Optimality check of last visited point */
int ool_conmin_is_optimal( ool_conmin_minimizer *M );

/*----------------------------------------------------------------------------*/
__END_DECLS

#endif /*__OOL_CONMIN_COMMON_H__*/
/*----------------------------------------------------------------------------*/
