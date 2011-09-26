/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * tools/ools_neval.c
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
 * since: May, 23, 2004
 *
 * $Id: tools_diff.c,v 1.2 2004/06/04 21:26:23 akiles Exp $
 *----------------------------------------------------------------------------*/

#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <ool/ool_tools_diff.h>

typedef struct {

   double(*F)( const gsl_vector*, void* );
   gsl_vector *X;
   void *PF;

   size_t i;

} tools_f_t;

/* Evaluate F(X), but with the prototype of gsl_function
 *----------------------------------------------------------------------------*/
double tools_f( double    x,
		void      *Pv )
{
   double xi;
   double fval;
   
   tools_f_t *P = (tools_f_t *) Pv;
   
   xi = gsl_vector_get( P->X, P->i );

   gsl_vector_set( P->X, P->i, x );
   fval = P->F( P->X, P->PF );

   gsl_vector_set( P->X, P->i, xi );
   
   return fval;

}

/* Evaluate numerically the gradient of f
   
   A good choice for eps is

   eps = (eta)^(1/3) * ||x||_inf

   where eta is the precision on the evaluation of f
 *----------------------------------------------------------------------------*/
int ool_diff_g( double(*F)( const gsl_vector*, void* ),
		const gsl_vector *X,
		void             *fparam,
		gsl_vector       *G,
		double            eps     )
{
   size_t ii;
   double xi, dfi;

   for (ii=0; ii < X->size; ii++){
      xi = gsl_vector_get( X, ii );

      gsl_vector_set( X, ii, xi + eps );
      dfi = F( X, fparam);

      gsl_vector_set( X, ii, xi - eps );
      dfi -= F( X, fparam );

      gsl_vector_set( G, ii, dfi/(2*eps) );
      gsl_vector_set( X, ii, xi );
   }

  return GSL_SUCCESS;
}

/* Evaluate numerically the gradient of f
   with auto selection of the step
   through the gsl_diff_central function
 *----------------------------------------------------------------------------*/
int ool_diff_g_auto( double(*F)( const gsl_vector*, void* ),
		     const gsl_vector *X,
		     void             *Fparam,
		     gsl_vector       *G       )
{
   size_t ii;
   double dfi, err, xi; 
   tools_f_t fparam;
   gsl_function func;

   func.function = &tools_f;
   func.params = &fparam;

   fparam.F = F;
   fparam.PF = Fparam;
   fparam.X = X;
   
   for (ii=0; ii < X->size; ii++){
      fparam.i = ii;
      xi = gsl_vector_get( X, ii );

      gsl_diff_central(&func, xi, &dfi, &err );

      gsl_vector_set( G, ii, dfi );
   }

   return GSL_SUCCESS;
}

/* Acceleration alloc for evaluation of hessian
 *----------------------------------------------------------------------------*/
ool_diff_Hv_accel * ool_diff_Hv_accel_alloc( size_t n )
{
   ool_diff_Hv_accel *a;
   
   a = (ool_diff_Hv_accel *) malloc(n * sizeof(ool_diff_Hv_accel));

   a->gradf1 = gsl_vector_alloc( n );
   a->gradf2 = gsl_vector_alloc( n );

   return a;
}

/* Acceleration free for evalutation of hessian
 *----------------------------------------------------------------------------*/
void ool_diff_Hv_accel_free( ool_diff_Hv_accel *a )
{
   gsl_vector_free(a->gradf1);
   gsl_vector_free(a->gradf2);

   free(a);
}

/* Evaluate numerically the hessian of f times a vector
 *----------------------------------------------------------------------------*/
int ool_diff_Hv( const ool_diff_Hv_accel *a,
		 void(*df)( const gsl_vector*, void*, gsl_vector* ),
		 const gsl_vector        *X,
		 void                    *fparam,
		 const gsl_vector        *V,
		 gsl_vector              *Hv,
		 double                   eps     )
{
   double xi, dfij;
   size_t i;
   gsl_vector *gradf1, *gradf2;

   if (a == NULL){
      gradf1 = gsl_vector_alloc( X->size );
      gradf2 = gsl_vector_alloc( X->size );
   }
   else{
      gradf1 = a->gradf1;
      gradf2 = a->gradf2;
   }

   for (i=0; i < X->size; i++){

      xi = gsl_vector_get( X, i );

      gsl_vector_set ( X, i, xi + eps );
      df( X, fparam, gradf1 );

      gsl_vector_set ( X, i, xi - eps );
      df( X, fparam, gradf2 );

      gsl_vector_set( X, i, xi );

      gsl_vector_sub( gradf1, gradf2 );
      gsl_blas_ddot ( gradf1, V, &dfij );
      
      gsl_vector_set( Hv, i, dfij/(2*eps) );

   }
   
   if ( a == NULL ){
      gsl_vector_free( gradf1 );
      gsl_vector_free( gradf2 );
   }

   return GSL_SUCCESS;
}

/*----------------------------------------------------------------------------*/
