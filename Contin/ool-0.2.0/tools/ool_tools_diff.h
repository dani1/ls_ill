/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * tools/ool_tools_neval.h
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
 * $Id: ool_tools_diff.h,v 1.2 2004/06/04 21:26:23 akiles Exp $
 *----------------------------------------------------------------------------*/

#ifndef __OOL_DIFF_H__
#define __OOL_DIFF_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_diff.h>

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


/* Hessian evaluation acceleration struct *
 *----------------------------------------------------------------------------*/
typedef struct{

   gsl_vector *gradf1;
   gsl_vector *gradf2;

} ool_diff_Hv_accel;

/* Eval numerically the gradient of f
 *----------------------------------------------------------------------------*/
int ool_diff_g( double(*f)( const gsl_vector*, void* ),
		const gsl_vector *X,
		void             *fparam,
		gsl_vector       *G,
		const double      eps     );


/* Eval numerically the gradient of f
 *----------------------------------------------------------------------------*/
int ool_diff_g_auto( double(*f)( const gsl_vector*, void* ),
		     const gsl_vector *X,
		     void             *fparam,
		     gsl_vector       *G       );

/* Acceleration alloc for evaluation of hessian
 *----------------------------------------------------------------------------*/
ool_diff_Hv_accel * ool_diff_Hv_accel_alloc( size_t n );

/* Acceleration free for evalutation of hessian
 *----------------------------------------------------------------------------*/
void ool_diff_Hv_accel_free( ool_diff_Hv_accel *a );

/* Eval numerically the hessian of f times a vector
 *----------------------------------------------------------------------------*/
int ool_diff_Hv( const ool_diff_Hv_accel *a,
		 void(*df)( const gsl_vector*, void*, gsl_vector* ),
		 const gsl_vector        *X,
		 void                    *fparam,
		 const gsl_vector        *V,
		 gsl_vector              *Hv,
		 const double             eps    );
/*----------------------------------------------------------------------------*/
__END_DECLS

#endif /*__OOL_DIFF_H__*/
/*----------------------------------------------------------------------------*/
