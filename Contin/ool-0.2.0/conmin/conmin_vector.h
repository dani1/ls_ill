/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * conmin/conmin_vector.h
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
 * since: Apr, 09, 2004
 *
 * $Id: conmin_vector.h,v 1.3 2005/04/20 19:38:49 biloti Exp $
 *----------------------------------------------------------------------------*/

#ifndef __CONMIN_VECTOR_H__
#define __CONMIN_VECTOR_H__

#include <ool/ool_conmin.h>

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
/*----------------------------------------------------------------------------*/

/* A = max( B, min( C, D ) ) 
 * Replace by inner_point */
void conmin_vector_maxofmin( const size_t nn, 
			     double *A, double *B,
			     double *C, double *D );

/* A = min( B, max( C, D ) ) */
void conmin_vector_minofmax( const size_t nn, 
			     double *A, double *B,
			     double *C, double *D );

/* Destine = Source */
void conmin_vector_memcpy( const size_t nn, double *D, double *S );

/* X = v */
void conmin_vector_set_all( const size_t nn, double *X, double vv );

/* X = 0 */
void conmin_vector_set_zero( const size_t nn, double *X );

/* dist = | X - Y |_2 */
double conmin_vector_dist( const size_t nn, double *X, double *Y );

/* dist = | X - Y |_inf */
double conmin_vector_dist_inf( const size_t nn, double *X, double *Y );

/*----------------------------------------------------------------------------*/
__END_DECLS

#endif /*__CONMIN_VECTOR_H__*/
/*----------------------------------------------------------------------------*/
