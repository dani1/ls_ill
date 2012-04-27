/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * conmin/ool_conmin_tests.h
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
 * since: Apr, 04, 2004
 *
 * $Id: ool_conmin_tests.h,v 1.1.1.1 2004/05/18 17:25:43 akiles Exp $
 *----------------------------------------------------------------------------*/

#ifndef __OOL_CONMIN_TESTS_H__
#define __OOL_CONMIN_TESTS_H__

#include <ool/ool_conmin_common.h>

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

/* Feasiblility tests
 *----------------------------------------------------------------------------*/

/* Test if x if feasible */
int ool_conmin_is_feasible( const ool_conmin_constraint *C, 
			    const gsl_vector            *X );

/* Test if L <= x <= u */
int ool_conmin_is_in_box( const ool_conmin_constraint *C, 
			  const gsl_vector            *X );

/* Test if g(x) <= 0 
 * Would be good to choose a better name!
 */
int ool_conmin_is_inequality_fasible( const ool_conmin_constraint *C, 
				      const gsl_vector            *X );

/* Test if h(x) = 0 ( indeed if |h(x)| < eq_tol )
 * Would be good to choose a better name! 
 */
int ool_conmin_is_equality_fasible( const ool_conmin_constraint *C, 
				    const gsl_vector            *X );

/*----------------------------------------------------------------------------*/

__END_DECLS

#endif /*__OOL_CONMIN_TESTS_H__*/

/*----------------------------------------------------------------------------*/



