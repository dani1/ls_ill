/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * conmin/conmin_tools.h
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
 * $Id: conmin_tools.h,v 1.4 2005/04/15 13:46:11 biloti Exp $
 *----------------------------------------------------------------------------*/

#ifndef __CONMIN_TOOLS_H__
#define __CONMIN_TOOLS_H__

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

/*----------------------------------------------------------------------------*/
/* Shrink vector V from the full space (n) to the reduced space (nind). */
void conmin_shrink( const size_t nind, gsl_vector_uint *Ind, gsl_vector *V );

/* Expands vector v from the reduced space (nind) to the full space (n). */
void conmin_expand( const size_t nind, gsl_vector_uint *Ind, gsl_vector *V );

/* Evaluate objective function from the reduced space */
double conmin_calc_f( ool_conmin_minimizer *M, 
		      const size_t          nind, 
		      gsl_vector_uint      *Ind,
		      gsl_vector           *X,
		      gsl_vector           *Xc );

/* Evaluate gradient from the reduced space */
void conmin_calc_g( ool_conmin_minimizer *M, 
		    const size_t          nind, 
		    gsl_vector_uint      *Ind,
		    gsl_vector           *X,
		    gsl_vector           *Xc,
		    gsl_vector           *G  );

/* Evaluate hessian times a vector from the reduced space */
void conmin_calc_Hv( ool_conmin_minimizer *M, 
		     const size_t          nind, 
		     gsl_vector_uint      *Ind,
		     gsl_vector           *X,
		     gsl_vector           *Xc,
		     gsl_vector           *V,
		     gsl_vector           *Hv   );

__END_DECLS

#endif /*__CONMIN_TOOLS_H__*/
/*----------------------------------------------------------------------------*/



