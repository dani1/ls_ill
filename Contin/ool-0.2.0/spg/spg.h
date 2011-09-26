/*----------------------------------------------------------------------------*
 * Open Optimization Library - Spectral Projected Gradient Method (header)
 * 
 * spg/spg.h
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
 * $Id: spg.h,v 1.3 2005/05/10 20:24:27 biloti Exp $
 *----------------------------------------------------------------------------*/

#ifndef __SPG_H__
#define __SPG_H__

#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

#include <ool/conmin_vector.h>
#include <ool/ool_conmin_common.h>
#include <ool/ool_conmin_spg.h>

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

/* spg state struct
 *----------------------------------------------------------------------------*/
typedef struct
{
  /* Space dimension */
   size_t n;

   /* non-monotone parameter */
   size_t M; 
   size_t m;

   /* Lower and Upper bounds */
   gsl_vector *L;
   gsl_vector *U;
   
   /* Armijo parameter */
   double alpha;
   
   int tail;
   /* Working vectors */
   double *f;
   
   gsl_vector *xx;
   gsl_vector *d;
   gsl_vector *s;
   gsl_vector *y;
   
} conmin_spg_state;

/* spg functions
 *----------------------------------------------------------------------------*/

static
void spg_line_search ( ool_conmin_minimizer *M );

void spg_proj( gsl_vector *L,
	       gsl_vector *U,
	       gsl_vector *x );

/*----------------------------------------------------------------------------*/
__END_DECLS

#endif /*__OOL_CONMIN_SPG_H__*/
/*----------------------------------------------------------------------------*/
