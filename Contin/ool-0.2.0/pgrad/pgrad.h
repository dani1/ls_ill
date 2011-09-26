/*----------------------------------------------------------------------------*
 * Open Optimization Library - Projected Gradient Method (header)
 * 
 * pgrad/pgrad.h
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
 * $Id: pgrad.h,v 1.8 2005/05/11 01:09:09 biloti Exp $
 *----------------------------------------------------------------------------*/

#ifndef __PGRAD_H__
#define __PGRAD_H__

#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include <ool/conmin_vector.h>
#include <ool/ool_conmin_common.h>
#include <ool/ool_conmin_pgrad.h>

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

/* pgrad state struct
 *----------------------------------------------------------------------------*/
typedef struct
{
  /* Space dimension */
  size_t n;

  /* Lower and Upper bounds */
  gsl_vector *L;
  gsl_vector *U;

  /* Working vectors */
  gsl_vector *xx;

} conmin_pgrad_state;

/* gradproj functions
 *----------------------------------------------------------------------------*/

static
void pgrad_line_search ( ool_conmin_minimizer *M );

void pgrad_proj( gsl_vector *L, gsl_vector *U, gsl_vector *X );
/*----------------------------------------------------------------------------*/
__END_DECLS

#endif /*__OOL_CONMIN_PGRAD_H__*/
/*----------------------------------------------------------------------------*/
