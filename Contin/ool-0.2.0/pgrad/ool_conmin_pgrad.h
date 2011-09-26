/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * pgrad/ool_conmin_pgrad.h
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
 * since: Jul, 02, 2004
 *
 * $Id: ool_conmin_pgrad.h,v 1.4 2005/05/18 22:38:11 damxtha Exp $
 *----------------------------------------------------------------------------*/

#ifndef __OOL_CONMIN_PGRAD_H__
#define __OOL_CONMIN_PGRAD_H__

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

extern const ool_conmin_minimizer_type *ool_conmin_minimizer_pgrad;

/* pgrad parameters struct
 *----------------------------------------------------------------------------*/
typedef struct
{
  
  double fmin;   /* default = -1e+99 */
  double tol;    /* default =  1e-4  */
  double alpha;  /* default =  1e-4  */
 
  double sigma1;   /* default = 0.1    */
  double sigma2;   /* default = 0.9    */

} ool_conmin_pgrad_parameters;

__END_DECLS

#endif /*__OOL_CONMIN_PGRAD_H__*/

/*----------------------------------------------------------------------------*/
