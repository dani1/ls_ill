/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * spg/ool_conmin_spg.h
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
 * $Id: ool_conmin_spg.h,v 1.2 2005/05/09 20:19:41 biloti Exp $
 *----------------------------------------------------------------------------*/

#ifndef __OOL_CONMIN_SPG_H__
#define __OOL_CONMIN_SPG_H__

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

extern const ool_conmin_minimizer_type *ool_conmin_minimizer_spg;

/* spg parameters struct
 *----------------------------------------------------------------------------*/
typedef struct
{

   double fmin;     /* default = -1e+99 */
   double tol;      /* default =  1e-4  */

   size_t M;        /* default = 1      */
   
   double alphamin; /* default = 1e-30  */
   double alphamax; /* default = 1e+30  */
   
   double gamma;    /* default = 1e-4   */
   
   double sigma1;   /* default = 0.1    */
   double sigma2;   /* default = 0.9    */
   
} ool_conmin_spg_parameters;

__END_DECLS

#endif /*__OOL_CONMIN_SPG_H__*/

/*----------------------------------------------------------------------------*/
