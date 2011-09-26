/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * gencan/ool_conmin_gencan.h
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
 * Ricardo Biloti
 * since: Feb 16th, 2004
 *
 * $Id: ool_conmin_gencan.h,v 1.6 2005/05/17 19:00:51 biloti Exp $
 *----------------------------------------------------------------------------*/

#ifndef __OOL_CONMIN_GENCAN_H__
#define __OOL_CONMIN_GENCAN_H__

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

extern const ool_conmin_minimizer_type *ool_conmin_minimizer_gencan; 

/* gencan parameters struct
 *----------------------------------------------------------------------------*/
typedef struct
{
   
   double epsgpen;    /* default = 1.0e-05 */
   double epsgpsn;    /* default = 1.0e-05 */

   double fmin;       /* default =-1.0e+99 */

   double udelta0;    /* default = -1      */

   double ucgmia;     /* default = -1      */
   double ucgmib;     /* default = -1      */

   int    cg_scre;    /* default = 1       */
   double cg_gpnf;    /* default = epsgpen */
   double cg_epsi;    /* default = 1.0e-01 */ 
   double cg_epsf;    /* default = 1.0e-05 */

   double cg_epsnqmp; /* default = 1.0e-04 */
   size_t cg_maxitnqmp; /* default = 5     */

   int    nearlyq;    /* default = 0       */
   double nint;       /* default = 2.0     */
   double next;       /* default = 2.0     */
   size_t mininterp;  /* default = 4       */
   size_t maxextrap;  /* default = 100     */

   int    trtype;     /* default = 0       */

   double eta;        /* default = 0.9     */
   double delmin;     /* default = 0.1     */

   double lspgmi;     /* default = 1.0e-10 */
   double lspgma;     /* default = 1.0e+10 */

   double theta;      /* default = 1.0e-06 */
   double gamma;      /* default = 1.0e-04 */ 
   double beta;       /* default = 0.5     */

   double sigma1;     /* default = 0.1     */
   double sigma2;     /* default = 0.9     */

   double epsrel;     /* default = 1.0e-07 */
   double epsabs;     /* default = 1.0e-10 */

   double infrel;     /* default = 1.0e+20 */
   double infabs;     /* default = 1.0e+99 */


} ool_conmin_gencan_parameters;

__END_DECLS

#endif /*__OOL_CONMIN_GENCAN_H__*/

/*----------------------------------------------------------------------------*/
