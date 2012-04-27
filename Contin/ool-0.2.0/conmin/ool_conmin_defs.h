/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * conmin/ool_conmin_defs.h
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
 * since: Apr, 27, 2004
 *
 * $Id: ool_conmin_defs.h,v 1.3 2005/05/07 19:52:30 biloti Exp $
 *----------------------------------------------------------------------------*/

#ifndef __OOL_CONMIN_DEFS_H__
#define __OOL_CONMIN_DEFS_H__

#include <gsl/gsl_errno.h>

enum {
   /* OOL error codes from GSL error codes */
   OOL_SUCCESS    = GSL_SUCCESS  ,  /* Success */
   OOL_FAILURE    = GSL_FAILURE  ,
   OOL_CONTINUE   = GSL_CONTINUE ,  /* Iterations has not converged yet */ 
   OOL_EDOM       = GSL_EDOM     ,  /* input domain error, e.g sqrt(-1) */
   OOL_ERANGE     = GSL_ERANGE   ,  /* output range error, e.g. exp(1e100) */
   OOL_EFAULT     = GSL_EFAULT   ,  /* invalid pointer */
   OOL_EINVAL     = GSL_EINVAL   ,  /* invalid argument supplied by user */
   OOL_EFAILED    = GSL_EFAILED  ,  /* generic failure */
   OOL_EFACTOR    = GSL_EFACTOR  ,  /* factorization failed */
   OOL_ESANITY    = GSL_ESANITY  ,  /* sanity check failed - shouldn't happen */
   OOL_ENOMEM     = GSL_ENOMEM   ,  /* malloc failed */
   OOL_EBADFUNC   = GSL_EBADFUNC ,  /* problem with user-supplied function */
   OOL_ERUNAWAY   = GSL_ERUNAWAY ,  /* iterative process is out of control */
   OOL_EMAXITER   = GSL_EMAXITER ,  /* exceeded max number of iterations */
   OOL_EZERODIV   = GSL_EZERODIV ,  /* tried to divide by zero */
   OOL_EBADTOL    = GSL_EBADTOL  ,  /* user specified an invalid tolerance */
   OOL_ETOL       = GSL_ETOL     ,  /* failed to reach the specified tolerance */
   OOL_EUNDRFLW   = GSL_EUNDRFLW ,  /* underflow */
   OOL_EOVRFLW    = GSL_EOVRFLW  ,  /* overflow  */
   OOL_ELOSS      = GSL_ELOSS    ,  /* loss of accuracy */
   OOL_EROUND     = GSL_EROUND   ,  /* failed because of roundoff error */
   OOL_EBADLEN    = GSL_EBADLEN  ,  /* matrix, vector lengths are not conformant */
   OOL_ENOTSQR    = GSL_ENOTSQR  ,  /* matrix not square */
   OOL_ESING      = GSL_ESING    ,  /* apparent singularity detected */
   OOL_EDIVERGE   = GSL_EDIVERGE ,  /* integral or series is divergent */
   OOL_EUNSUP     = GSL_EUNSUP   ,  /* requested feature is not supported by the hardware */
   OOL_EUNIMPL    = GSL_EUNIMPL  ,  /* requested feature not (yet) implemented */
   OOL_ECACHE     = GSL_ECACHE   ,  /* cache limit exceeded */
   OOL_ETABLE     = GSL_ETABLE   ,  /* table limit exceeded */
   OOL_ENOPROG    = GSL_ENOPROG  ,  /* iteration is not making progress towards solution */
   OOL_ENOPROGJ   = GSL_ENOPROGJ ,  /* jacobian evaluations are not improving the solution */
   OOL_ETOLF      = GSL_ETOLF    ,  /* cannot reach the specified tolerance in F */
   OOL_ETOLX      = GSL_ETOLX    ,  /* cannot reach the specified tolerance in X */
   OOL_ETOLG      = GSL_ETOLG    ,  /* cannot reach the specified tolerance in gradient */
   OOL_EOF        = GSL_EOF      ,  /* end of file */

   /* OOL-specific error codes */
   OOL_UNBOUNDEDF =  1101,  /* Lower unbounded function           */
   OOL_INFEASIBLE =  1102,  /* Infeasible point                   */
   OOL_FINNERIT   =  1103,  /* Too many inner iterations          */
   OOL_FLSEARCH   =  1104,  /* Line search failed                 */
   OOL_FDDIR      =  1105   /* Unable to find a descent direction */
};

#endif /*__OOL_CONMIN_DEFS_H__*/

/*----------------------------------------------------------------------------*/



