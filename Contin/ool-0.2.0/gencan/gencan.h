/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * gencan/gencan.h
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
 * since: Apr 03rd, 2004
 *
 * $Id: gencan.h,v 1.10 2005/05/17 18:10:05 biloti Exp $
 *----------------------------------------------------------------------------*/

#ifndef __GENCAN_H__
#define __GENCAN_H__

#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include <ool/ool_conmin_common.h>
#include <ool/conmin_vector.h>
#include <ool/conmin_tools.h>
#include <ool/ool_conmin_gencan.h>

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

/* gencan state struct
 *----------------------------------------------------------------------------*/
typedef struct
{
   /* Space dimension */
   size_t n;
   
   /* Lower and Upper bounds */
   gsl_vector *L;
   gsl_vector *U;

   /* Numerical lower and upper bounds */
   double *near_l;
   double *near_u;

  /* Dimesion of reduced space */
  size_t nind;

  /* Index of reduced space */
  gsl_vector_uint *Ind;

  double xeucn;   /* Euclidean norm of X */ 
  double xsupn;   /* Infinite  norm of X */ 

  /* Working vectors */
  gsl_vector *S;
  gsl_vector *Y;
  gsl_vector *D;

  /* ometa2 = ( 1 - eta )^2 */
  double ometa2;

  /* epsgpen2 = epspgen^2 */
  double epsgpen2;

  /* Working variables */
  double lambda;
  double sts;
  double sty;
  double sinf;

  /* Reduced gradient related variables                   */
  double gpeucn2; /* Euclidean norm of "gpi"              */
  double gpsupn;  /* Infinite  norm of "gpi"              */
  double gieucn2; /* Internal gradient euclidean norm     */

  double acgeps;
  double bcgeps;
  double gpsupn0;
  double gpeucn20;

  /* working vector */
  gsl_vector *Xtrial;

  /* tnls state
   *--------------------------------------------------------------------------*/

  /* tnls working vectors */
  gsl_vector *tnls_Xtemp;

  /* Maximum step */
  double tnls_amax;

  /* Conjungate Gradient state 
   *--------------------------------------------------------------------------*/
  /* cg working vectors */
  gsl_vector *cg_W;
  gsl_vector *cg_R;
  gsl_vector *cg_D;
  gsl_vector *cg_Sprev;

  double cg_delta;

} conmin_gencan_state;

/* gencan functions
 *----------------------------------------------------------------------------*/
int gencan_prepare_iteration( ool_conmin_minimizer *M );

void gencan_projected_gradient( conmin_gencan_state *st,
				gsl_vector          *X,
				gsl_vector          *gradient );

int gencan_are_close( size_t  nn, 
		      double  alpha,
		      double *d, 
		      double *x, 
		      double  epsrel,
		      double  epsabs  );

int gencan_actual_iterate( ool_conmin_minimizer         *M, 
			   conmin_gencan_state          *st,
			   ool_conmin_gencan_parameters *P );

int gencan_spgls( ool_conmin_minimizer         *M, 
		  conmin_gencan_state          *st,
		  ool_conmin_gencan_parameters *P   );

int gencan_tnls( ool_conmin_minimizer         *M, 
		 conmin_gencan_state          *st,
		 ool_conmin_gencan_parameters *P  );

/* cg declarations
 *----------------------------------------------------------------------------*/

/* Return values */
#define GENCAN_CG_TRUST_REGION   1101 /* Achieved boundary of trust region    */
#define GENCAN_CG_BOUNDARY       1102 /* Convergence to the boundary          */
#define GENCAN_CG_S              1104 /* Stopping with s = s_k  such that 
				       * <g,s_k> <=-theta|g|_2 |s_k|_2 and 
				       * <g,s_{k+1}> > -theta|g|_2|s_{k+1}|_2 */
#define GENCAN_CG_CLOSE_ITERATE  1105 /* Very similar consecutive iterates    */
#define GENCAN_CG_INSUF_PROG     1106 /* Insufficent progress on the quadratic
					 model during cg_maxitnqmp iterations */
#define GENCAN_TNLS_MAXEXTRAP    1107 /* Maximum number of extrapolations
					 exceeded */


int gencan_cg( ool_conmin_minimizer         *M, 
	       conmin_gencan_state          *st,
	       ool_conmin_gencan_parameters *P  );

/*----------------------------------------------------------------------------*/
__END_DECLS

#endif /*__OOL_CONMIN_GENCAN_H__*/

/*----------------------------------------------------------------------------*/
