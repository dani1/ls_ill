/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * conmin/conmin_vector.c
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
 * since: Mar, 28, 2004
 *
 * $Id: conmin_vector.c,v 1.5 2005/05/12 23:16:13 biloti Exp $
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <conmin_vector.h>

/* A = max( B, min( C, D ) )
 *----------------------------------------------------------------------------*/
void conmin_vector_maxofmin( const size_t nn, 
			     double *A, double *B,
			     double *C, double *D )
{
  size_t ii;

  for( ii = 0; ii < nn; ii++ )
    A[ii] = GSL_MAX( B[ii], GSL_MIN( C[ii], D[ii] ) );
}

/* A = min( B, max( C, D ) )
 *----------------------------------------------------------------------------*/
void conmin_vector_minofmax( const size_t nn, 
			     double *A, double *B,
			     double *C, double *D )
{
  size_t ii;

  for( ii = 0; ii < nn; ii++ )
    A[ii] = GSL_MIN( B[ii], GSL_MAX( C[ii], D[ii] ) );
}

/* Destine = Source
 *----------------------------------------------------------------------------*/
void conmin_vector_memcpy( const size_t nn, double *D, double *S )
{
  const size_t size = nn * sizeof( double );

  memcpy( (void*)D, (void*)S, size );
}

/* X = v
 *----------------------------------------------------------------------------*/
void conmin_vector_set_all( const size_t nn, double *X, double vv )
{
  size_t ii;

  for( ii = 0; ii < nn; ii++ )
     X[ii] = vv;
}

/* X = 0
 *----------------------------------------------------------------------------*/
void conmin_vector_set_zero( const size_t nn, double *X )
{
  size_t ii;

  /* This was formely done with memset, but that was wrong since the
   * double zero must not be represent as zero bytes on any machine */
  for( ii = 0; ii < nn; ii++ )
     X[ii] = 0.0;
}

/* dist = | X - Y |_2
 *----------------------------------------------------------------------------*/
double conmin_vector_dist( const size_t nn, double *X, double *Y )
{
   double scale = 0.0;
   double ssq = 1.0;
   size_t ii;
   
   if (nn == 1) {
      return fabs ( X[0] - Y[0] );
   }
   
   for (ii = 0; ii < nn; ii++) {
      double x = X[ii] - Y[ii];

      
      if (x != 0.0) {
	 const double ax = fabs(x);
	 
	 if (scale < ax) {
	    ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
	    scale = ax;
	 } else {
	    ssq += (ax / scale) * (ax / scale);
	 }
      }
   }

   return scale * sqrt(ssq);
}

/* dist = | X - Y |_inf
 *----------------------------------------------------------------------------*/
double conmin_vector_dist_inf( const size_t nn, double *X, double *Y )
{
  size_t ii;
  double dist = 0;
  
  for( ii = 0; ii < nn; ii++ )
    dist = GSL_MAX( dist, fabs( X[ii] - Y[ii] ) );

  return dist;
}

/*----------------------------------------------------------------------------*/
