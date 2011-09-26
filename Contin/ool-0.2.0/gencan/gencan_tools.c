/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * gencan/gencan_tools.c
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
 * $Id: gencan_tools.c,v 1.7 2005/05/13 22:08:00 biloti Exp $
 *----------------------------------------------------------------------------*/

#include<gencan.h>

/* Compute continuous projected gradient infinite and euclidian norm, internal
 * gradient euclidian-norm, and store in nind the number of free variables
 * and in array p its identifiers */
void gencan_projected_gradient( conmin_gencan_state *st,
				gsl_vector          *X,
				gsl_vector          *gradient )
{
  /* Direct access to vector data */
  size_t *ind = st->Ind->data;

  double *x = X->data;
  double *g = gradient->data;

  double *l = st->L->data;
  double *u = st->U->data;

  /* Internal variables */
  size_t ii;
  size_t nn;
  size_t nind;

  double gpsupn;
  double gpeucn2;
  double gieucn2;

  nn = X->size;

  nind    = 0;
  gpsupn  = 0;
  gpeucn2 = 0;
  gieucn2 = 0;

  for( ii = 0; ii < nn; ii++ )
    {
      double gpi;
      double gpiabs;
      double gpi2;
      double aux, auxmax;
      
      auxmax = GSL_MAX( *l, *x - *g );
      aux    = GSL_MIN( *u, auxmax  );

      gpi    = aux - *x;
      gpiabs = fabs( gpi );
      gpi2   = gsl_pow_2( gpi );

      gpsupn   = GSL_MAX( gpsupn, gpiabs );
      gpeucn2 += gpi2;

      if ( *x > *l && *x < *u ) 
	{
	  gieucn2  += gpi2;
	  ind[nind] = ii;
	  nind++;
	}

      x++;
      u++;
      l++;
      g++;
    }

  st->nind    = nind;
  st->gpsupn  = gpsupn;
  st->gpeucn2 = gpeucn2;
  st->gieucn2 = gieucn2;
}

/* Test if two consecutive iterates are too close
 *----------------------------------------------------------------------------*/
int gencan_are_close( size_t  nn, 
		      double  alpha,
		      double *d, 
		      double *x, 
		      double  epsrel,
		      double  epsabs  )
{
   size_t ii;
   double aux;

   int samep = 1;

   for ( ii = 0; ii < nn && samep; ii++){
      
      aux = epsrel * fabs( x[ii] );

      if ( fabs(alpha * d[ii] ) > GSL_MAX( aux, epsabs ) )
	 samep = 0;
   }
   
   return samep;
}

/*----------------------------------------------------------------------------*/
