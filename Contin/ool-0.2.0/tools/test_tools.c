/*----------------------------------------------------------------------------*
 * Open Optimization Library - Example Programs
 * 
 * tools/test_tools.c
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
 * since: May, 23, 2004
 *
 * $Id: test_tools.c,v 1.4 2004/06/04 21:26:23 akiles Exp $
 *----------------------------------------------------------------------------*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <ool/ool_tools_diff.h>

typedef struct 
{
   gsl_matrix_view A;
   gsl_vector_view b;
   double c;

   double eps;

   /* working space */
   gsl_vector *Ax_plus_b;

} test_parameters;

ool_diff_Hv_accel *a;

/*----------------------------------------------------------------------------*/
static double test_f  ( const gsl_vector *X, void *vP )
{
   double f;
   test_parameters *P = (test_parameters*) vP;

   gsl_vector_memcpy( P->Ax_plus_b, &(P->b.vector) );
   gsl_blas_dgemv( CblasNoTrans, 0.5, &(P->A.matrix),
		   X, 1.0, P->Ax_plus_b );
   
   gsl_blas_ddot( X, P->Ax_plus_b, &f );
   
   return f + P->c;
}

static void test_df ( const gsl_vector *X, void* vP, gsl_vector* G )
{
   test_parameters *P = (test_parameters*) vP;

   gsl_vector_memcpy( G, &(P->b.vector) );
   gsl_blas_dgemv( CblasNoTrans, 1.0, &(P->A.matrix), X, 1.0, G );
}


static void test_fdf( const gsl_vector *X, void* vP, double *f, gsl_vector *G )
{
   (*f) = test_f( X, vP );

   test_df( X, vP, G );

}

static void test_Hv ( const gsl_vector *X, void *vP, 
		      const gsl_vector *V, gsl_vector *Hv )
{
   test_parameters *P = (test_parameters*) vP;

   gsl_blas_dgemv( CblasNoTrans, 1.0, &(P->A.matrix), V, 0.0, Hv);

}

/*----------------------------------------------------------------------------*/
static void test_ndf ( const gsl_vector *X, void* vP, gsl_vector *G )
{
   test_parameters *P = (test_parameters*) vP;
   
   ool_diff_g( &test_f, X, vP, G, P->eps );
}


static void test_ndf2 ( const gsl_vector *X, void* vP, gsl_vector *G )
{
   test_parameters *P = (test_parameters*) vP;
   
   ool_diff_g_auto( &test_f, X, vP, G );
}

static void test_nHv ( const gsl_vector *X, void *vP, 
		       const gsl_vector *V, gsl_vector *Hv )
{
  test_parameters *P = (test_parameters*) vP;

  const double eps = P->eps;

  ool_diff_Hv( a, &test_df, X, vP, V, Hv, eps );
}

static void test_nHv1 ( const gsl_vector *X, void *vP, 
			const gsl_vector *V, gsl_vector *Hv )
{
  test_parameters *P = (test_parameters*) vP;

  const double eps = P->eps;

  ool_diff_Hv( a, &test_ndf, X, vP, V, Hv, eps );
}

static void test_nHv2 ( const gsl_vector *X, void *vP, 
			const gsl_vector *V, gsl_vector *Hv )
{
  test_parameters *P = (test_parameters*) vP;

  const double eps = P->eps;

  ool_diff_Hv( a, &test_ndf2, X, vP, V, Hv, eps );
}

/*----------------------------------------------------------------------------*/
int main( void )
{
  size_t nn = 5;
  double A[] = { 55,  -3,  35,  33,  23,
		 -3,  70,  63,   6, -34,
		 35,  63, 167,  -3,  21,
		 33,   6,  -3,  84,   8,
		 23, -34,  21,   8, 154};
  double b[] = { -2,  10,   7,   3,   3};

  gsl_vector *X;
  gsl_vector *G1;
  gsl_vector *G2;
  gsl_vector *V;
  
  double f;
  size_t i;
  
  test_parameters P;
  
  /* Test function paramters */
  P.A = gsl_matrix_view_array( A, nn, nn );
  P.b = gsl_vector_view_array( b, nn );
  P.c = 2;

  P.eps = 2.0e-5;

  P.Ax_plus_b = gsl_vector_alloc( nn );

  /* Working vectors */
  X = gsl_vector_alloc( nn );
  G1 = gsl_vector_alloc( nn );
  G2 = gsl_vector_alloc( nn );
  V = gsl_vector_alloc( nn );

  for (i=0; i<nn; i++){
     gsl_vector_set( X, i, 2.0*i + 10.1 );
     gsl_vector_set( V, i, 3.0*i -  5.1 );
  }

  /* Evaluating the test function */
  f = test_f( X, (void*)(&P) );
  printf( "Function value: %.6e\n", f);

  /* Analytical gradient */
  test_df( X, (void*)(&P), G2 );
  printf( "Analytical gradient:\n" );
  gsl_vector_fprintf( stdout, G2, "% 12.4f" );

  /* Numerical gradient */
  test_ndf( X, (void*)(&P), G1 );
  printf( "\nDifference between analytical and numerical gradient:\n"
	  "(with choosen difference step)\n");
  gsl_vector_sub( G1, G2 );
  gsl_vector_fprintf( stdout, G1, "% 12.4e" );
  printf( "Relative error = %6.4e\n\n",
	  gsl_blas_dnrm2(G1)/gsl_blas_dnrm2(G2) );

  /* Numerical gradient */
  test_ndf2( X, (void*)(&P), G1 );
  printf( "Difference between analytical and numerical gradient):\n"
	  "(with automaticaly set difference step)\n");
  gsl_vector_sub( G1, G2 );
  gsl_vector_fprintf( stdout, G1, "% 12.4e" );
  printf( "Relative error = %6.4e\n\n",
	  gsl_blas_dnrm2(G1)/gsl_blas_dnrm2(G2) );

  /* Hessian test */

  a = ool_diff_Hv_accel_alloc( nn );

  test_Hv(  X, (void*)(&P), V, G2 );
  fprintf( stdout, "Analytical Hessian times vector:\n" );
  gsl_vector_fprintf( stdout, G2, "%12.4e" );

  test_nHv(  X, (void*)(&P), V, G1 );
  printf( "\nDifference between analytical and numerical hessian times vector:\n" 
	  "(from analytical gradient)\n");
  gsl_vector_sub( G1, G2 );
  gsl_vector_fprintf( stdout, G1, "%12.4e");
  printf( "Relative error = %6.4e\n\n",
	  gsl_blas_dnrm2(G1)/gsl_blas_dnrm2(G2) );

  test_nHv1(  X, (void*)(&P), V, G1 );
  printf( "Difference between analytical and numerical hessian times vector:\n"
	  "(from numerical gradient with choosen difference step)\n");
  gsl_vector_sub( G1, G2 );
  gsl_vector_fprintf( stdout, G1, "%12.4e");
  printf( "Relative error = %6.4e\n\n", 
	  gsl_blas_dnrm2(G1)/gsl_blas_dnrm2(G2) );

  test_nHv2(  X, (void*)(&P), V, G1 );
  printf( "Difference between analytical and numerical hessian times vector:\n"
	  "(from numerical gradient with automatically set difference step)\n" );
  gsl_vector_sub( G1, G2 );
  gsl_vector_fprintf( stdout, G1, "%12.4e");
  printf( "Relative error = %6.4e\n\n",
	  gsl_blas_dnrm2(G1)/gsl_blas_dnrm2(G2) );

  ool_diff_Hv_accel_free( a );

  gsl_vector_free( P.Ax_plus_b );
  gsl_vector_free( V );
  gsl_vector_free( G1 );
  gsl_vector_free( G2 );
  gsl_vector_free( X );

  return GSL_SUCCESS;
}

/*----------------------------------------------------------------------------*/
