/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"

#include "error.h"
#include "psmath.h"

using namespace PsMath_NS;

#define EPSILON 1.0e-6

/* ----------------------------------------------------------------------
   Calculate inverse matrix
   http://en.wikipedia.org/wiki/Invertible_matrix
---------------------------------------------------------------------- */
int PsMath_NS::Matrix_Inverse_3D(double inv_M[3][3], double M[3][3]) 
{
	double a, b, c, d, e, f, g, h, i;
	double A, B, C, D, E, F, G, H, I;	

	a = M[0][0];
	b = M[0][1];
	c = M[0][2];
	d = M[1][0];
	e = M[1][1];
	f = M[1][2];
	g = M[2][0];
	h = M[2][1];
	i = M[2][2];

	double detA = Matrix_Det_3D(M);

	if (fabs(detA) < EPSILON) return 0;

	double inv_detA = 1.0 / detA;

	A = e*i - f*h;
	D = -(b*i - c*h);
	G = (b*f - c*e);
	B = -(d*i - f*g);
	E = (a*i - c*g);
	H = -(a*f - c*d);
	C = (d*h - e*g);
	F = -(a*h - b*g);
	I = (a*e - b*d);
	
	inv_M[0][0] = A * inv_detA;
	inv_M[0][1] = D * inv_detA;
	inv_M[0][2] = G * inv_detA;
	inv_M[1][0] = B * inv_detA;
	inv_M[1][1] = E * inv_detA;
	inv_M[1][2] = H * inv_detA;
	inv_M[2][0] = C * inv_detA;
	inv_M[2][1] = F * inv_detA;
	inv_M[2][2] = I * inv_detA;

	return 1;
}

/* ----------------------------------------------------------------------
   Calculate determinant of a matrix
   http://en.wikipedia.org/wiki/Invertible_matrix
---------------------------------------------------------------------- */

double PsMath_NS::Matrix_Det_3D(double M[3][3]) 
{
	double det;	
	double a, b, c, d, e, f, g, h, i;
	
	a = M[0][0];
	b = M[0][1];
	c = M[0][2];
	d = M[1][0];
	e = M[1][1];
	f = M[1][2];
	g = M[2][0];
	h = M[2][1];
	i = M[2][2];

	det = a*(e*i - f*h) - b*(i*d - f*g) + c*(d*h - e*g);

	return det;
}

/* ----------------------------------------------------------------------
   Matrix A (3*3) * Matrix B (3*3) 
   http://en.wikipedia.org/wiki/Matrix_multiplication
---------------------------------------------------------------------- */

void PsMath_NS::Matrix_Prod_3D(double AB[3][3], double A[3][3], double B[3][3]) 
{
	int i, j, k, m, n, p;

	m = n = p = 3;

	double temp;
	for (i = 0; i < n; i++)
	for (j = 0; j < p; j++) {
		temp = 0.0;
		for (k = 0; k < m; k++) {
			temp += A[i][k]*B[k][j];
		}
		AB[i][j] = temp;
	}
}

/* ----------------------------------------------------------------------
   Matrix A (3*3) * Matrix B (3*1) 
   http://en.wikipedia.org/wiki/Matrix_multiplication
---------------------------------------------------------------------- */

void PsMath_NS::Matrix_Prod_3D(double AB[3], double A[3][3], double B[3]) 
{
	int i, j, k, m, n, p;

	m = n = 3;

	double temp;
	for (i = 0; i < n; i++) {
		temp = 0.0;
		for (k = 0; k < m; k++) {
			temp += A[i][k]*B[k];
		}
		AB[i] = temp;
	}
}

/* ---------------------------------------------------------------------- */

double PsMath_NS::Vec_Norm2(double A[3]) 
{
	double norm;

	norm = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
	norm = sqrt(norm);

	return norm;
}

/* ----------------------------------------------------------------------
   Vector A cross_prod Vector B  
   http://en.wikipedia.org/wiki/Cross_product
---------------------------------------------------------------------- */

double PsMath_NS::Vec_Dot_Prod_3D(double A[3], double B[3]) 
{
	double ans;

	ans = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
	
	return ans;
}

/* ----------------------------------------------------------------------
   Vector A cross_prod Vector B  
   http://en.wikipedia.org/wiki/Cross_product
---------------------------------------------------------------------- */

void PsMath_NS::Vec_Cross_Prod_3D(double C[3], double A[3], double B[3]) 
{
	
	C[0] = A[1]*B[2] - A[2]*B[1];
	C[1] = A[2]*B[0] - A[0]*B[2];
	C[2] = A[0]*B[1] - A[1]*B[0];
}

/* ---------------------------------------------------------------------- */

void PsMath_NS::Matrix_Print_3D(double A[3][3], const char *str)
{
	fprintf(stdout, "Matrix %s = \n", str);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			fprintf(stdout, "%f ", A[i][j]);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "\n");		 
}

/* ---------------------------------------------------------------------- */

void PsMath_NS::Vec_Print_3D(double A[3], const char *str) 
{
	fprintf(stdout, "Vector %s = ", str);
	for (int i = 0; i < 3; i++) {
		fprintf(stdout, "%f ", A[i]); 
	}
	fprintf(stdout, "\n");
}

