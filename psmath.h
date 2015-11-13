/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_PSMATH_H
#define PS_PSMATH_H

namespace PsMath_NS {
	// det = |A[3][3]|
	double Matrix_Det_3D(double [3][3]);

	// inv_A[3][3] = (A[3][3])^-1
	int Matrix_Inverse_3D(double [3][3], double [3][3]);

	// C[3][3] = A[3][3] * B[3][3]
	void Matrix_Prod_3D(double [3][3], double [3][3], double [3][3]);

	// C[3] = A[3][3] * B[3]
	void Matrix_Prod_3D(double [3], double [3][3], double [3]);

	// |A[3]| norm2
	double Vec_Norm2(double [3]);

	// C[3] = A[3] corss_prod B[3]
	void Vec_Cross_Prod_3D(double [3], double [3], double [3]);

	double Vec_Dot_Prod_3D(double [3], double [3]);

    //--------------- debug purpose --------------------
	// A[3][3] and its name "A"
	void Matrix_Print_3D(double [3][3], const char [3]);

	// A[3] and its name "A"
	void Vec_Print_3D(double [3], const char [3]);

	
}

#endif
