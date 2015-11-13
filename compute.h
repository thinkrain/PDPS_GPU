/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_COMPUTE_H
#define PS_COMPUTE_H

#include "stdio.h"
#include "pointers.h"

namespace PDPS_NS {

class Compute : protected Pointers {
 public:
	char *name;
    char *style;
	int gid;
	int groupbit;

	double scalar;                        // computed global scalar
	double *vector;                       // computed global vector
	double **array;                       // computed global array

	// flag
	int invoked_flag;                     // non-zero if invoked or accessed this step, 0 if not
	int invoked_scalar;                   // last timestep on which compute_scalar() was invoked
	int invoked_vector;
	int invoked_array;

	int scalar_flag;                      // 1 = calculate scalar  0 = not calculate scalar
	int vector_flag; 
	int array_flag;

	int size_vector;
	int size_array_rows;
	int size_array_columns;

	double dof;                           // degrees-of-freedom for temperature

	int comm_forward;                     // size of forward communication (0 if none)
	int comm_reverse;                     // size of reverse communication (0 if none)

	Compute(class PDPS *, int, char **);
	virtual ~Compute();	

	virtual void init() = 0;
    virtual double compute_scalar() {return 0.0;}
    virtual void compute_vector() {}
	virtual void compute_array() {}
	//void compute_temp();
	void compute_ke();
	void compute_pe();

protected:
	//int extra_dof;
};

}

#endif

