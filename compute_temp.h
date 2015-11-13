/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(temp,ComputeTemp)

#else

#ifndef PS_COMPUTE_TEMP_H
#define PS_COMPUTE_TEMP_H

#include "compute.h"

namespace PDPS_NS {

class ComputeTemp : public Compute {
public:
	ComputeTemp(class PDPS *, int, char **);
	virtual ~ComputeTemp();
	void init();
	double compute_scalar();
	//void compute_vector();

protected:
	int fix_dof;
	double tfactor;

	virtual void dof_compute();
};

}

#endif
#endif
