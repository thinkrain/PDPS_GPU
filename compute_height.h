/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(height, ComputeHeight)

#else

#ifndef PS_COMPUTE_HEIGHT_H
#define PS_COMPUTE_HEIGHT_H

#include "compute.h"

namespace PDPS_NS {

class ComputeHeight : public Compute {
public:
	ComputeHeight(class PDPS *, int, char **);
	virtual ~ComputeHeight();
	void init();
	//double compute_scalar();
	double compute_scalar();
	double findhigh();
	double createtop();

protected:

	double height, heightcut;
	double xtemp, ytemp, z0;
	double rho_ref;
	//void virial_compute(int);
};

}

#endif
#endif
