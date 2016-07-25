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
	double compute_scalar();
	double findhigh();

protected:

	double height, heightcut;		//  cutoff length to compute height
	double xtemp, ytemp, z0;		//  position to compute height
	double rho_ref;					//  reference rho at the limit to compute height

};

}

#endif
#endif
