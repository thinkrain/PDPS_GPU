/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(level, ComputeLevel)

#else

#ifndef PS_COMPUTE_LEVEL_H
#define PS_COMPUTE_LEVEL_H

#include "compute.h"

namespace PDPS_NS {

class ComputeLevel : public Compute {
public:
	ComputeLevel(class PDPS *, int, char **);
	virtual ~ComputeLevel();
	void init();
	//double compute_scalar();
	double compute_scalar();
	double findhigh();
	double createtop();

protected:

	int frequency;
	double level, level_pre, levelgap;

	//void virial_compute(int);
};

}

#endif
#endif
