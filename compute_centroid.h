/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(centroid, ComputeCentroid)

#else

#ifndef PS_COMPUTE_CENTROID_H
#define PS_COMPUTE_CENTROID_H

#include "compute.h"

namespace PDPS_NS {

class ComputeCentroid : public Compute {
public:
	ComputeCentroid(class PDPS *, int, char **);
	virtual ~ComputeCentroid();
	void init();
	//double compute_scalar();
	void compute_vector();

protected:

	double centroid[3];

	//void virial_compute(int);
};

}

#endif
#endif
