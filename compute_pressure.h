/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(pressure,ComputePressure)

#else

#ifndef PS_COMPUTE_PRESSURE_H
#define PS_COMPUTE_PRESSURE_H

#include "compute.h"

namespace PDPS_NS {

class ComputePressure : public Compute {
public:
	ComputePressure(class PDPS *, int, char **);
	virtual ~ComputePressure();
	void init();
	double compute_scalar();
	//void compute_vector();

protected:
	int dim;
	double boltz,nktv2p,volume;
	Compute *temperature;
	char *id_temp;
	int cid;
	double virial[6];

	void virial_compute(int);
};

}

#endif
#endif
