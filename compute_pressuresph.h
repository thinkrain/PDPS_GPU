/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(pressuresph, ComputePressuresph)

#else

#ifndef PS_COMPUTE_PRESSURESPH_H
#define PS_COMPUTE_PRESSURESPH_H

#include "compute.h"

namespace PDPS_NS {

class ComputePressuresph : public Compute {
public:
	ComputePressuresph(class PDPS *, int, char **);
	virtual ~ComputePressuresph();
	void init();
	//double compute_scalar();
	double compute_scalar();
	double findpressure();


protected:

	int cubic_flag, quintic_flag;
	double pressuresph, h, soundspeed, rho_ref;
	double xtemp, ytemp, ztemp;
	double a2D, a3D, wfsum;
	//void virial_compute(int);
};

}

#endif
#endif
