/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(setforce,FixSetForce)

#else

#ifndef PS_FIX_SET_FORCE_H
#define PS_FIX_SET_FORCE_H

#include "fix.h"

namespace PDPS_NS {

class FixSetForce : public Fix {
public:
	FixSetForce(class PDPS *, int, char **);
	~FixSetForce();
	int setmask();
	void post_force();
	//void init();
	void setup();
	//double compute_vector(int);
	//double memory_usage();

 private:
	double fx, fy, fz;
	int flag_x, flag_y, flag_z;
};

}

#endif
#endif
