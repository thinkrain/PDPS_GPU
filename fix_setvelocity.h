/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(setvelocity,FixSetVelocity)

#else

#ifndef PS_FIX_SET_VELOCITY_H
#define PS_FIX_SET_VELOCITY_H

#include "fix.h"

namespace PDPS_NS {

class FixSetVelocity : public Fix {
public:
	FixSetVelocity(class PDPS *, int, char **);
	~FixSetVelocity();
	int setmask();
	void post_force();
	//void init();
	void setup();
	//double compute_vector(int);
	//double memory_usage();

 private:
	 double vx, vy, vz;
	int flag_x, flag_y, flag_z;
};

}

#endif
#endif
