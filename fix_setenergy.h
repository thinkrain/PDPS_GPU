/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(setenergy,FixSetEnergy)

#else

#ifndef PS_FIX_SET_ENERGY_H
#define PS_FIX_SET_ENERGY_H

#include "fix.h"

namespace PDPS_NS {

class FixSetEnergy : public Fix {
public:
	FixSetEnergy(class PDPS *, int, char **);
	~FixSetEnergy();
	int setmask();
	void post_force();
	//void init();
	void setup();
	//double compute_vector(int);
	//double memory_usage();

 private:
	 double energy;
};

}

#endif
#endif
