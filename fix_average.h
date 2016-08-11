/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator/*
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(average,FixAverage)

#else

#ifndef PS_FIX_AVERAGE_H
#define PS_FIX_AVERAGE_H

#include "fix.h"

namespace PDPS_NS {

class FixAverage : public Fix {

public:
	FixAverage(class PDPS *, int, char **);
	~FixAverage();
	int setmask();
	void init();
	void setup();
	virtual void post_force();
	
protected:
	int gid;                  // group id
	int force_flag, force_x, force_y, force_z;
};
}

#endif
#endif
