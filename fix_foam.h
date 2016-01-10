/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(foam,FixFoam)

#else

#ifndef PS_FIX_FOAM_H
#define PS_FIX_FOAM_H

#include "fix.h"

namespace PDPS_NS {

class FixFoam : public Fix {

public:
	FixFoam(class PDPS *, int, char **);
	~FixFoam();
	int setmask();
	void init();
	void setup();
	virtual void post_force();
	
protected:
	int rid, tid, newgid, refgid;                  // region id, type id, new group id, reference group id
	int region_flag;          // region enable flag
	double frequency;
	int direction;
	double seed;
	double gap;
	int count;
	double radius_bubble;
	double mass_bubble;
	int level_flag, neighbor_flag;
	

};
}

#endif
#endif
