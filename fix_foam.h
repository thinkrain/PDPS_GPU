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
	int rid, tid, newgid, refgid;           // region id, type id, new group id, reference group id
	int level_flag, neighbor_flag;			//	The way to delete particle, based on liquid level or neighbor particle numbers
	

};
}

#endif
#endif
