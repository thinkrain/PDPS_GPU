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
	int level_flag, neighbor_flag, pressure_flag;			//	The way to delete particle, based on liquid level or neighbor particle numbers
	double radius_initial;					//	inital value of radius, to prevent early deletion
	int neighbor_delete;					//  minimun number of neighbor particles to delete
	double frequency, frequency0, rho_ref, v_ref;				//  frequency to do deletion, reference density, reference velocity, for pressure_flag
	int step_last, step_next;				//	record the last and next step to fix

};
}

#endif
#endif
