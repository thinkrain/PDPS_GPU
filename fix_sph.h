/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(sph,FixSPH)

#else

#ifndef PS_FIX_SPH_H
#define PS_FIX_SPH_H

#include "fix.h"

namespace PDPS_NS {

class FixSPH : public Fix {
public:
	FixSPH(class PDPS *, int, char **);
	virtual ~FixSPH() {}
	int setmask();
	virtual void init();
	virtual void initial_integrate();
	virtual void final_integrate();
	virtual void setup_pre_force(int);


private:
	class NeighList *list;
protected:
	double dtv,dtf;
	int sphere_flag;

};

}

#endif
#endif
