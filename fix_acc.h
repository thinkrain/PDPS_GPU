/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(acc,FixAcc)

#else

#ifndef PS_FIX_ACC_H
#define PS_FIX_ACC_H

#include "fix.h"

namespace PDPS_NS {

class FixAcc : public Fix {

public:
	FixAcc(class PDPS *, int, char **);
	~FixAcc();
	int setmask();
	void init();
	void setup();
	virtual void post_force();
	
protected:
	double xacc, yacc, zacc;
	int rid;                  // region id
	int region_flag;          // region enable flag
};
}

#endif
#endif
