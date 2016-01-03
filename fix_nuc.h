/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nuc,FixNuc)

#else

#ifndef PS_FIX_NUC_H
#define PS_FIX_NUC_H

#include "fix.h"

namespace PDPS_NS {

class FixNuc : public Fix {

public:
	FixNuc(class PDPS *, int, char **);
	~FixNuc();
	int setmask();
	void init();
	void setup();
	virtual void post_force();
	
protected:
	int rid, tid, newgid;                  // region id, type id, new group id
	int region_flag;          // region enable flag
	double frequency;
	int direction;
	double seed;
	double gap;
	int count;
	double radius_bubble;
	double mass_bubble;

};
}

#endif
#endif
