/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(andersen, FixAndersen)

#else

#ifndef PS_FIX_ANDERSON_H
#define PS_FIX_ANDERSON_H

#include "fix.h"

namespace PDPS_NS {

class FixAndersen : public Fix {
public:
	FixAndersen(class PDPS *, int, char **);
	~FixAndersen();
	int setmask();
	virtual void init();
	virtual void end_of_step();

protected:
	int **setflag;
	double temperature;
	class Pair *pair;
	double probability, Gamma;
	double muij, kbT;
	int seed;
	class RanMars *random;
};

}

#endif
#endif
