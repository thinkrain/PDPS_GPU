/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_RANMARS_H
#define PS_RANMARS_H

#include "pointers.h"

namespace PDPS_NS {

class RanMars : protected Pointers {
public:
	RanMars(class PDPS *, int);
	~RanMars();
	double uniform();
	double gaussian();

private:
	int seed,save;
	double second;
	double *u;
	int i97,j97;
	double c,cd,cm;
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid seed for Marsaglia random # generator

The initial seed for this random number generator must be a positive
integer less than or equal to 900 million.

*/
