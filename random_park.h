/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_RANPARK_H
#define PS_RANPARK_H

#include "pointers.h"

namespace PDPS_NS {

class RanPark : protected Pointers {
  friend class Set;
public:
	RanPark(class PDPS *, int);
	double uniform();
	double gaussian();
	void reset(int);
	void reset(int, double *);
	int state();
 
private:
	int seed,save;
	double second;
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid seed for Park random # generator

The initial seed for this random number generator must be a positive
integer.

*/
