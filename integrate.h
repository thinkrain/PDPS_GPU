/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_INTEGRATE_H
#define PS_INTEGRATE_H

#include "pointers.h"

namespace PDPS_NS {

class Integrate : protected Pointers {
public:
	char *style;

	Integrate(class PDPS *, int , char **);
	virtual ~Integrate();
	void init();
	virtual void setup() = 0;
	virtual void run(int) = 0;

protected:
	int eflag, vflag;                             // flags for energy/virial computation
	int virial_style;                             // compute virial explicityly or implicitly

	void ev_setup();
	void ev_set(bigint);

};

}

#endif
