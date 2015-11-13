/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nve/sphere, FixNVESphere)

#else

#ifndef PS_FIX_NVE_SPHERE_H
#define PS_FIX_NVE_SPHERE_H

#include "fix_nve.h"

namespace PDPS_NS {

class FixNVESphere : public FixNVE {
public:
	FixNVESphere(class PDPS *, int, char **);
	virtual ~FixNVESphere() {}
	int setmask();
	virtual void init();
	virtual void initial_integrate();
	virtual void final_integrate();
	//virtual void reset_dt();

protected:
	//double dtv, dtf;
	int mass_require;
};

}

#endif
#endif
