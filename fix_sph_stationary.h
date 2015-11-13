/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(sph/stationary, FixSPH_STATIONARY)

#else

#ifndef PS_FIX_SPH_STATIONARY_H
#define PS_FIX_SPH_STATIONARY_H

#include "fix.h"

namespace PDPS_NS {

	class FixSPH_STATIONARY : public Fix {
	public:
		FixSPH_STATIONARY(class PDPS *, int, char **);
		virtual ~FixSPH_STATIONARY() {}
		int setmask();
		virtual void init();
		virtual void initial_integrate();
		virtual void final_integrate();
		//virtual void setup_pre_force(int);
		//virtual void reset_dt();

	private:
		class NeighList *list;
	protected:
		double dtv, dtf;
		int mass_require;

		class Pair *pair;
	};

}

#endif
#endif
