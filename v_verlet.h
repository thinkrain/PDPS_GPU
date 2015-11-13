/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef INTEGRATE_CLASS

IntegrateStyle(v_verlet,V_Verlet)

#else

#ifndef PS_V_VERLET_H
#define PS_V_VERLET_H

#include "integrate.h"

namespace PDPS_NS {

class V_Verlet : public Integrate {
public:
	V_Verlet(class PDPS *, int, char **);
	~V_Verlet() {}
	void init();
	void setup();
	//void setup_minimal(int);
	void run(int);
	//void cleanup();

protected:
	int triclinic;                    // 0 if domain is orthog, 1 if triclinic
	int torqueflag,erforceflag;
	int e_flag,rho_flag;

	void force_clear();

private:
	int procid;
};

}

#endif
#endif
