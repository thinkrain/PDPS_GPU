/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(volume,FixVolume)

#else

#ifndef PS_FIX_VOLUME_H
#define PS_FIX_V0LUME_H

#include "fix.h"

namespace PDPS_NS {

class FixVolume : public Fix {
public:
	FixVolume(class PDPS *, int, char **);
	~FixVolume();
	int setmask();
	void post_force();
	//void init();
	void setup();
	//double compute_vector(int);
	//double memory_usage();
	int gravity_flag, initial_flag, region_flag, temperature_flag;
	int g_id, gasbit, rid;
	double radius_origin;
	
	

 private:
	
	double ratio;
	double xacc;
	double yacc;
	double zacc;
	double expansion;
	int count;
	double T_liq, T_boil, Latent, rho_liq, rho_bub;

};

}

#endif
#endif
