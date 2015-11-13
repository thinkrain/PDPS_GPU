/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/lbm, FixWallLBM)

#else

#ifndef PS_FIX_WALL_LBM_H
#define PS_FIX_WALL_LBM_H

#include "fix.h"

namespace PDPS_NS {

class FixWallLBM : public Fix {
public:
	FixWallLBM(class PDPS *, int, char **);
	virtual ~FixWallLBM();
	int setmask();
	void init();
	void setup();
	void pre_force();

protected:
	int rot_flag;             // by default, rotation is considered
	
	double sij_min;           // minimum separation distance
	double phi;
	double gamma;
	double Vpcb, ln_Vpcb, Vpcb_cube_root;
	

	double cut, rneigh;                           // for particle-wall interaction cutoff and rneigh

	int nwalls;                                   // number of walls
	int wall[6];                                  // 0: xlo; 1: xhi .... 5: zhi
	int wall_style;
	double coords0[6];                            // store box dimension

private:
	int *wall_rid;                                // region id for the wall
	int *wall_flag;                               // flag for wall 

	int drag_flag;                                // drag force
	double mu;

	int ncollisions, ncollisions_total;
	double liquid_volume, liquid_volume_total;    // liquid volume 

	void pre_force_lbm(int, double, double *, int);

	double fA();
	double fB();
	double fC();
};
}

#endif
#endif
