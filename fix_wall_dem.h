/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/dem,FixWallDEM)

#else

#ifndef PS_FIX_WALL_DEM_H
#define PS_FIX_WALL_DEM_H

#include "fix.h"

namespace PDPS_NS {

class FixWallDEM : public Fix {
public:
	FixWallDEM(class PDPS *, int, char **);
	virtual ~FixWallDEM();
	int setmask();
	void init();
	void setup();
	void pre_force();

protected:
	int rot_flag;             // by default, rotation is considered
	double kn, Cn;            // normal force
	double kt, Ct;            // tangential force
	double mu;                // friction coefficient
	double e;                 // restitution coefficient
	double cut, rneigh;       // for particle-wall interaction cutoff and rneigh

	int nwalls;               // number of walls
	int wall[6];              // 0: xlo; 1: xhi .... 5: zhi
	int wall_style;
	double coords0[6];        // store box dimension


private:
	int *wall_rid;            // region id for the wall
	int *wall_flag;           // flag for wall 

	//double **vel_wall;        // velocity of the wall
	//double **omega;           // rotation speed of the wall

	class PairList *pair_list;
	int tbsize;

	void set_pairlist();
	void tag2str(char *, int, int);

	void pre_force_dem_lsd(int, int, double, double *, int);
};

}

#endif
#endif
