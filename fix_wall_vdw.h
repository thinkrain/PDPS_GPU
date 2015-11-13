/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/vdw, FixWallVDW)

#else

#ifndef PS_FIX_WALL_VDW_H
#define PS_FIX_WALL_VDW_H

#include "fix.h"

namespace PDPS_NS {

	class FixWallVDW : public Fix {
	public:
		FixWallVDW(class PDPS *, int, char **);
		virtual ~FixWallVDW();
		int setmask();
		void init();
		void setup();
		void pre_force();

	protected:
		int rot_flag;             // by default, rotation is considered
		double Ha;                // Ha: Hamaker constant; 
		double lambda, b;         // lambda: dipole interaction wavelength; b: constant
		double delta_min;         // delta_min: minimum separation distance to avoid sigularity
		double cut, rneigh;       // for particle-wall interaction cutoff and rneigh

		int nwalls;               // number of walls
		int wall[6];              // 0: xlo; 1: xhi .... 5: zhi
		int wall_style;
		double coords0[6];        // store box dimension

	private:
		int *wall_rid;            // region id for the wall
		int *wall_flag;           // flag for wall 

		void pre_force_dem_lsd(int, double, double *, int);
	};
}

#endif
#endif
