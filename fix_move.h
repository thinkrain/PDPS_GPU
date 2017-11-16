/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator/*

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(move, FixMove)

#else

#ifndef PS_FIX_MOVE_H
#define PS_FIX_MOVE_H

#include "fix.h"

namespace PDPS_NS {

	class FixMove : public Fix {

	public:
		FixMove(class PDPS *, int, char **);
		~FixMove();
		int setmask();
		void init();
		void setup();
		virtual void post_force();

	protected:
		int gid;                  // group id
		int move_style, once_flag;
		double kx, ky;
		int startstep, processstep, endstep, direction;
		double axis_h, axis1_r, axis1_w, axis2_r, axis2_w, thita, width, ang0_1, ang0_2, cen1_x, cen1_y;	// HV
		double move_v, move_dx, move_dy, move_dz;	// CONSTANT
	};
}

#endif
#endif
