/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/reflect,FixWallReflect)

#else

#ifndef PS_FIX_WALL_REFLECT_H
#define PS_FIX_WALL_REFLECT_H

#include "fix.h"

namespace PDPS_NS {

class FixWallReflect : public Fix {
public:
	FixWallReflect(class PDPS *, int, char **);
	virtual ~FixWallReflect();
	int setmask();
	void init();
	void pre_integrate();
	void post_integrate();

protected:
	int nwalls;               // number of walls
	int wall_style;           // wall_flag = 0: BOUNCE_BACK  
                              //           = 1: SPECTRUM_REFLECT
						      //           = 2: MAXWELL
	                          //           = 3: CUSTOM
	int wall[6];              // 0: xlo; 1: xhi .... 5: zhi

	int *wall_rid;            // region id for the wall
	int *wall_flag;           // flag for wall  


	int x_flag, v_flag, f_flag;
	double coords0[6];        // store box dimension
	double x_scale[6][3], v_scale[6][3], f_scale[6][3];  // scale factor for vel and force at each wall

	double **x_old;
	int xold_flag;
	int nlocal_max;
};

}

#endif
#endif
