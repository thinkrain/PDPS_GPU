/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "region_block.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

RegionBlock::RegionBlock(PDPS *ps, int narg, char **arg) : Region(ps, narg, arg)
{
	point_flag = 0;
	line_flag = 0;
	plane_flag = 0;
	volume_flag = 1;

	xlo = atof(arg[2]);
	xhi = atof(arg[3]);
	ylo = atof(arg[4]);
	yhi = atof(arg[5]);
	zlo = atof(arg[6]);
	zhi = atof(arg[7]);

    extent_lo[0] = xlo;
    extent_hi[0] = xhi;
    extent_lo[1] = ylo;
    extent_hi[1] = yhi;
    extent_lo[2] = zlo;
    extent_hi[2] = zhi;

	extent_le[0] = xhi - xlo;
	extent_le[1] = yhi - ylo;
	extent_le[2] = zhi - zlo;

	extent_xlo = extent_lo[0];
	extent_xhi = extent_hi[0];
	extent_xle = extent_le[0];
	extent_ylo = extent_lo[1];
	extent_yhi = extent_hi[1];
	extent_yle = extent_le[1];
	extent_zlo = extent_lo[2];
	extent_zhi = extent_hi[2];
	extent_zle = extent_le[2];

	volume = extent_xle * extent_yle * extent_zle;
}

/* ---------------------------------------------------------------------- */

RegionBlock::~RegionBlock()
{
 
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegionBlock::inside(double *x0)
{
	double x, y, z;
	x = x0[0];
	y = x0[1];
	z = x0[2];
	if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi) {
		return 1;
	}

	return 0;
}

/* ---------------------------------------------------------------------- */

double RegionBlock::find_distance(double *x0)
{
	return 0.0;
}
