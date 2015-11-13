/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "domain.h"
#include "error.h"
#include "region_intersect.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

RegionIntersect::RegionIntersect(PDPS *ps, int narg, char **arg) : Region(ps, narg, arg)
{
	point_flag = 0;
	line_flag = 0;
	plane_flag = 0;
	volume_flag = 1;

	nregions = narg - 2;
	rid_list = new int[nregions];

	int iarg = 2;
	nregions = 0;
	int rid;
	Region **regions = domain->regions;
	while (iarg < narg) {
		rid = domain->find_region(arg[iarg]);
		if (rid == -1) error->all(FLERR, "Illegal region id");
		point_flag = regions[rid]->point_flag;
		line_flag = regions[rid]->line_flag;
		plane_flag = regions[rid]->plane_flag;
		volume_flag = regions[rid]->volume_flag;
		rid_list[nregions] = rid;
		iarg++;
		nregions++;
	}

	if (point_flag + line_flag + plane_flag + volume_flag > 1) {
		error->all(FLERR, "Regions have to be the same type");
	}
	
	for (int i = 0; i < 3; i++) {
		extent_lo[i] = regions[rid_list[0]]->extent_lo[i];
		extent_hi[i] = regions[rid_list[0]]->extent_hi[i];
	}
	for (int i = 1; i < nregions; i++) {
		for (int j = 0; j < 3; j++) {
			extent_lo[j] = MAX(extent_lo[j], regions[rid_list[i]]->extent_lo[j]);
			extent_hi[j] = MIN(extent_hi[j], regions[rid_list[i]]->extent_hi[j]);
		}
	}

	for (int i = 0; i < 3; i++) {
		extent_le[i] = extent_hi[i] - extent_lo[i];
	}

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

RegionIntersect::~RegionIntersect()
{
	delete[] rid_list;
	rid_list = NULL;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegionIntersect::inside(double *x0)
{
	int inside_flag = 1;
	int rid;
	Region **regions = domain->regions;
	for (int i = 0; i < nregions; i++) {
		rid = rid_list[i];
		inside_flag = regions[rid]->inside(x0);
		if (inside_flag == 0) break;
	}

	return inside_flag;
}

/* ---------------------------------------------------------------------- */

double RegionIntersect::find_distance(double *x0)
{
	return 0.0;
}
