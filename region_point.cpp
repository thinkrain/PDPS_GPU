/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "error.h"
#include "output.h"
#include "phy_const.h"
#include "region_point.h"

using namespace PDPS_NS;
using namespace PhyConst;

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

RegionPoint::RegionPoint(PDPS *ps, int narg, char **arg) : Region(ps, narg, arg)
{
	if ((narg - 2) % 4 != 0) error->all(FLERR, "Illegal polygon points");
	nps = (narg - 2) / 4;

	point_flag = 1;
	line_flag = 0;
	plane_flag = 0;
	volume_flag = 0;

	coords = new double *[nps];
	p_name = new char *[nps];
	for (int i = 0; i < nps; i++) {
		coords[i] = new double[3];
		p_name[i] = new char[128];	
	}

	int iarg = 2;

	int ip = 0;
	while (iarg < narg) {
		strcpy(p_name[ip], arg[iarg]);
		coords[ip][0] = atof(arg[iarg+1]);
		coords[ip][1] = atof(arg[iarg+2]);
		coords[ip][2] = atof(arg[iarg+3]);
		iarg += 4;
		ip++;
	}
}

/* ---------------------------------------------------------------------- */

RegionPoint::~RegionPoint()
{
	
	
}

/* ---------------------------------------------------------------------- */

int RegionPoint::inside(double *x0)
{
	int counter;

	for (int i = 0; i < nps; i++) {
		counter = 0;
		for (int j = 0; j < 3; j++) {
			if (fabs(x0[j] - coords[i][j]) < EPSILON) counter++;
		}
		if (counter == 3) return 1;
	}

	return 0;
}

/* ---------------------------------------------------------------------- */

double RegionPoint::find_distance(double *x0)
{


	return 0.0;
}

/* ---------------------------------------------------------------------- */

int RegionPoint::find_point(char * str)
{
	int pid = -1;

	for (int i = 0; i < nps; i++) {
		if (strcmp(str, p_name[i]) == 0) {
			pid = i;
			break;
		}
	}


	return pid;
}
