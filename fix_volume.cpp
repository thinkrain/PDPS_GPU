/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "error.h"
#include "fix_volume.h"
#include "force.h"
#include "memory.h"
#include "pair.h"
#include "particle.h"
#include "domain.h"
#include "phy_const.h"

using namespace PDPS_NS;
using namespace FixConst;
using namespace PhyConst;

enum{DRAG_STOKES,BUOYANCY,CUSTOM};

/* ---------------------------------------------------------------------- */

FixVolume::FixVolume(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 5) error->all(FLERR,"Illegal fix setforce command");
	if (particle->sphere_flag == 0) error->all(FLERR, "Illegal particle type");
	xacc = yacc = zacc = 0.0;
	gravity_flag = 0;
	int iarg;
	iarg = 3;
	while (iarg < narg) {
		if (!strcmp(arg[iarg],"gravity")) {
			gravity_flag = 1;
			if (!strcmp(arg[iarg+1],"x")) xacc = atof(arg[iarg+2]);
			if (!strcmp(arg[iarg+1],"y")) yacc = atof(arg[iarg+2]);
			if (!strcmp(arg[iarg+1],"z")) zacc = atof(arg[iarg+2]);
			iarg += 3;
		}
		else if (!strcmp(arg[iarg], "region")) {
			rid = domain->find_region(arg[iarg+1]);
			if (rid == -1) error->all(FLERR, "Cannot find the region id");
			region_flag = 1;
			iarg += 2;
		}
		else error->all(FLERR, "Illegal command option");
	} 
}

/* ---------------------------------------------------------------------- */

FixVolume::~FixVolume()
{
	
}

/* ---------------------------------------------------------------------- */

int FixVolume::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixVolume::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixVolume::post_force()
{
	double **x = particle->x;
	double *density = particle->density;
	double *radius = particle->radius;
	double xlo,xhi,ylo,yhi,zlo,zhi;
	int *type = particle->type;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	double expansion, pressure;
	double xdepth, ydepth, zdepth;
	xlo = domain->boxlo[0];									// by now only work for box domain
	xhi = domain->boxhi[0];
	ylo = domain->boxlo[1];
	yhi = domain->boxhi[1];
	zlo = domain->boxlo[2];
	zhi = domain->boxhi[2];
	
	if (gravity_flag == 1){
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				xdepth = xhi - x[i][0];
				ydepth = yhi - x[i][1];
				zdepth = zhi - x[i][2];
				pressure = sqrt(xdepth * xacc + xdepth * xacc + ydepth * yacc + ydepth * yacc + zdepth * zacc + zdepth * zacc); 
				expansion = 1.0 / pressure;
				density[i] = density[i] / expansion;
				radius[i] = radius[i] * exp(espansion, 3);
			}
		}
	}
}

/* ---------------------------------------------------------------------- */

