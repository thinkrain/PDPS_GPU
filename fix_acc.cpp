/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "domain.h"
#include "error.h"
#include "fix_acc.h"
#include "particle.h"
#include "region.h"

#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAcc::FixAcc(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 5) {
		error->all(FLERR,"Illegal fix acc command");
	}
	
	xacc = yacc = zacc = 0.0;

	region_flag = 0;

	int iarg;
	iarg = 3;
	while (iarg < narg) {
		if (!strcmp(arg[iarg],"gravity")) {
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

FixAcc::~FixAcc()
{

}

/* ---------------------------------------------------------------------- */

int FixAcc::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixAcc::init()
{

}

/* ---------------------------------------------------------------------- */

void FixAcc::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixAcc::post_force()
{
	double **x = particle->x;
	double **f = particle->f;
	double *mass = particle->mass;
	double *rmass = particle->rmass;
	int rmass_flag = particle->rmass_flag;
	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	double massone;

	int inside_flag;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			if (region_flag == 1) {
				inside_flag = domain->regions[rid]->inside(x[i]);
				if (inside_flag == 0) continue;
			}
			if (rmass_flag) massone = rmass[i];
			else massone = mass[type[i]];
			f[i][0] += massone*xacc;
			f[i][1] += massone*yacc;
			f[i][2] += massone*zacc;
		}
	}
}
