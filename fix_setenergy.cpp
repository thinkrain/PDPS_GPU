/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "error.h"
#include "fix_setenergy.h"
#include "memory.h"
#include "particle.h"

using namespace PDPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSetEnergy::FixSetEnergy(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 4) error->all(FLERR,"Illegal fix setenergy command");
	energy = atof(arg[3]);
	
}

/* ---------------------------------------------------------------------- */

FixSetEnergy::~FixSetEnergy()
{

}

/* ---------------------------------------------------------------------- */

int FixSetEnergy::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixSetEnergy::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixSetEnergy::post_force()
{

	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	double *e = particle->e;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			e[i] = energy;
		}
	}
}

