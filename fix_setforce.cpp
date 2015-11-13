/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "error.h"
#include "fix_setforce.h"
#include "memory.h"
#include "particle.h"

using namespace PDPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSetForce::FixSetForce(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 6) error->all(FLERR,"Illegal fix setforce command");

	fx = fy = fz = 0.0;

	fx = atof(arg[3]);
	fy = atof(arg[4]);
	fz = atof(arg[5]);
}

/* ---------------------------------------------------------------------- */

FixSetForce::~FixSetForce()
{

}

/* ---------------------------------------------------------------------- */

int FixSetForce::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixSetForce::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixSetForce::post_force()
{
	double **f = particle->f;
	double *mass = particle->mass;
	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			f[i][0] = fx;
			f[i][1] = fy;
			f[i][2] = fz;
		}
	}
}

