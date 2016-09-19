/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "error.h"
#include "fix_setvelocity.h"
#include "memory.h"
#include "particle.h"

using namespace PDPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSetVelocity::FixSetVelocity(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 6) error->all(FLERR,"Illegal fix setvelocity command");

	vx = vy = vz = 0.0;
	if(strcmp(arg[3], "save") != 0){
		flag_x = 1;
		vx = atof(arg[3]);
	}
	else
		flag_x = 0;
	if (strcmp(arg[4], "save") != 0){
		flag_y = 1;
		vy = atof(arg[4]);
	}
	else
		flag_y = 0;
	if (strcmp(arg[5], "save") != 0){
		flag_z = 1;
		vz = atof(arg[5]);
	}
	else
		flag_z = 0;

}

/* ---------------------------------------------------------------------- */

FixSetVelocity::~FixSetVelocity()
{

}

/* ---------------------------------------------------------------------- */

int FixSetVelocity::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixSetVelocity::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixSetVelocity::post_force()
{
	double **v = particle->v;
	double *mass = particle->mass;
	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			if (flag_x == 1)
				v[i][0] = vx;
			if (flag_y == 1)
				v[i][1] = vy;
			if (flag_z == 1)
				v[i][2] = vz;
		}
	}
}

