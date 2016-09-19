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
#include "fix_move.h"
#include "particle.h"
#include "region.h"
#include "group.h"
#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;

enum{ RTI };

/* ---------------------------------------------------------------------- */

FixMove::FixMove(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 5) {
		error->all(FLERR,"Illegal fix move command");
	}
	
	
	int iarg;
	iarg = 3;
	gid = group->find_group(arg[1]);
	while (iarg < narg) {
		if (!strcmp(arg[iarg],"RTI")) {
			move_style = RTI;
			kx = atof(arg[4]);
			ky = atof(arg[5]);
			iarg += 3;
		}
		else error->all(FLERR, "Illegal command option");
	}
	once_flag = 0;
}

/* ---------------------------------------------------------------------- */

FixMove::~FixMove()
{

}

/* ---------------------------------------------------------------------- */

int FixMove::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixMove::init()
{

}

/* ---------------------------------------------------------------------- */

void FixMove::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixMove::post_force()
{
	double **x = particle->x;
	double **f = particle->f;
	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	double *radius = particle->radius;
	double dz;
	if (once_flag == 0){
		if (move_style == RTI){
			for (int i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
					dz = radius[i] / 2.0 * (1 - cos(kx * x[i][0])) * (1 - cos(ky * x[i][1]));
					x[i][2] -= dz;
				}
			}


		}
		once_flag = 1;
	}



}
    