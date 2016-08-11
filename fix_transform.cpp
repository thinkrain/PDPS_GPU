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
#include "create_particle.h"
#include "domain.h"
#include "error.h"
#include "fix_transform.h"
#include "particle.h"
#include "particle_type.h"
#include "region.h"
#include "random_mars.h"
#include "update.h"
#include "group.h"
#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;
enum{UP, DOWN, FRONT, BACK, LEFT, RIGHT};

/* ---------------------------------------------------------------------- */

FixTransform::FixTransform(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 10) {
		error->all(FLERR,"Illegal fix Transform command");
	}
	

	int iarg;
	tid = atof(arg[3]);
	int ngid = group->find_group(arg[4]);
	newgid = group->bitmask[ngid];
	ratio = atof(arg[5]);
	radius_bubble = atof(arg[6]);
	mass_bubble = atof(arg[7]);
	rho_bubble = atof(arg[8]);
	seed = atoi(arg[9]);
	//iarg = 10;
	random = NULL;
	random = new RanMars(ps, seed);
	//else error->all(FLERR, "Illegal command option");
	fixed_flag = 0;


}

/* ---------------------------------------------------------------------- */

FixTransform::~FixTransform()
{

}

/* ---------------------------------------------------------------------- */

int FixTransform::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixTransform::init()
{

}

/* ---------------------------------------------------------------------- */

void FixTransform::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixTransform::post_force()
{
	double **x = particle->x;
	double **f = particle->f;
	double *mass = particle->mass;
	double *rmass = particle->rmass;
	double *radius = particle->radius;
	double *rho = particle->rho;
	double *density = particle->density;
	double *volume = particle->volume;
	int rmass_flag = particle->rmass_flag;
	int *mask = particle->mask;
	int *bitmask = group->bitmask;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	int *tag = particle->tag;
	double temp;

	int inside_flag;
	if (fixed_flag == 0){
		//double randnum;
		int count = 0;
		int count1, count2;
		int q1, q2;
		q1 = int(ratio / (ratio + 1) * 10);
		q2 = int(ratio / (ratio + 1) * 100) - q1 * 10;
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				//randnum = random->uniform();
				//		generate a new particle, i.e. transfer a existing particle to the desired type
				count++;
				count1 = count % 10;
				count2 = (count / 10) % 10;
				if (count1 < q1 || count1 == q1 && count2 < q2){
					//	double coord[3];
					//	coord[0] = x[i][0];
					//	coord[1] = x[i][1];
					//	coord[2] = x[i][2];
					//	int k = i % gennum;
					//	int cou = 0;
					//	int l = 0;
					//	for (l = 0; cou < k; l++){
					//		if (mask[l] & gengid)
					//			cou++;
					//	}
					//	x[i][0] = x[l - 1][0];
					//	x[i][1] = x[l - 1][1];
					//	x[i][2] = x[l - 1][2];

					mask[i] |= newgid;
					type[i] = tid;
					radius[i] = radius_bubble;
					rmass[i] = mass_bubble;
					rho[i] = rho_bubble;
					density[i] = rho_bubble;
					volume[i] = 4.0 / 3.0 * 3.1416 * radius[i] * radius[i] * radius[i];
					group->glocal[newgid] = group->glocal[newgid] + 1;
					group->gparticles[newgid] = group->gparticles[newgid] + 1;
					group->glocal[groupbit] = group->glocal[groupbit] - 1;
					group->gparticles[groupbit] = group->gparticles[groupbit] - 1;


				}

			}
		}
		fixed_flag = 1;
	}

	
	


	
}
