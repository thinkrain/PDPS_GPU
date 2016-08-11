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
#include "fix_nuc.h"
#include "particle.h"
#include "particle_type.h"
#include "region.h"
#include "random_park.h"
#include "update.h"
#include "group.h"
#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;
enum{UP, DOWN, FRONT, BACK, LEFT, RIGHT};

/* ---------------------------------------------------------------------- */

FixNuc::FixNuc(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 12) {
		error->all(FLERR,"Illegal fix Nuc command");
	}
	
	region_flag = 0;

	int iarg;
	tid = atof(arg[3]);
	int ngid = group->find_group(arg[4]);
	newgid = group->bitmask[ngid];
	frequency = atoi(arg[5]);
	gengid = group->find_group(arg[6]);
	radius_bubble = atof(arg[7]);
	mass_bubble = atof(arg[8]);
	rho_bubble = atof(arg[9]);
	iarg = 10;
	if (!strcmp(arg[iarg], "up"))
		direction = UP;
	else if (!strcmp(arg[iarg], "down"))
		direction = DOWN;
	else if (!strcmp(arg[iarg], "front"))
		direction = FRONT;
	else if (!strcmp(arg[iarg], "back"))
		direction = BACK;
	else if (!strcmp(arg[iarg], "left"))
		direction = LEFT;
	else if (!strcmp(arg[iarg], "right"))
		direction = RIGHT;
	else error->all(FLERR, "Illegal command option");
	seed = atoi(arg[11]);
	count = 0;

}

/* ---------------------------------------------------------------------- */

FixNuc::~FixNuc()
{

}

/* ---------------------------------------------------------------------- */

int FixNuc::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixNuc::init()
{

}

/* ---------------------------------------------------------------------- */

void FixNuc::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixNuc::post_force()
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
	int rised = 0;
	int gennum = group->glocal[gengid];
	//	balance the speed while using more processors
	if (update->ntimestep == 0){
		int rank, nproc;
		MPI_Comm_rank(mworld, &rank);
		MPI_Comm_size(mworld, &nproc);
		count = seed / nproc * rank;
	}
	

	if (update->ntimestep % 10 == 0){

		for (int j = 0; j < nlocal; j++) {
			int i = (j + update->ntimestep * seed) % nlocal;
			if (mask[i] & groupbit) {
				if (mask[i] & newgid)
					continue;
				count++;
				//		generate a new particle, i.e. transfer a existing particle to the desired type
				if (count > seed){
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
					rised++;
					
					if (direction == UP){
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
					count = 0;

				}
				if (rised >= frequency)		//  control the maximum generation speed
					break;

			}
		}

	}

	
}
