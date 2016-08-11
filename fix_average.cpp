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
#include "fix_average.h"
#include "particle.h"
#include "region.h"
#include "group.h"
#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAverage::FixAverage(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 5) {
		error->all(FLERR,"Illegal fix average command");
	}
	
	
	int iarg;
	iarg = 3;
	force_flag = 0;
	force_x = 0;
	force_y = 0;
	force_z = 0;
	gid = group->find_group(arg[1]);
	while (iarg < narg) {
		if (!strcmp(arg[iarg],"force")) {
			iarg += 1;
			force_flag = 1;
			if (!strcmp(arg[iarg], "x")){
				iarg += 1;
				force_x = 1;
			}
			if (!strcmp(arg[iarg], "y")){
				iarg += 1;
				force_y = 1;
			}
			if (!strcmp(arg[iarg], "z")){
				iarg += 1;
				force_z = 1;
			}
		}
		else error->all(FLERR, "Illegal command option");
	}

}

/* ---------------------------------------------------------------------- */

FixAverage::~FixAverage()
{

}

/* ---------------------------------------------------------------------- */

int FixAverage::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixAverage::init()
{

}

/* ---------------------------------------------------------------------- */

void FixAverage::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixAverage::post_force()
{
	double **x = particle->x;
	double **f = particle->f;
	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;

	if (force_flag == 1){
		double xforce_sum, yforce_sum, zforce_sum, xforce_ave, yforce_ave, zforce_ave;
		xforce_sum = yforce_sum = zforce_sum = 0.0;
		int n_ave;
		if (force_x){
			for (int i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
					xforce_sum += f[i][0];
				}
			}
			MPI_Allreduce(&xforce_sum, &xforce_ave, 1, MPI_DOUBLE, MPI_SUM, mworld);
			xforce_ave /= group->gparticles[gid];
			for (int i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
						f[i][0] = xforce_ave;
				}
			}
		}
		
		if (force_y){
			for (int i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
					yforce_sum += f[i][1];
				}
			}
			MPI_Allreduce(&yforce_sum, &yforce_ave, 1, MPI_DOUBLE, MPI_SUM, mworld);
			yforce_ave /= group->gparticles[gid];
			for (int i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
					f[i][1] = yforce_ave;
				}
			}
		}

		if (force_z){
			for (int i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
					zforce_sum += f[i][2];
				}
			}
			MPI_Allreduce(&zforce_sum, &zforce_ave, 1, MPI_DOUBLE, MPI_SUM, mworld);
			zforce_ave /= group->gparticles[gid];
			for (int i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
					f[i][2] = zforce_ave;
				}
			}
		}



	}


}
