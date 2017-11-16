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
#include "parallel.h"

using namespace PDPS_NS;
using namespace FixConst;

#define PI 3.1416

enum{ RTI, HV, CONSTANT };

/* ---------------------------------------------------------------------- */

FixMove::FixMove(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 5) {
		error->all(FLERR, "Illegal fix move command");
	}


	int iarg;
	iarg = 3;
	gid = group->find_group(arg[1]);
	while (iarg < narg) {
		if (!strcmp(arg[iarg], "RTI")) {
			move_style = RTI;
			kx = atof(arg[4]);
			ky = atof(arg[5]);
			iarg += 3;
		}
		else if (!strcmp(arg[iarg], "HV")) {
			move_style = HV;
			axis_h = atof(arg[4]);
			axis1_r = atof(arg[5]);
			axis1_w = atof(arg[6]);
			axis2_r = atof(arg[7]);
			axis2_w = atof(arg[8]);
			thita = atof(arg[9]);
			width = atof(arg[10]);
			startstep = atoi(arg[11]);
			processstep = atof(arg[12]);
			endstep = atoi(arg[13]);
			ang0_1 = atof(arg[14]);
			ang0_2 = atof(arg[15]);
			cen1_x = atof(arg[16]);
			cen1_y = atof(arg[17]);
			direction = atoi(arg[18]);
			iarg += 16;
		}
		else if (!strcmp(arg[iarg], "CONSTANT")){
			move_style = CONSTANT;
			startstep = atoi(arg[4]);
			endstep = atoi(arg[5]);
			move_v = atof(arg[6]);
			move_dx = atof(arg[7]);
			move_dy = atof(arg[8]);
			move_dz = atof(arg[9]);
			iarg += 7;
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
	double *hlocal = particle->hlocal;		//	use this variable as the distance towards axis
	double *volume = particle->volume;	// use this variable as the angle towards axis
	double dz;
	double ang1, ang2, axis1_x, axis1_y, axis2_x, axis2_y;
	if (move_style == RTI){
		if (once_flag == 0){
			for (int i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
					dz = radius[i] / 2.0 * (1 - cos(kx * x[i][0])) * (1 - cos(ky * x[i][1]));
					x[i][2] -= dz;
				}
			}
			once_flag = 1;

		}
	}
	else if (move_style == HV){
		double time = (update->ntimestep - startstep) * update->dt;
		if (time >= 0 && update->ntimestep < endstep){
			double process = (update->ntimestep - startstep) * 1.0 / processstep;
			double processtime = processstep * update->dt;
			double endprocess = (update->ntimestep - endstep + processstep) * 1.0 / processstep;
			if (process < 1){
				ang1 = ang0_1 + axis1_w * process / 2.0 * time;
				ang2 = ang0_2 + axis2_w * process / 2.0 * time;
			}
			else if (endprocess < 0){
				ang1 = ang0_1 + axis1_w * (time - processtime / 2.0);
				ang2 = ang0_2 + axis2_w * (time - processtime / 2.0);
			}
			else{
				ang1 = ang0_1 + axis1_w * (time - processtime / 2.0 - endprocess * endprocess * processtime / 2.0);
				ang2 = ang0_2 + axis2_w * (time - processtime / 2.0 - endprocess * endprocess * processtime / 2.0);
			}
			axis1_x = cen1_x + axis1_r * cos(ang1);
			axis1_y = cen1_y + axis1_r * sin(ang1);
			axis2_x = axis1_x + axis2_r * cos(ang2);
			axis2_y = axis1_y + axis2_r * sin(ang2);
			if (once_flag == 0){
				for (int i = 0; i < nlocal; i++){
					if (mask[i] & groupbit){
						hlocal[i] = (x[i][0] - axis1_x) * (x[i][0] - axis1_x) + (x[i][1] - axis1_y) * (x[i][1] - axis1_y);
						hlocal[i] = sqrt(hlocal[i]);
						if (hlocal[i] < 0.0001)
							volume[i] = 0;
						else{
							if (x[i][0] > axis1_x){
								if (x[i][1] > axis1_y)
									volume[i] = asin((x[i][1] - axis1_y) / hlocal[i]);
								else
									volume[i] = asin((x[i][1] - axis1_y) / hlocal[i]) + 2 * PI;
							}
							else
								volume[i] = PI - asin((x[i][1] - axis1_y) / hlocal[i]);
						}
					}
				}
				once_flag = 1;
			}
			for (int i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
					if (direction == 1){
						x[i][0] = axis1_x + hlocal[i] * cos(ang1 + ang2 + volume[i]);
						x[i][1] = axis1_y + hlocal[i] * sin(ang1 + ang2 + volume[i]);
						//if (particle->tag[i] == 10){
						//	printf("T=%d procid = %d", update->ntimestep,parallel->procid);
						//	printf("axis1_x = %f\n hlocal[i]=%f\n ang1=%f\n ang2=%f\n volume[i]=%f\n", axis1_x, hlocal[i], ang1, ang2, volume[i]);
						//	printf("x[%d][0] = %f x[%d][1] = %f\n", i, x[i][0], i, x[i][1]);
						//}
					}
					else{
						x[i][0] = axis1_x - hlocal[i] * cos(ang1 + ang2 + volume[i]);
						x[i][1] = axis1_y - hlocal[i] * sin(ang1 + ang2 + volume[i]);
					}

				}
			}

		}	// time >= 0

	}	//	move_style == HV
	else if (move_style == CONSTANT){
		if (update->ntimestep >= startstep && update->ntimestep <= endstep){
			double dx, dy, dz;
			dx = move_v * move_dx * update->dt;
			dy = move_v * move_dy * update->dt;
			dz = move_v * move_dz * update->dt;
			for (int i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
					x[i][0] += dx;
					x[i][1] += dy;
					x[i][2] += dz;
				}
			}
		}	// ntimestep 

	}	// move_style == CONSTANT

}
