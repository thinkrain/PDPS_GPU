/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "domain.h"
#include "error.h"
#include "fix_wall_reflect.h"
#include "memory.h"
#include "region.h"
#include "particle.h" 
#include "psmath.h"
#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;
using namespace PsMath_NS;

enum{SPECULAR, BOUNCE_BACK, BOUNCE_FORWARD, CUSTOM};
enum{BOUND, REGION};

/* ---------------------------------------------------------------------- */

FixWallReflect::FixWallReflect(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 6) error->all(FLERR,"Illegal fix wall/reflect command");

	nwalls = 0;

	nlocal_max = 0;
	x_old = NULL;
	xold_flag = 0;

	x_flag = v_flag = f_flag = 0;

	wall_rid = NULL;
	wall_flag = NULL;

	// reflection style
	// flag = 1: specular
	// flag = 2: fully inverse
	// flag = 0: none
	int iarg = 4;
	if (!strcmp(arg[3], "specular")) {
		wall_style = SPECULAR;
		x_flag = 1;
		v_flag = 1;
		f_flag = 0;
	}
	else if (!strcmp(arg[3],"bounce_back")) {
		wall_style = BOUNCE_BACK;
		x_flag = 2;
		v_flag = 2;
		f_flag = 0;
	}
	else if (!strcmp(arg[3],"bounce_forward")) {
		wall_style = BOUNCE_FORWARD;
		x_flag = 1;
		v_flag = 2;
		f_flag = 0;
	}
	else if (!strcmp(arg[3],"custom")) {
		wall_style == CUSTOM;
		while (iarg < narg) {
			if (!strcmp(arg[iarg],"pos")) {
				if (!strcmp(arg[iarg+1], "specular")) x_flag = 1;
				else if (!strcmp(arg[iarg+1], "inverse")) x_flag = 2;
				else if (!strcmp(arg[iarg+1], "none")) x_flag = 0;
				else error->all(FLERR, "Illega fix wall/reflect force option");
			}
			else if (!strcmp(arg[iarg],"vel")) {
				if (!strcmp(arg[iarg+1], "specular")) v_flag = 1;
				else if (!strcmp(arg[iarg+1], "specular")) v_flag = 2;
				else if (!strcmp(arg[iarg+1], "none")) v_flag = 0;
				else error->all(FLERR, "Illega fix wall/reflect force option");
			}
			else if (!strcmp(arg[iarg],"force")) {
				if (!strcmp(arg[iarg+1], "specular")) f_flag = 1;
				else if (!strcmp(arg[iarg+1], "inverse")) f_flag = 2;
				else if (!strcmp(arg[iarg+1], "none")) f_flag = 0;
				else error->all(FLERR, "Illega fix wall/reflect force option");
			}
			else error->all(FLERR, "Illegal fix wall/reflect custom options");
			iarg += 2;
		} 
	}  // else if (!strcmp(arg[3],"custom"))
	else error->all(FLERR, "Illegal fix wall/reflect style");

	// parse options

	int nwalls_initial = narg - 4;

	wall_rid = new int[nwalls_initial];
	wall_flag = new int[nwalls_initial];

	while (iarg < narg) {
		if (!strcmp(arg[iarg], "xlo")) {
			wall[nwalls] = 0;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "xhi")) {
			wall[nwalls] = 1;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "ylo")) {
			wall[nwalls] = 2;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "yhi")) {
			wall[nwalls] = 3;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		} 
		else if (!strcmp(arg[iarg], "zlo")) {
			if (domain->dim == 2) error->all(FLERR, "It is a 2D simulation");
			wall[nwalls] = 4;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "zhi")) {
			if (domain->dim == 2) error->all(FLERR, "It is a 2D simulation");
			wall[nwalls] = 5;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg+1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "region")) {
			wall_flag[nwalls] = REGION;
			wall_rid[nwalls] = domain->find_region(arg[iarg+1]);
			if (wall_rid[nwalls] == -1) error->all(FLERR, "Cannot find the region id");
			nwalls++;
		}
		else error->all(FLERR, "Illegal fix wall/reflect options");
		iarg += 2;
	} // while (iarg < narg)
}

/* ---------------------------------------------------------------------- */

FixWallReflect::~FixWallReflect()
{
	if (xold_flag == 1) {
		memory->destroy(x_old);
	}
}

/* ---------------------------------------------------------------------- */

void FixWallReflect::init()
{
	// fully inversed on the position needs pos of paritcle at last time step
	if (x_flag == 2) xold_flag = 1; 	
}

/* ---------------------------------------------------------------------- */

int FixWallReflect::setmask()
{
	int mask = 0;
	mask |= POST_INTEGRATE;
	mask |= PRE_INTEGRATE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallReflect::pre_integrate()
{
	// full inverse needs particle pos at last time step
	if (xold_flag == 1) {
		int nlocal = particle->nlocal;
		double **x = particle->x;

		if (nlocal > nlocal_max) {
			nlocal_max = nlocal;
			memory->destroy(x_old);
			memory->create(x_old, nlocal_max, 3, "FixWallReflect: x_old");
		}

		for (int i = 0; i < nlocal; i++)
		for (int j = 0; j < 3; j++) {
			x_old[i][j] = x[i][j];
		}
	}
}

/* ---------------------------------------------------------------------- */

void FixWallReflect::post_integrate()
{
	int dim, side;
	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	double dt = update->dt;

	int rid;
	double vel_wall[3], dist;
	Region **regions = domain->regions;

	double inside_vec[3];
	double z_vec[3];
	double n[3], pro_coords[3];
	double dot_pro;

	// loop for all walls
	for (int iwall = 0; iwall < nwalls; iwall++) {
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				if (wall_flag[iwall] == BOUND) {
					dim = wall[iwall] / 2;
					side = wall[iwall] % 2;
					// accross the lower side
					if (side == 0) {
						if (x[i][dim] < coords0[iwall]) {
							if (x_flag == 1) {                 // specular
								x[i][dim] = coords0[iwall] + (coords0[iwall] - x[i][dim]);
							}
							else if (x_flag == 2) {            // fully inversed
								double x_new = coords0[iwall] + (coords0[iwall] - x[i][dim]);
								for (int j = 0; j < domain->dim; j++) {
									if (j != dim) {
										x[i][j] += (x_old[i][j] - x[i][j]) / (x_old[i][dim] - x[i][dim])*(x_new - x[i][dim]);
									}
								}
								x[i][dim] = x_new;
							}
							// inverse velocity components
							if (v_flag == 1) {
								v[i][dim] = -v[i][dim];    // specular
							}
							else if (v_flag == 2) {
								for (int j = 0; j < domain->dim; j++) v[i][j] = -v[i][j]; // fully inversed
							}
							if (f_flag == 1) {
								f[i][dim] = -f[i][dim];    // specular
							}
							else if (f_flag == 2) {
								for (int j = 0; j < domain->dim; j++) f[i][j] = -f[i][j]; // fully inversed
							}
						} // if (x[i][dim] < coords0[iwall])
						
					} // if (side == 0) {
					// accross the upper side
					else if (side == 1) {
						if (x[i][dim] > coords0[iwall]) {
							if (x_flag == 1) {
								x[i][dim] = coords0[iwall] - (x[i][dim] - coords0[iwall]);
							}
							else if (x_flag == 2) {
								double x_new = coords0[iwall] - (x[i][dim] - coords0[iwall]);
								for (int j = 0; j < domain->dim; j++) {
									if (j != dim) {
										x[i][j] += (x_old[i][j] - x[i][j]) / (x[i][dim] - x_old[i][dim])*(x[i][dim] - x_new);
									}
								}
								x[i][dim] = x_new;
							}
							// inverse velocity components
							if (v_flag == 1) {
								v[i][dim] = -v[i][dim];    // specular
							}
							else if (v_flag == 2) {
								for (int j = 0; j < domain->dim; j++) v[i][j] = -v[i][j]; // fully inversed
							}
							// inverse force components
							if (f_flag == 1) {
								f[i][dim] = -f[i][dim];    // specular
							}
							else if (f_flag == 2) {
								for (int j = 0; j < domain->dim; j++) f[i][j] = -f[i][j]; // fully inversed
							}
						} // if (x[i][dim] > coords0[iwall])
					
					} // else if (side == 1)
				
				
				} // if (wall_flag[iwall] == BOUND) {
				else if (wall_flag[iwall] == REGION) {
					rid = wall_rid[iwall];
					for (int j = 0; j < 3; j++) {
						vel_wall[j] = 0.0;
					}
					if (regions[rid]->dynamic_flag) {
						for (int j = 0; j < 3; j++) vel_wall[j] = regions[rid]->v_coords[0][j];
					}
					// 2D is different from 3D
					if (domain->dim == 2) {
						z_vec[0] = 0.0;
						z_vec[1] = 0.0;
						z_vec[2] = 1.0;
						if (x_flag == 1) {
							Vec_Cross_Prod_3D(inside_vec, z_vec, regions[rid]->normal);
							dist = regions[rid]->find_projection_point(pro_coords, n, particle->x[i]);
							dot_pro = Vec_Dot_Prod_3D(n, inside_vec);
							if (dot_pro < 0) {
								for (int j = 0; j < 3; j++) {
									x[i][j] = x[i][j] + (-2 * dist * n[j]);
								}
							}
						}
						else if (x_flag == 2) {
							// needs further development
						}
						if (v_flag == 1) {
							// needs further development
						}
						else if (v_flag == 2) {
							for (int j = 0; j < 3; j++) v[i][j] = 2*vel_wall[j] - v[i][j]; // fully inversed
						}
						// inverse force components
						if (f_flag == 1) {
							// needs further development
						}
						else if (f_flag == 2) {
							// needs further development
						}
					}
				}
			} // if (mask[i] & groupbit)
		} // for (i = 0; i < nlocal; i++)
	} // for (int iwall = 0; iwall < nwalls; iwall++)
}
