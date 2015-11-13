/* ----------------------------------------------------------------------
< Particle Dynamics Parallel Simulator (PDPS) >
Copyright(C) <2014>  <Author: Lingqi Yang>
Email: ly2282@columbia.edu

This program is free software : you can redistribute it and / or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "domain.h"
#include "error.h"
#include "force.h"
#include "fix_wall_vdw.h"
#include "memory.h"
#include "parallel.h"
#include "particle.h" 
#include "region.h"

using namespace PDPS_NS;
using namespace FixConst;
#define EPSILON 1.0e-10

enum{BOUND, REGION};

/* ---------------------------------------------------------------------- */

FixWallVDW::FixWallVDW(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 5) error->all(FLERR, "Illegal fix wall/vdw command");

	wall_rid = NULL;
	wall_flag = NULL;

	nwalls = 0;

	// parse options
	int iarg = 3;

	if ((narg - 8) % 2 != 0) error->all(FLERR, "Illegal fix wall/vdw option");
	int nwalls_initial = (narg - 8) / 2;

	wall_rid = new int[nwalls_initial];
	wall_flag = new int[nwalls_initial];

	while (iarg < narg - 5) {
		if (!strcmp(arg[iarg], "xlo")) {
			wall[nwalls] = 0;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg + 1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "xhi")) {
			wall[nwalls] = 1;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg + 1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "ylo")) {
			wall[nwalls] = 2;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg + 1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "yhi")) {
			wall[nwalls] = 3;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg + 1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "zlo")) {
			if (domain->dim == 2) error->all(FLERR, "It is a 2D simulation");
			wall[nwalls] = 4;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg + 1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "zhi")) {
			if (domain->dim == 2) error->all(FLERR, "It is a 2D simulation");
			wall[nwalls] = 5;
			wall_flag[nwalls] = BOUND;
			coords0[nwalls] = atof(arg[iarg + 1]);
			nwalls++;
		}
		else if (!strcmp(arg[iarg], "region")) {
			wall_flag[nwalls] = REGION;
			wall_rid[nwalls] = domain->find_region(arg[iarg + 1]);
			if (wall_rid[nwalls] == -1) error->all(FLERR, "Cannot find the region id");
			nwalls++;
		}
		else error->all(FLERR, "Illegal fix wall/vdw option");

		iarg += 2;
	} // while (iarg < narg)


	if (narg - iarg != 5) error->all(FLERR, "Illegal fix wall/vdw command");

	rot_flag = 1;
	b = atof(arg[iarg++]);
	lambda = atof(arg[iarg++]);
	delta_min = atof(arg[iarg++]);
	Ha = atof(arg[iarg++]);
	cut = atof(arg[iarg++]);
}

/* ---------------------------------------------------------------------- */

FixWallVDW::~FixWallVDW()
{

}

/* ---------------------------------------------------------------------- */

void FixWallVDW::init()
{

}

/* ---------------------------------------------------------------------- */

void FixWallVDW::setup()
{

}

/* ---------------------------------------------------------------------- */

int FixWallVDW::setmask()
{
	int mask = 0;
	mask |= PRE_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallVDW::pre_force()
{
	int i, j, iwall, itag;
	int table, index;

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	int *mask = particle->mask;
	int *tag = particle->tag;
	int radius_flag = particle->radius_flag;
	double *radius = particle->radius;
	int nlocal = particle->nlocal;

	int rid, inside_flag;
	int dim, side;
	double dist, sep, n[3];

	for (iwall = 0; iwall < nwalls; iwall++) {
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				if (wall_flag[iwall] == BOUND) {
					dim = wall[iwall] / 2;
					side = wall[iwall] % 2;
					for (int j = 0; j < 3; j++) n[j] = 0.0;
					if (side == 0) {
						dist = x[i][dim] - coords0[iwall];
						n[dim] = 1.0;
					}
					else {
						dist = coords0[iwall] - x[i][dim];
						n[dim] = -1.0;
					}
				}
				else if (wall_flag[iwall] == REGION) {
					rid = wall_rid[iwall];
					dist = domain->regions[rid]->find_interaction_distance(n, x[i]);
				}

				if (radius_flag) sep = fabs(dist) - radius[i];
				else sep = fabs(dist) - cut;

				// check if the particle is within the interaction range
				if (radius_flag == 0 && fabs(sep) < 0.25*cut) {
					pre_force_dem_lsd(i, sep, n, iwall);
				}
				else if (radius_flag == 1 && fabs(sep) < 0.5*radius[i]) {
					pre_force_dem_lsd(i, sep, n, iwall);
				}
			} // if (mask[i] & groupbit)
		} // for (int i = 0; i < nlocal; i++)
	} // for (iwall = 0; iwall < nwalls; iwall++)
}

/* ----------------------------------------------------------------------
dem: force
------------------------------------------------------------------------- */

void FixWallVDW::pre_force_dem_lsd(int iparticle, double sij, double *n, int iwall)
{
	double fpair, Reff, sijsq;
	double nx, ny, nz;                             // unit vector along normal direction

	double **x = particle->x;
	double **f = particle->f;
	double *radius = particle->radius;

	nx = n[0];
	ny = n[1];
	nz = n[2];

	if (particle->radius_flag) Reff = radius[iparticle];
	else Reff = cut / 2.0;

	if (sij < delta_min) {
		sij = delta_min;
	}
	sijsq = sij * sij;
	
	fpair = -Ha * Reff / (6 * sijsq) * (1 - 1.0 / (1.0 + lambda / (b*sij)));

	f[iparticle][0] += fpair * nx;
	f[iparticle][1] += fpair * ny;
	f[iparticle][2] += fpair * nz;
}
