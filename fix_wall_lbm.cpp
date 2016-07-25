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
#include "fix_wall_lbm.h"
#include "memory.h"
#include "parallel.h"
#include "particle.h" 
#include "phy_const.h"
#include "psmath.h"
#include "region.h"
#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;
using namespace PhyConst;
using namespace PsMath_NS;

#define EPSILON 1.0e-10

enum{ BOUND, REGION };

/* ---------------------------------------------------------------------- */

FixWallLBM::FixWallLBM(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 8) error->all(FLERR, "Illegal fix wall/lbm command");

	wall_rid = NULL;
	wall_flag = NULL;

	nwalls = 0;
	drag_flag = 0;

	phi = 0.0;

	liquid_volume = liquid_volume_total = 0.0;

	// parse commands

	if ((narg - 7) % 2 != 0) error->all(FLERR, "Illegal fix wall/lbm option");
	int nwalls_initial = (narg - 6) / 2;

	wall_rid = new int[nwalls_initial];
	wall_flag = new int[nwalls_initial];

	int iarg = 3;
	while (iarg < narg - 4) {
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
		else if (!strcmp(arg[iarg], "every")) {
			nevery = atoi(arg[iarg+1]);
		}
		else if (!strcmp(arg[iarg], "file")) {
			int n = strlen(arg[iarg+1]) + 1;
			fname = new char[n];
			strcpy(fname, arg[iarg+1]);
			if (procid == 0) {
				file = fopen(fname, "w");
				if (file == NULL) {
					char str[128];
					sprintf(str, "Cannot open file %s", arg[iarg + 1]);
					error->one(FLERR, str);
				}
				fprintf(file, "step ncollisions liquid_volume\n");
			}
		}
		else if (!strcmp(arg[iarg], "drag")) {
			drag_flag = 1;
			mu = atof(arg[iarg+1]);
		}
		else error->all(FLERR, "Illegal fix wall/lbm option");
		iarg += 2;
	} // while (iarg < narg)


	if (narg - iarg != 4) error->all(FLERR, "Illegal fix wall/lbm command");

	gamma = atof(arg[iarg++]);
	Vpcb = atof(arg[iarg++]);
	ln_Vpcb = log(Vpcb);
	Vpcb_cube_root = pow(Vpcb, 1.0 / 3.0);
	sij_min = atof(arg[iarg++]);
	double R1, R12, scut;
	R1 = atof(arg[iarg++]);
	R12 = R1;
	scut = (1 + phi / 2.0) * R12 * Vpcb_cube_root;

	cut = R1 + scut;
}

/* ---------------------------------------------------------------------- */

FixWallLBM::~FixWallLBM()
{
	delete[] wall_rid;
	wall_rid = NULL;

	delete[] wall_flag;
	wall_flag = NULL;
}

/* ---------------------------------------------------------------------- */

void FixWallLBM::init()
{

}

/* ---------------------------------------------------------------------- */

void FixWallLBM::setup()
{

}

/* ---------------------------------------------------------------------- */

int FixWallLBM::setmask()
{
	int mask = 0;
	mask |= PRE_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallLBM::pre_force()
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
	double R12, scut;

	liquid_volume = liquid_volume_total = 0.0;
	ncollisions = ncollisions_total = 0;
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
				
				R12 = radius[i];
				scut = (1 + phi / 2.0) * R12 * Vpcb_cube_root;

				// check if the particle is within the separation cut-off
				if (radius_flag == 1 && sep < scut) {
					ncollisions += 1;
					liquid_volume += Vpcb * R12 * R12 *R12;
					pre_force_lbm(i, sep, n, iwall);
				}
			} // if (mask[i] & groupbit)
		} // for (int i = 0; i < nlocal; i++)
	} // for (iwall = 0; iwall < nwalls; iwall++)

	MPI_Allreduce(&ncollisions, &ncollisions_total, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&liquid_volume, &liquid_volume_total, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (procid == 0 && file) {
		if (update->ntimestep % nevery == 0) {
			fprintf(file, "%d %d %g\n", update->ntimestep, ncollisions_total, liquid_volume_total);
			fflush(file);
		}
	}
}

/* ----------------------------------------------------------------------
   Calculate LBM force
------------------------------------------------------------------------- */

void FixWallLBM::pre_force_lbm(int iparticle, double sep, double *n, int iwall)
{
	double drijn, drijnx, drijny, drijnz;                // normal displacement
	double drijtx, drijty, drijtz;                       // tangential displacement
	double vijx, vijy, vijz;                             // relative velocity: vij = vj - vi
	double vijn, vijnx, vijny, vijnz;                    // relative velocity along normal direction
	double vijt, vijt_inv, vijtx, vijty, vijtz;          // relative velocity along tangential direction
	double omegainij[3], omegajnij[3];                   // omega_i cross nij
	double fvijnx, fvijny, fvijnz;                        // normal force
	double fvijtx, fvijty, fvijtz;                       // tangential force
	double fpair, R12, Sij;
	double nx, ny, nz;                                   // unit vector along normal direction
	double tx, ty, tz;
	double A, B, C;
	double Li[3], Lj[3];                                 // vector from the center of particle to the contact point 

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	double **omega = particle->omega;
	double *radius = particle->radius;
	Region **regions = domain->regions;

	nx = n[0];
	ny = n[1];
	nz = n[2];

	R12 = radius[iparticle];            // harmonice mean radius:  1 / R12 = 2/(1/R1 + 1/R2) 
	                                    // assume cylinder's wall's curvature is inifinity
	                                    // This may not be accurate, needs further attention
	drijn = -sep;
	if (sep < sij_min) sep = sij_min;
	Sij = sep / R12;

	// may need to consider wall translational velocity in the future
	double vel_wall[3];
	for (int i = 0; i < 3; i++) vel_wall[i] = 0.0;

	// needs to be ruther investigated, because this drag force may not be suitable for particle-wall interaction
	if (drag_flag) {
		if (wall_flag[iwall] == REGION) {
			int rid = wall_rid[iwall];
			if (regions[rid]->dynamic_flag) {
				for (int j = 0; j < 3; j++) vel_wall[j] = regions[rid]->v_coords[0][j];
			}
		}
		// vij = vi - vj + (wi cross Ri) - (wj cross Rj) (relative velocity)
		vijx = v[iparticle][0] - vel_wall[0];
		vijy = v[iparticle][1] - vel_wall[1];
		vijz = v[iparticle][2] - vel_wall[2];

		if (particle->radius_flag) {
			for (int i = 0; i < 3; i++) Li[i] = (radius[iparticle] - drijn)*(-n[i]);
			Vec_Cross_Prod_3D(omegainij, omega[iparticle], Li);
			vijx += omegainij[0];
			vijy += omegainij[1];
			vijz += omegainij[2];
			if (wall_flag[iwall] == REGION) {
				int rid = wall_rid[iwall];
				if (regions[rid]->rotate_flag) {
					double pro_coord[3], nn[3];
					regions[rid]->find_projection_point(pro_coord, nn, x[iparticle]);
					regions[rid]->rot_axis->find_vector(Lj, pro_coord);
					Vec_Cross_Prod_3D(omegajnij, regions[rid]->omega, Lj);
					// the direction of the center of the region to the contact point is opposite of n
					vijx += -omegajnij[0];
					vijy += -omegajnij[1];
					vijz += -omegajnij[2];
				}
			}
		}

		// |vijn| = vij . nij
		vijn = vijx * nx + vijy * ny + vijz * nz;
		// vijn = |vijn| . nij
		vijnx = vijn*nx;
		vijny = vijn*ny;
		vijnz = vijn*nz;

		// vijt = vij - (vij . nij) . nij
		vijtx = vijx - vijnx;
		vijty = vijy - vijny;
		vijtz = vijz - vijnz;
		vijt = sqrt(vijtx*vijtx + vijty*vijty + vijtz*vijtz);
		if (fabs(vijt) < EPSILON) vijt_inv = 0.0;
		else vijt_inv = 1.0 / vijt;

		// tij = tij/|tij|
		tx = vijtx * vijt_inv;
		ty = vijty * vijt_inv;
		tz = vijtz * vijt_inv;
		
		fvijnx = -6 * PI * mu * vijnx * R12 * R12 / sep;
		fvijny = -6 * PI * mu * vijny * R12 * R12 / sep;
		fvijnz = -6 * PI * mu * vijnz * R12 * R12 / sep;

		fvijtx = -6 * PI * mu * R12 * vijtx * (8.0 / 15 * log(R12 / sep) + 0.9588);
		fvijty = -6 * PI * mu * R12 * vijty * (8.0 / 15 * log(R12 / sep) + 0.9588);
		fvijtz = -6 * PI * mu * R12 * vijtz * (8.0 / 15 * log(R12 / sep) + 0.9588);
	}

	A = fA();
	B = fB();
	C = fC();

	fpair = -PI * R12 * gamma * (exp(A * Sij + B) + C);

	f[iparticle][0] += fpair * nx;
	f[iparticle][1] += fpair * ny;
	f[iparticle][2] += fpair * nz;

	if (drag_flag) {
		f[iparticle][0] += fvijnx + fvijtx;
		f[iparticle][1] += fvijny + fvijty;
		f[iparticle][2] += fvijnz + fvijtz;
	}
}

/* ----------------------------------------------------------------------
   Check the paper for details:
   "2010 Mixing characteristics of wet granular matter in a bladed mixer"
   "1998 Numerical simulation of cohesive powder behavior in a fluidized bed"
------------------------------------------------------------------------- */

double FixWallLBM::fA()
{
	double res;

	res = -1.9 * pow(Vpcb, -0.51);

	return res;
}

/* ---------------------------------------------------------------------- */

double FixWallLBM::fB()
{
	double res;

	res = (-0.016*ln_Vpcb - 0.76)*phi*phi - 0.12*ln_Vpcb + 1.2;

	return res;
}

/* ---------------------------------------------------------------------- */

double FixWallLBM::fC()
{
	double res;

	res = 0.013 * ln_Vpcb + 0.18;

	return res;
}
