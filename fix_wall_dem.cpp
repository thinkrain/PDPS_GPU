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
#include "fix_wall_dem.h"
#include "memory.h"
#include "neighbor.h"
#include "pair.h"
#include "pair_dem_lsd.h"
#include "parallel.h"
#include "particle.h" 
#include "pair_list.h"
#include "psmath.h"
#include "region.h"
#include "update.h"

#include "output.h"

using namespace PDPS_NS;
using namespace FixConst;
using namespace PsMath_NS;
#define EPSILON 1.0e-10

enum{LSD};
enum{BOUND, REGION};

/* ---------------------------------------------------------------------- */

FixWallDEM::FixWallDEM(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 12) error->all(FLERR,"Illegal fix wall/dem command");

	nwalls = 0;

	wall_rid = NULL;
	wall_flag = NULL;


	if (!strcmp(arg[3], "lsd")) {
		wall_style = LSD;
	}
	else error->all(FLERR, "Illegal fix wall/dem style");

	// parse options
	int iarg = 4;
	
	if ((narg-7-4) % 2 != 0) error->all(FLERR, "Illegal fix wall/dem option");
	int nwalls_initial = (narg-7-4) / 2;

	wall_rid = new int[nwalls_initial];
	wall_flag = new int[nwalls_initial];
	//vel_wall = new double*[nwalls_initial];
	/*for (int i = 0; i < nwalls_initial; i++) {
		vel_wall[i] = new double[3];
		for (int j = 0; j < 3; j++) vel_wall[i][j] = 0.0;
	}*/
	
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
		else break;
		iarg += 2;
	} // while (iarg < narg)


	if (narg - iarg != 7) error->all(FLERR, "Illegal fix wall/dem command");

	rot_flag = 1;
	e = atof(arg[iarg++]);
	kn = atof(arg[iarg++]);
	Cn = atof(arg[iarg++]);
	kt = atof(arg[iarg++]);
	Ct = atof(arg[iarg++]);
	mu = atof(arg[iarg++]);
	cut = atof(arg[iarg++]);

	rneigh = cut + neighbor->rskin;

	pair_list = NULL;
	tbsize = 10000;
	pair_list = new PairList(ps, tbsize);
	pair_list->ifields = 1;
	pair_list->dfields = 3;
}

/* ---------------------------------------------------------------------- */

FixWallDEM::~FixWallDEM()
{
	delete[] wall_rid;
	wall_rid = NULL;
	
	delete[] wall_flag;
	wall_flag = NULL;

	delete pair_list;
	pair_list = NULL;
}

/* ---------------------------------------------------------------------- */

void FixWallDEM::init()
{
	int size = particle->nlocal;
	pair_list->init_hash(size);
}

/* ---------------------------------------------------------------------- */

void FixWallDEM::setup()
{
	set_pairlist();
}

/* ---------------------------------------------------------------------- */

int FixWallDEM::setmask()
{
	int mask = 0;
	mask |= PRE_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallDEM::pre_force()
{
	int iwall, itag;

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	int *mask = particle->mask;
	int *tag = particle->tag;
	int radius_flag = particle->radius_flag;
	double *radius = particle->radius;
	int nlocal = particle->nlocal;
	double dt = update->dt;

	char str[BYTES];
	int ipair;

	int nflag = neighbor->nflag;

	if (nflag == 1) {
		set_pairlist();
	}

	PairList::HashPair **hash_pair = pair_list->hash_pair;

	int rid;
	int dim, side;
	double dist, n[3];

	for (iwall = 0; iwall < nwalls; iwall++) {
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				itag = tag[i];
				tag2str(str, itag, iwall);
				ipair = pair_list->find_hash(str);
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

				// check if the particle is within the interaction range
				if (radius_flag == 0 && fabs(dist) < cut) {
					if (ipair == -1) ipair = pair_list->insert_hash(str);
					pre_force_dem_lsd(ipair, i, cut-fabs(dist), n, iwall);
				}
				else if (radius_flag == 1 && fabs(dist) < radius[i]) {
					if (ipair == -1) ipair = pair_list->insert_hash(str);
					pre_force_dem_lsd(ipair, i, radius[i]-fabs(dist), n, iwall);
				}
				else {
					if (ipair > -1) {
						pair_list->set_zero(ipair);
					}
				}
			} // if (mask[i] & groupbit)
		} // for (int i = 0; i < nlocal; i++)
	} // for (iwall = 0; iwall < nwalls; iwall++)
}

/* ----------------------------------------------------------------------
   only reset pairlist when nflag = 1
------------------------------------------------------------------------- */

void FixWallDEM::set_pairlist()
{
	int ipair;
	int itag;
	int dim, side;
	int iwall;
	char str[BYTES];

	int nlocal = particle->nlocal;
	int radius_flag = particle->radius_flag;
	int *tag = particle->tag;
	int *mask = particle->mask;
	double **x = particle->x;
	double *radius = particle->radius;

	// initialize # of undetermined pairs to send and recv
	pair_list->set_hash();

	double dist;
	int rid;
	double n[3];

	for (iwall = 0; iwall < nwalls; iwall++) {
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				itag = tag[i];
				tag2str(str, itag, iwall);
				ipair = pair_list->find_hash(str);
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

				// check if the particle is within the interaction range
				if (radius_flag == 0 && fabs(dist) < cut) {
					if (pair_list->nbuilds_total == 1) ipair = pair_list->insert_hash(str);
					else ipair = pair_list->set_pair(str);
				}
				else if (radius_flag == 1 && fabs(dist) < radius[i]) {
					if (pair_list->nbuilds_total == 1 && ipair == -1) ipair = pair_list->insert_hash(str);
					else ipair = pair_list->set_pair(str);
				}
				else {
					if (ipair > -1) {
						pair_list->set_zero(ipair);
					}
				}
			} // if (mask[i] & groupbit)
		} // for (int i = 0; i < nlocal; i++)
	} // for (iwall = 0; iwall < nwalls; iwall++)

	if (pair_list->nbuilds_total > 1 && parallel->nprocs > 1) pair_list->exchange_pair();
}

/* ----------------------------------------------------------------------
   str = itag-fix%name-Wall%id-Side%
------------------------------------------------------------------------- */

void FixWallDEM::tag2str(char *str, int itag, int iwall)
{
	sprintf(str, "%d-fix%s-Wall%d", itag, name, iwall);
}

/* ----------------------------------------------------------------------
   dem: force
------------------------------------------------------------------------- */

void FixWallDEM::pre_force_dem_lsd(int ipair, int iparticle, double delta, double *n, int iwall)
{
	int table, index;

	double drijn, drijnx, drijny, drijnz;                // normal displacement
	double drijtx, drijty, drijtz;                       // tangential displacement
	double vijx, vijy, vijz;                             // relative velocity: vij = vj - vi
	double vijn, vijnx, vijny, vijnz;                    // relative velocity along normal direction
	double vijt, vijt_inv, vijtx, vijty, vijtz;          // relative velocity along tangential direction
	double fijn, fijnx, fijny, fijnz;                    // normal force
	double fijt_a[3], fijt, fijtx, fijty, fijtz;         // tangential force	
	double nx, ny, nz;                                   // unit vector along normal direction
	double tx, ty, tz;                                   // unit vector along tangential direction
	double omegainij[3], omegajnij[3];                   // omega_i cross nij
	double torqueij[3];                                  // torque
	double Li[3], Lj[3];                                 // vector from the center of particle to the contact point 

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	double *radius = particle->radius;
	double **omega = particle->omega;
	double **torque = particle->torque;

	double dt = update->dt;

	Region **regions = domain->regions;

	PairList::HashPair **hash_pair = pair_list->hash_pair;
	table = TABLE(ipair);
	index = INDEX(ipair);

	nx = n[0];
	ny = n[1];
	nz = n[2];

	// Overlap
	drijn = delta;
	drijnx = drijn * nx;
	drijny = drijn * ny;
	drijnz = drijn * nz;
	
	// may need to consider wall translational velocity in the future
	double vel_wall[3];
	for (int i = 0; i < 3; i++) vel_wall[i] = 0.0;

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

	if (iparticle == 0 && update->ntimestep > 172700 && iwall == 5) {
		iparticle = iparticle;
	}

	// vijt = vij - (vij . nij) . nij
	vijtx = vijx - vijnx;
	vijty = vijy - vijny;
	vijtz = vijz - vijnz;
	vijt = sqrt(vijtx*vijtx + vijty*vijty + vijtz*vijtz);
	if (fabs(vijt) < EPSILON) vijt_inv = 0.0;
	else vijt_inv = 1.0/vijt;

	// tij = tij/|tij|
	tx = vijtx * vijt_inv;
	ty = vijty * vijt_inv;
	tz = vijtz * vijt_inv;
		
	// calculate drijt
	// first time to contact
	int contact_flag = hash_pair[table][index].ivalues[0];
	double temp;
	if (contact_flag == 0) {
		if (fabs(vijn) < EPSILON) temp = dt;
		else temp = MIN(fabs(drijn/vijn),dt);
		drijtx = vijtx * temp;
		drijty = vijty * temp;
		drijtz = vijtz * temp;
		hash_pair[table][index].ivalues[0] = 1;
	}
	// during the same contact
	else {
		drijtx = hash_pair[table][index].dvalues[0];
		drijty = hash_pair[table][index].dvalues[1];
		drijtz = hash_pair[table][index].dvalues[2];
		
		// update the tang. disp. for the next time step
		drijtx = drijtx + vijtx*dt;
		drijty = drijty + vijty*dt;
		drijtz = drijtz + vijtz*dt;

		temp = drijtx*nx + drijty*ny + drijtz*nz;

		drijtx = drijtx - temp*nx;
		drijty = drijty - temp*ny;
		drijtz = drijtz - temp*nz;
	}

	hash_pair[table][index].dvalues[0] = drijtx;
	hash_pair[table][index].dvalues[1] = drijty;
	hash_pair[table][index].dvalues[2] = drijtz;
	
	// fijn = kn*delta_rijn - Cn*vijn
	fijnx = kn*drijnx - Cn*vijnx;
	fijny = kn*drijny - Cn*vijny;
	fijnz = kn*drijnz - Cn*vijnz;

	fijn = sqrt(fijnx*fijnx + fijny*fijny + fijnz*fijnz);
	
	// tangential force
	// fijt = -kt*delta_rijt - Ct*vijt
	fijtx = -kt * drijtx - Ct * vijtx;
	fijty = -kt * drijty - Ct * vijty;
	fijtz = -kt * drijtz - Ct * vijtz;

	fijt = sqrt(fijtx*fijtx + fijty*fijty + fijtz*fijtz); 
	temp = fabs(mu*fijn);
	// |fijt| > mu*|fnij|
	if (fijt > temp) {
		fijtx = -temp * tx;
		fijty = -temp * ty;
		fijtz = -temp * tz;
	}

	f[iparticle][0] += (fijnx + fijtx);
	f[iparticle][1] += (fijny + fijty);
	f[iparticle][2] += (fijnz + fijtz);

	if (particle->torque_flag) {
		fijt_a[0] = fijtx;
		fijt_a[1] = fijty;
		fijt_a[2] = fijtz;
		Vec_Cross_Prod_3D(torqueij, Li, fijt_a);

		torque[iparticle][0] += torqueij[0];
		torque[iparticle][1] += torqueij[1];
		torque[iparticle][2] += torqueij[2];
	}
}
