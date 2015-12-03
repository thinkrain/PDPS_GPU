/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "domain.h"
#include "error.h"
#include "fix_addforce.h"
#include "force.h"
#include "memory.h"
#include "pair.h"
#include "particle.h"
#include "phy_const.h"
#include "region.h"

using namespace PDPS_NS;
using namespace FixConst;
using namespace PhyConst;

#define EPSILON 1.0e-10

enum{DRAG_STOKES, DRAG_GENERAL, DRAG_FELICE, BUOYANCY, CUSTOM};

/* ---------------------------------------------------------------------- */

FixAddForce::FixAddForce(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 5) error->all(FLERR,"Illegal fix setforce command");

	fx = fy = fz = 0.0;

	coupled = 0;

	voidage = vol_solid = vol_solid_total = NULL;

	if (strcmp(arg[3], "drag/stokes") == 0) {
		force_style = DRAG_STOKES;
		mu = atof(arg[4]); 
		if (strcmp(arg[5], "coupled") == 0)
			coupled = 1;
	}
	else if (strcmp(arg[3], "buoyancy") == 0) {
		force_style = BUOYANCY;
		if (strcmp(arg[4], "x") == 0) {
			g_dim = 0;
		}
		else if (strcmp(arg[4], "y") == 0) {
			g_dim = 1;
		}
		else if (strcmp(arg[4], "z") == 0) {
			g_dim = 2;
		}
		else error->all(FLERR, "Illegal fix addforce buoyancy options");
		rho = atof(arg[5]);
		g = atof(arg[6]);
	}
	else if (!strcmp(arg[3], "drag/general")) {
		force_style = DRAG_GENERAL;
		mu = atof(arg[4]);
		rho = atof(arg[5]);
		eta = atof(arg[6]);
	}
	else if (!strcmp(arg[3], "drag/felice")) {
		force_style = DRAG_FELICE;
		mu = atof(arg[4]);
		rho = atof(arg[5]);
		cle[2] = atof(arg[6]);
		cle[0] = 0.0;
		cle[1] = 0.0;
		if (!strcmp(arg[7], "region")) {
			rid = domain->find_region(arg[8]);
			if (rid == -1) error->all(FLERR, "Cannot find region id");
		}
		cell[0] = cell[1] = 1;
		cell[2] = static_cast<int> (domain->regions[rid]->extent_le[2] / cle[2]);
		voidage = new double[cell[2]];
		vol_solid = new double[cell[2]];
		vol_solid_total = new double[cell[2]];
	}
	else if (strcmp(arg[3], "custom") == 0) {
		force_style = CUSTOM;
	}
	else error->all(FLERR, "Illegal fix addforce keyword");
}

/* ---------------------------------------------------------------------- */

FixAddForce::~FixAddForce()
{
	
}

/* ---------------------------------------------------------------------- */

int FixAddForce::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddForce::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixAddForce::post_force()
{
	if (force_style == DRAG_STOKES) {
		add_drag_stokes();
	}
	else if (force_style == DRAG_GENERAL) {
		add_drag_general();
	}
	else if (force_style == DRAG_FELICE) {
		add_drag_felice();
	}
	else if (force_style == BUOYANCY) {
		add_buoyancy();
	}
}

/* ---------------------------------------------------------------------- */

void FixAddForce::add_drag_stokes()
{
	double **v = particle->v;
	double **f = particle->f;
	double **x = particle->x;
	double *mass = particle->mass;
	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;

	double R; 
	int pair_id, itype;
	double coeff = -6 * PI * mu;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			itype = type[i];
			pair_id = force->type2pair[itype][itype];
			R = 0.5 * force->pair[pair_id]->cut[itype][itype];
			f[i][0] += coeff * R * v[i][0];
			f[i][1] += coeff * R * v[i][1];
			f[i][2] += coeff * R * v[i][2];
		}
	}

	if (coupled = 1){
		int *ilist, *jlist, *numneigh, **firstneigh;
		int inum, jnum, jtype, ii, jj, i, j;
		double imass, h, ih, ihsq, wf;
		double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
		inum = neighbor->neighlist->inum;
		ilist = neighbor->neighlist->ilist;
		numneigh = neighbor->neighlist->numneigh;
		firstneigh = neighbor->neighlist->firstneigh;


			// add density at each particle via kernel function overlap
		for (ii = 0; ii < inum; ii++) {
			i = ilist[ii];
			if (mask[i] & groupbit){
				xtmp = x[i][0];
				ytmp = x[i][1];
				ztmp = x[i][2];
				itype = type[i];
				jlist = firstneigh[i];
				jnum = numneigh[i];

				for (jj = 0; jj < jnum; jj++) {
					j = jlist[jj];

					jtype = type[j];
					delx = xtmp - x[j][0];
					dely = ytmp - x[j][1];
					delz = ztmp - x[j][2];
					rsq = delx * delx + dely * dely + delz * delz;
					pair_id = force->type2pair[itype][itype];
					h = force->pair[pair_id]->cut[itype][itype];
					if (rsq < h * h) {
						ih = 1.0 / h;
						ihsq = ih * ih;

						if (domain->dim == 3) {
							/*
							// Lucy kernel, 3d
							r = sqrt(rsq);
							wf = (h - r) * ihsq;
							wf =  2.0889086280811262819e0 * (h + 3. * r) * wf * wf * wf * ih;
							*/

							// quadric kernel, 3d
							wf = 1.0 - rsq * ihsq;
							wf = wf * wf;
							wf = wf * wf;
							wf = 2.1541870227086614782e0 * wf * ihsq * ih;
						}
						else {
							// Lucy kernel, 2d
							//r = sqrt(rsq);
							//wf = (h - r) * ihsq;
							//wf = 1.5915494309189533576e0 * (h + 3. * r) * wf * wf * wf;

							// quadric kernel, 2d
							wf = 1.0 - rsq * ihsq;
							wf = wf * wf;
							wf = wf * wf;
							wf = 1.5915494309189533576e0 * wf * ihsq;
							//wf = 0.9 * wf * ihsq;
						}
						R = 0.5 * force->pair[pair_id]->cut[itype][itype];
						f[j][0] -= coeff * R * v[i][0];
						f[j][1] -= coeff * R * v[i][1];
						f[j][2] -= coeff * R * v[i][2];

					}

				}
			}
		

		}
	}
}

/* ---------------------------------------------------------------------- */
// This part needs further investigation

void FixAddForce::add_drag_general()
{
	double **v = particle->v;
	double **f = particle->f;
	double *mass = particle->mass;
	int *mask = particle->mask;
	double *radius = particle->radius;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	double Cd;

	int rid = domain->find_region("vessel");
	if (rid == -1) error->all(FLERR, "Cannot find the region\n");
	double **v_coords = domain->regions[rid]->v_coords;

	double Re, coeff, xi;
	double vel_f, v_f[3];
	double radius_sq;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			if (particle->radius_flag) {
				for (int j = 0; j < 3; j++) v_f[j] = v[i][j] - v_coords[0][j];
				vel_f = (v_f[0]*v_f[0] + v_f[1]*v_f[1] + v_f[2]*v_f[2]);
				vel_f = sqrt(vel_f);
				Re = eta * rho * vel_f * 2 * radius[i] / mu;
				if (Re > EPSILON) {
					if (Re < 0.5) Cd = 24 / Re;
					else {
						Cd = (0.63 + 4.8 / sqrt(Re));
						Cd = Cd * Cd;
					}
					xi = 3.7 - 0.65*exp(-(1.5 - log10(Re))*(1.5 - log10(Re)) / 2);
					// Fd = 0.5 * rho * vel_relative^2 * Cd * pi * R^2
					radius_sq = radius[i] * radius[i];
					coeff = 0.5 * rho * PI * radius_sq * pow(eta, 1 - xi);
					f[i][0] += (-coeff * Cd * v_f[0] * vel_f);
					f[i][1] += (-coeff * Cd * v_f[1] * vel_f);
					f[i][2] += (-coeff * Cd * v_f[2] * vel_f);
				}
			}
		}
	}
}

/* ---------------------------------------------------------------------- */
// This part needs further investigation

void FixAddForce::add_drag_felice()
{
	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	double *mass = particle->mass;
	int *mask = particle->mask;
	double *radius = particle->radius;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	double Cd;
	
	double **coords = domain->regions[rid]->coords;
	double **v_coords = domain->regions[rid]->v_coords;

	compute_voidage();

	double Re, coeff, xi;
	double vel_f, v_f[3];
	double radius_sq;
	int cid, inside_flag;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			if (particle->radius_flag) {
				for (int j = 0; j < 3; j++) v_f[j] = v[i][j] - v_coords[0][j];
				vel_f = (v_f[0] * v_f[0] + v_f[1] * v_f[1] + v_f[2] * v_f[2]);
				vel_f = sqrt(vel_f);
				inside_flag = domain->regions[rid]->inside(x[i]);
				if (inside_flag == 0) {
					continue;
				}
				cid = static_cast<int> ((x[i][2]-coords[0][2]) / cle[2]);
				if (cid == cell[2]) cid--;
				if (cid < 0 || cid > cell[2]) {
					error->warning(FLERR, "Check code: compute_voidage()");
				}
				Re = voidage[cid] * rho * vel_f * 2 * radius[i] / mu;
				if (Re > EPSILON) {
					if (Re < 0.5) Cd = 24 / Re;
					else {
						Cd = (0.63 + 4.8 / sqrt(Re));
						Cd = Cd * Cd;
					}
					xi = 3.7 - 0.65*exp(-(1.5 - log10(Re))*(1.5 - log10(Re)) / 2);
					// Fd = 0.5 * rho * vel_relative^2 * Cd * pi * R^2
					radius_sq = radius[i] * radius[i];
					coeff = 0.5 * rho * PI * radius_sq * pow(voidage[cid], 1 - xi);
					f[i][0] += (-coeff * Cd * v_f[0] * vel_f);
					f[i][1] += (-coeff * Cd * v_f[1] * vel_f);
					f[i][2] += (-coeff * Cd * v_f[2] * vel_f);
				}
			}
		}
	}
}

/* ---------------------------------------------------------------------- */

void FixAddForce::compute_voidage()
{
	int *mask = particle->mask;
	double **x = particle->x;
	double *radius = particle->radius;
	int nlocal = particle->nlocal;
	double **coords;

	for (int i = 0; i < cell[2]; i++) {
		voidage[i] = 0.0;
		vol_solid[i] = 0.0;
		vol_solid_total[i] = 0.0;
	}

	coords = domain->regions[rid]->coords;
	int cid, inside_flag;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			inside_flag = domain->regions[rid]->inside(x[i]);
			if (inside_flag == 0) {
				continue;
			}
			cid = static_cast<int> ((x[i][2]-coords[0][2]) / cle[2]);
			if (cid == cell[2]) cid--;
			if (cid < 0 || cid > cell[2]) {
				error->warning(FLERR, "Check code: compute_voidage()");
			}
			vol_solid[cid] += 4.0 / 3 * PI * radius[i] * radius[i] * radius[i];
		}
	}
	MPI_Allreduce(&vol_solid[0], &vol_solid_total[0], cell[2], MPI_DOUBLE, MPI_SUM, mworld);

	double vol_cell = PI * domain->regions[rid]->radius * domain->regions[rid]->radius * cle[2];
	double vol_last_cell;
	for (int i = 0; i < cell[2]; i++) {
		if (i < cell[2] - 1) voidage[i] = 1 - vol_solid_total[i] / vol_cell;
		else {
			vol_last_cell = vol_cell / cle[2] * (domain->regions[rid]->extent_le[2] - (cell[2] - 1)*cle[2]);
			voidage[i] = 1 - vol_solid_total[i] / vol_cell;
		}
		// for some special case, when a very large big particle's center is inside one cell
		if (voidage[i] < 0.0) voidage[i] = 0.0;   
	}
}

/* ---------------------------------------------------------------------- */

void FixAddForce::add_buoyancy()
{
	double **v = particle->v;
	double **f = particle->f;
	double *mass = particle->mass;
	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;

	double R; 
	int pair_id, itype;
	double coeff = 4.0 / 3 * PI * rho *g;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			itype = type[i];
			pair_id = force->type2pair[itype][itype];
			R = 0.5 * force->pair[pair_id]->cut[itype][itype];
	//		R = particle->radius[itype];
			f[i][g_dim] += coeff * R * R * R;
		}
	}
}
