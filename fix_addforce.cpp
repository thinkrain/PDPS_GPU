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
#include "group.h"
#include "memory.h"
#include "pair.h"
#include "particle.h"
#include "phy_const.h"
#include "region.h"
#include "pair.h"

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
		if (strcmp(arg[5], "coupled") == 0){
			coupled = 1;
			int itype = atoi(arg[6]);
			int jtype = atoi(arg[7]);
			cutoff = atof(arg[8]);
		}
			
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
		rho_ref = atof(arg[5]);
		g = atof(arg[6]);
		
	}
	else if (!strcmp(arg[3], "drag/general")) {
		force_style = DRAG_GENERAL;
		mu = atof(arg[4]);
		rho_ref = atof(arg[5]);
		eta = atof(arg[6]);
		if (strcmp(arg[7], "coupled") == 0){
			coupled = 1;
		}
	}
	else if (!strcmp(arg[3], "drag/felice")) {
		force_style = DRAG_FELICE;
		lgid = group->find_group(arg[4]);
		lgroupbit = group->bitmask[lgid];
		mu = atof(arg[5]);
		rho_ref = atof(arg[6]);
		voi_ref = atof(arg[7]);
		if (!strcmp(arg[8], "voiset")){
			voiset_flag = 1;
			voiset = atof(arg[9]);
			if (strcmp(arg[10], "coupled") == 0){
				coupled = 1;
			}
		}
		else if (!strcmp(arg[8], "voiporo")){
			voiset_flag = 0;
			if (strcmp(arg[9], "coupled") == 0){
				coupled = 1;
				h = atof(arg[10]);
			}
		}
		else
			error->all(FLERR, "Illegal fix addforce drag/felice options");

			
//		cle[2] = atof(arg[6]);
//		cle[0] = 0.0;
//		cle[1] = 0.0;
//		if (strcmp(arg[7], "region") == 0) {
//			rid = domain->find_region(arg[8]);
//			if (rid == -1) error->all(FLERR, "Cannot find region id");
//		}
//		cell[0] = cell[1] = 1;
//		cell[2] = static_cast<int> (domain->regions[rid]->extent_le[2] / cle[2]);
//		voidage = new double[cell[2]];
//		vol_solid = new double[cell[2]];
//		vol_solid_total = new double[cell[2]];
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
	double *radius = particle->radius;
	double *volume = particle->volume;

	int itype;
	double coeff = -6 * PI * mu;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			itype = type[i];
			f[i][0] += coeff * radius[i] * v[i][0];
			f[i][1] += coeff * radius[i] * v[i][1];
			f[i][2] += coeff * radius[i] * v[i][2];
		}
	}

	if (coupled == 1){
		int *ilist, *jlist, *numneigh, **firstneigh;
		int inum, jnum, jtype, ii, jj, i, j;
		double imass, h, q, wf;
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
					if (mask[i] == mask[j])
						continue;
					jtype = type[j];
					delx = xtmp - x[j][0];
					dely = ytmp - x[j][1];
					delz = ztmp - x[j][2];
					rsq = delx * delx + dely * dely + delz * delz;
					//h = force->pair[pair_id]->cut[itype][itype];
					h = radius[i] / 2.0;
					if (rsq < radius[i] * radius[i]) {

						q = sqrt(rsq) / h;

						if (q < 1)
							wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
						else
							wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
						if (domain->dim == 3)
							wf = wf * 1.0 / PI / h / h / h;
						else
							wf = wf * 10.0 / 7.0 / PI / h / h;

						f[j][0] -= coeff * radius[i] * v[i][0] * volume[i] * wf;
						f[j][1] -= coeff * radius[i] * v[i][1] * volume[i] * wf;
						f[j][2] -= coeff * radius[i] * v[i][2] * volume[i] * wf;

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
	double **x = particle->x;
	double **f = particle->f;
	double *mass = particle->mass;
	int *mask = particle->mask;
	double *radius = particle->radius;
	double *volume = particle->volume;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	double Cd;

//	int rid = domain->find_region("vessel");
//	if (rid == -1) error->all(FLERR, "Cannot find the region\n");
//	double **v_coords = domain->regions[rid]->v_coords;

	double Re, coeff, xi;
	double vel_f, v_f[3];
	double radius_sq;
	if (coupled == 0){
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				//		if (particle->radius_flag) {
				for (int j = 0; j < 3; j++)
					v_f[j] = v[i][j];
				//	v_f[j] = v[i][j] - v_coords[0][j];
				vel_f = (v_f[0] * v_f[0] + v_f[1] * v_f[1] + v_f[2] * v_f[2]);
				vel_f = sqrt(vel_f);
				Re = eta * rho_ref * vel_f * 2 * radius[i] / mu;
				if (Re > EPSILON) {
					if (Re < 0.5) Cd = 24 / Re;
					else {
						Cd = (0.63 + 4.8 / sqrt(Re));
						Cd = Cd * Cd;
					}
					//				xi = 3.7 - 0.65*exp(-(1.5 - log10(Re))*(1.5 - log10(Re)) / 2);
					// Fd = 0.5 * rho_ref * vel_relative^2 * Cd * pi * R^2
					radius_sq = radius[i] * radius[i];
					//				coeff = 0.5 * rho_ref * PI * radius_sq * pow(eta, 1 - xi);
					coeff = 0.5 * rho_ref * PI * radius_sq;
					f[i][0] += (-coeff * Cd * v_f[0] * vel_f);
					f[i][1] += (-coeff * Cd * v_f[1] * vel_f);
					f[i][2] += (-coeff * Cd * v_f[2] * vel_f);
				}
				//		}
			}
		}
	}
	else if (coupled == 1){
		int *ilist, *jlist, *numneigh, **firstneigh;
		int inum, jnum, jtype, ii, jj, i, j, itype;
		double imass, h, q, wf;
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

				for (int j = 0; j < 3; j++)
					v_f[j] = v[i][j];
				//	v_f[j] = v[i][j] - v_coords[0][j];
				vel_f = (v_f[0] * v_f[0] + v_f[1] * v_f[1] + v_f[2] * v_f[2]);
				vel_f = sqrt(vel_f);
				Re = eta * rho_ref * vel_f * 2 * radius[i] / mu;
				if (Re > EPSILON) {
					if (Re < 0.5) Cd = 24 / Re;
					else {
						Cd = (0.63 + 4.8 / sqrt(Re));
						Cd = Cd * Cd;
					}
					//				xi = 3.7 - 0.65*exp(-(1.5 - log10(Re))*(1.5 - log10(Re)) / 2);
					// Fd = 0.5 * rho_ref * vel_relative^2 * Cd * pi * R^2
					radius_sq = radius[i] * radius[i];
					//				coeff = 0.5 * rho_ref * PI * radius_sq * pow(eta, 1 - xi);
					coeff = 0.5 * rho_ref * PI * radius_sq;
					f[i][0] += (-coeff * Cd * v_f[0] * vel_f);
					f[i][1] += (-coeff * Cd * v_f[1] * vel_f);
					f[i][2] += (-coeff * Cd * v_f[2] * vel_f);
				}

				for (jj = 0; jj < jnum; jj++) {
					j = jlist[jj];
					if (mask[i] == mask[j])
						continue;
					jtype = type[j];
					delx = xtmp - x[j][0];
					dely = ytmp - x[j][1];
					delz = ztmp - x[j][2];
					rsq = delx * delx + dely * dely + delz * delz;
					//h = force->pair[pair_id]->cut[itype][itype];
					h = radius[i] / 2.0;
					if (rsq < radius[i] * radius[i]) {

						q = sqrt(rsq) / h;

						if (q < 1)
							wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
						else
							wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
						if (domain->dim == 3)
							wf = wf * 1.0 / PI / h / h / h;
						else
							wf = wf * 10.0 / 7.0 / PI / h / h;

						f[j][0] -= (-coeff * Cd * v_f[0] * vel_f) * volume[i] * wf;
						f[j][1] -= (-coeff * Cd * v_f[1] * vel_f) * volume[i] * wf;
						f[j][2] -= (-coeff * Cd * v_f[2] * vel_f) * volume[i] * wf;

					}

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
	double *rho = particle->rho;
	double *poro = particle->poro;
	double *volume = particle->volume;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	double Cd;
	
//	double **coords = domain->regions[rid]->coords;
//	double **v_coords = domain->regions[rid]->v_coords;

	//compute_voidage();

	double Re, coeff, xi;
	double vel_f, v_f[3];
	double radius_sq;
	int cid, inside_flag;
	if (coupled == 0){
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				for (int j = 0; j < 3; j++) v_f[j] = v[i][j];// -v_coords[0][j];
				vel_f = (v_f[0] * v_f[0] + v_f[1] * v_f[1] + v_f[2] * v_f[2]);
				vel_f = sqrt(vel_f);
//				inside_flag = domain->regions[rid]->inside(x[i]);
//				if (inside_flag == 0) {
//					continue;
//				}
				//	cid = static_cast<int> ((x[i][2]-coords[0][2]) / cle[2]);
				//	if (cid == cell[2]) cid--;
				//	if (cid < 0 || cid > cell[2]) {
				//		error->warning(FLERR, "Check code: compute_voidage()");
				//	}
				Re = rho_ref * vel_f * 2 * radius[i] / mu;
				if (Re > EPSILON) {
					if (Re < 0.5) Cd = 24 / Re;
					else {
						Cd = (0.63 + 4.8 / sqrt(Re));
						Cd = Cd * Cd;
					}
					voi = rho[i] / rho_ref;
					xi = 3.7 - 0.65*exp(-(1.5 - log10(Re))*(1.5 - log10(Re)) / 2);
					// Fd = 0.5 * rho_ref * vel_relative^2 * Cd * pi * R^2
					radius_sq = radius[i] * radius[i];
					if (voi > voi_ref){
						coeff = 0.5 * rho_ref * PI * radius_sq * pow(voi, -xi);
						f[i][0] += (-coeff * Cd * v_f[0] * vel_f);
						f[i][1] += (-coeff * Cd * v_f[1] * vel_f);
						f[i][2] += (-coeff * Cd * v_f[2] * vel_f);
					}
				}
			}
		}
	}
	else if (coupled == 1){
		int *ilist, *jlist, *numneigh, **firstneigh;
		int inum, jnum, jtype, ii, jj, i, j, itype;
		double imass, q, wf, wfsum;
		double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
		inum = neighbor->neighlist->inum;
		ilist = neighbor->neighlist->ilist;
		numneigh = neighbor->neighlist->numneigh;
		firstneigh = neighbor->neighlist->firstneigh;

		double v_ref[3];
	
		// add density at each particle via kernel function overlap
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				xtmp = x[i][0];
				ytmp = x[i][1];
				ztmp = x[i][2];
				itype = type[i];
				jlist = firstneigh[i];
				jnum = numneigh[i];
				v_ref[0] = 0.0;
				v_ref[1] = 0.0;
				v_ref[2] = 0.0;
				wfsum = 0.0;
				for (jj = 0; jj < jnum; jj++) {
					j = jlist[jj];
					if (mask[j] & lgroupbit){
						jtype = type[j];
						delx = xtmp - x[j][0];
						dely = ytmp - x[j][1];
						delz = ztmp - x[j][2];
						rsq = delx * delx + dely * dely + delz * delz;
						//h = force->pair[pair_id]->cut[itype][itype];
						//h = radius[i] / 2.0;
						if (rsq < 4 * h * h) {

							q = sqrt(rsq) / h;

							if (q < 1)
								wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
							else
								wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
							if (domain->dim == 3)
								wf = wf * 1.0 / PI / h / h / h;
							else
								wf = wf * 10.0 / 7.0 / PI / h / h;

							v_ref[0] += v[j][0] * wf;
							v_ref[1] += v[j][1] * wf;
							v_ref[2] += v[j][2] * wf;
							wfsum += wf;

						}
					}
				}
				if (wfsum > EPSILON){
					v_ref[0] = v_ref[0] / wfsum;
					v_ref[1] = v_ref[1] / wfsum;
					v_ref[2] = v_ref[2] / wfsum;

					for (j = 0; j < 3; j++) v_f[j] = v[i][j] - v_ref[j];// -v_coords[0][j];
					vel_f = (v_f[0] * v_f[0] + v_f[1] * v_f[1] + v_f[2] * v_f[2]);
					vel_f = sqrt(vel_f);
					//				inside_flag = domain->regions[rid]->inside(x[i]);
					//		if (inside_flag == 0) {
					//				continue;
					//			}
					//	cid = static_cast<int> ((x[i][2]-coords[0][2]) / cle[2]);
					//	if (cid == cell[2]) cid--;
					//	if (cid < 0 || cid > cell[2]) {
					//		error->warning(FLERR, "Check code: compute_voidage()");
					//	}
					Re = rho_ref * vel_f * 2 * radius[i] / mu;
					if (Re > EPSILON) {
						if (Re < 0.5) Cd = 24 / Re;
						else {
							Cd = (0.63 + 4.8 / sqrt(Re));
							Cd = Cd * Cd;
						}
						//voi = 0.8;
						//voi = rho[i] / rho_ref;
						voi = 1 - poro[i];
						if (voiset_flag == 1)
							voi = voiset;
						xi = 3.7 - 0.65*exp(-(1.5 - log10(Re))*(1.5 - log10(Re)) / 2);
						// Fd = 0.5 * rho_ref * vel_relative^2 * Cd * pi * R^2
						radius_sq = radius[i] * radius[i];
						if (voi > voi_ref){
							coeff = 0.5 * rho_ref * PI * radius_sq * pow(voi, -xi);
							f[i][0] += (-coeff * Cd * v_f[0] * vel_f);
							f[i][1] += (-coeff * Cd * v_f[1] * vel_f);
							f[i][2] += (-coeff * Cd * v_f[2] * vel_f);

							for (jj = 0; jj < jnum; jj++) {
								j = jlist[jj];
								if (mask[j] & lgroupbit){
									jtype = type[j];
									delx = xtmp - x[j][0];
									dely = ytmp - x[j][1];
									delz = ztmp - x[j][2];
									rsq = delx * delx + dely * dely + delz * delz;
									//h = force->pair[pair_id]->cut[itype][itype];
									//h = radius[i] / 2.0;
									if (rsq < 4 * h * h) {

										q = sqrt(rsq) / h;

										if (q < 1)
											wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
										else
											wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
										if (domain->dim == 3)
											wf = wf * 1.0 / PI / h / h / h;
										else
											wf = wf * 10.0 / 7.0 / PI / h / h;

										f[j][0] -= (-coeff * Cd * v_f[0] * vel_f)  * wf / wfsum;
										f[j][1] -= (-coeff * Cd * v_f[1] * vel_f)  * wf / wfsum;
										f[j][2] -= (-coeff * Cd * v_f[2] * vel_f)  * wf / wfsum;

									}
								}


							}
						}



					}	// Re number

				}	//  if there is sph particle inside the cutoff distance

			}	// groupbit
		}	//	nlocal


	}	// coupled == 1
	
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
	double coeff = 4.0 / 3.0 * PI * rho_ref *g;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			R = particle->radius[i];
		//	if (R < EPSILON){
		//		itype = type[i];
		//		pair_id = force->type2pair[itype][itype];
		//		particle->radius[i] = 0.5 * force->pair[pair_id]->cut[itype][itype];
		//		R = particle->radius[i];
		//	}
		
			f[i][g_dim] += coeff * R * R * R;
		}
	}
}
