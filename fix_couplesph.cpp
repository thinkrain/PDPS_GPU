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
#include "neighbor.h"
#include "neigh_list.h"
#include "fix_couplesph.h"
#include "particle.h"
#include "region.h"
#include "group.h"
#include "update.h"
#include "parallel.h"

using namespace PDPS_NS;
using namespace FixConst;

#define PI 3.1416

/* ---------------------------------------------------------------------- */

FixCouplesph::FixCouplesph(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 5) {
		error->all(FLERR,"Illegal fix couplesph command");
	}
	
	poro_flag = 0;
	poroset_flag = 0;
	cubic_flag = 0;
	quintic_flag = 0;
	comm_forward = 3;
	comm_reverse = 3;


	sgid = group->find_group(arg[1]);
	sgroupbit = group->bitmask[sgid];
	lgid = group->find_group(arg[3]);
	lgroupbit = group->bitmask[lgid];
	rho_ref = atof(arg[4]);
	if (strcmp(arg[5], "Cubic") == 0)
		cubic_flag = 1;
	else if (strcmp(arg[5], "Quintic") == 0)
		quintic_flag = 1;
	else
		error->all(FLERR, "Wrong Kernel function");
	h = atof(arg[6]);
	if (strcmp(arg[7], "poro") == 0)
		poro_flag = 1;
	else if (strcmp(arg[7], "poroset") == 0){
		poroset_flag = 1;
		poroset = atof(arg[8]);
	}
	if (cubic_flag == 1){
		a2D = 10.0 / 7.0 / PI / h / h;
		a3D = 1.0 / PI / h / h / h;
	}
	else if (quintic_flag == 1){
		a2D = 7.0 / 4.0 / PI / h / h;
		a3D = 21.0 / 16.0 / PI / h / h / h;
	}

}

/* ---------------------------------------------------------------------- */

FixCouplesph::~FixCouplesph()
{

}

/* ---------------------------------------------------------------------- */

int FixCouplesph::setmask()
{
	int mask = 0;
	mask |= POST_FORCE;
	return mask;
}

/* ---------------------------------------------------------------------- */

void FixCouplesph::init()
{

}

/* ---------------------------------------------------------------------- */

void FixCouplesph::setup()
{
	post_force();
}

/* ---------------------------------------------------------------------- */

void FixCouplesph::post_force()
{
	int i, j, ii, jj, inum, jnum, itype, jtype;
	double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

	int *ilist, *jlist, *numneigh, **firstneigh;
	double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc;
	double rsq, tmp, wfd, delVdotDelR, mu, deltaE;

	double **v = particle->vest;
	double **x = particle->x;
	double **f = particle->f;
	double *rho = particle->rho;
	int *mask = particle->mask;
	double *mass = particle->mass;
	double *de = particle->de;
	double *drho = particle->drho;
	double *poro = particle->poro;
	double *radius = particle->radius;
	double *volume = particle->volume;
	double *hlocal = particle->hlocal;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	double wf, q, rij_inv;
	//  int newton_pair = force->newton_pair;

	// check consistency of pair coefficients

	inum = neighbor->neighlist->inum;
	ilist = neighbor->neighlist->ilist;
	numneigh = neighbor->neighlist->numneigh;
	firstneigh = neighbor->neighlist->firstneigh;

	//  update the force of ghost particles
	parallel->forward_force();
	//	calcuate the coupling force on other phase
	double weight;
	double fx, fy, fz;
	for (ii = 0; ii < nlocal; ii++) {
		i = ilist[ii];

		itype = type[i];
		if (mask[i] & sgroupbit){
			jlist = firstneigh[i];
			jnum = numneigh[i];
			imass = mass[itype];
			xtmp = x[i][0];
			ytmp = x[i][1];
			ztmp = x[i][2];
			weight = 0.0;
			fx = 0.0;
			fy = 0.0;
			fz = 0.0;
			for (jj = 0; jj < jnum; jj++) {
				j = jlist[jj];
				//     j &= NEIGHMASK;
				jtype = type[j];
				if (mask[j] & lgroupbit){
					jmass = mass[jtype];
					delx = xtmp - x[j][0];
					dely = ytmp - x[j][1];
					delz = ztmp - x[j][2];
					rsq = delx * delx + dely * dely + delz * delz;
					if (rsq < 4 * h * h) {
						q = sqrt(rsq) / h;

						if (cubic_flag == 1){
							if (q < 1)
								wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
							else
								wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
						}
						else if (quintic_flag == 1)
							wf = (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0)     * (2 * q + 1);

						if (domain->dim == 3)
							wf = wf * a3D;
						else
							wf = wf * a2D;

						fx += f[j][0] * wf;
						fy += f[j][1] * wf;
						fz += f[j][2] * wf;
						weight += wf;
					}		//  rsq < cutsq



				}		//  setflag[j][j]

			}		//  jnum

			if (weight > 0.001)
			{

				if (poroset_flag == 1){
					fx = fx / poroset;
					fy = fy / poroset;
					fz = fz / poroset;
				}
				if (poro_flag == 1){
					fx = fx / (1 - poro[i]);
					fy = fy / (1 - poro[i]);
					fz = fz / (1 - poro[i]);

				}

				f[i][0] += fx / weight * volume[i] * rho_ref / jmass;
				f[i][1] += fy / weight * volume[i] * rho_ref / jmass;
				f[i][2] += fz / weight * volume[i] * rho_ref / jmass;

				


			}

		} // setflag[i][i]
	} // inum  
	//  if vflag_fdotr) virial_fdotr_compute();

	// eleminate the influence of previous result
	for (i = nlocal; i < nlocal + particle->nghost; i++)
	{
		f[i][0] = 0.0;
		f[i][1] = 0.0;
		f[i][2] = 0.0;
	}

	//  add the coupled force back on fluid particle
	for (ii = 0; ii < nlocal; ii++) {
		i = ilist[ii];

		itype = type[i];
		if (mask[i] & sgroupbit){
			jlist = firstneigh[i];
			jnum = numneigh[i];
			imass = mass[itype];
			xtmp = x[i][0];
			ytmp = x[i][1];
			ztmp = x[i][2];
			weight = 0.0;

			for (jj = 0; jj < jnum; jj++) {
				j = jlist[jj];
				//     j &= NEIGHMASK;
				jtype = type[j];
				if (mask[j] & lgroupbit){
					jmass = mass[jtype];
					delx = xtmp - x[j][0];
					dely = ytmp - x[j][1];
					delz = ztmp - x[j][2];
					rsq = delx * delx + dely * dely + delz * delz;
					if (rsq < 4 * h * h) {
						q = sqrt(rsq) / h;

						if (cubic_flag == 1){
							if (q < 1)
								wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
							else
								wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
						}
						else if (quintic_flag == 1)
							wf = (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0)     * (2 * q + 1);

						if (domain->dim == 3)
							wf = wf * a3D;
						else
							wf = wf * a2D;

						fx += f[j][0] * wf;
						fy += f[j][1] * wf;
						fz += f[j][2] * wf;
						weight += wf;
					}		//  rsq < cutsq



				}		//  setflag[j][j]

			}		//  jnum

			if (weight > 0.001){
				for (jj = 0; jj < jnum; jj++) {
					j = jlist[jj];
					//     j &= NEIGHMASK;
					jtype = type[j];
					jmass = mass[jtype];
					if (mask[j] & lgroupbit){
						delx = xtmp - x[j][0];
						dely = ytmp - x[j][1];
						delz = ztmp - x[j][2];
						rsq = delx * delx + dely * dely + delz * delz;
						if (rsq < 4 * h * h) {
							q = sqrt(rsq) / h;

							if (cubic_flag == 1){
								if (q < 1)
									wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
								else
									wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
							}
							else if (quintic_flag == 1)
								wf = (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0)     * (2 * q + 1);

							if (domain->dim == 3)
								wf = wf * a3D;
							else
								wf = wf * a2D;

							f[j][0] -= f[i][0] * wf / weight;
							f[j][1] -= f[i][1] * wf / weight;
							f[j][2] -= f[i][2] * wf / weight;

						}		//  rsq < cutsq



					}		//  setflag[j][j]

				}		//  jnum
			}


		}

	}
	parallel->reverse_comm();

	

}
