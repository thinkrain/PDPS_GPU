/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "parallel.h"
#include "error.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair_sph_rhosum.h"
#include "particle.h"
#include "random_mars.h"
#include "update.h"

using namespace PDPS_NS;


#define DELTA 1
#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairSPH_RHOSUM::PairSPH_RHOSUM(PDPS *ps) : Pair(ps)
{
	first = 1;
	allocated = 0;
	cutr = NULL;
	cutrsq = NULL;
}

/* ---------------------------------------------------------------------- */

PairSPH_RHOSUM::~PairSPH_RHOSUM()
{
	if (allocated) {
		memory->destroy(setflag);
		memory->destroy(cutrsq);

		memory->destroy(cutr);

	}

}

/* ---------------------------------------------------------------------- */

void PairSPH_RHOSUM::allocate()
{
	allocated = 1;
	int n = particle->ntypes;

	memory->create(setflag, n + 1, n + 1, "pair:setflag");
	for (int i = 1; i <= n; i++)
	for (int j = i; j <= n; j++)
		setflag[i][j] = 0;
	memory->create(cutrsq, n + 1, n + 1, "pair:cutrsq");

	memory->create(cutr, n + 1, n + 1, "pair:cutr");


}

/* ----------------------------------------------------------------------
Compute force for all paritcles
------------------------------------------------------------------------- */

void PairSPH_RHOSUM::compute(int eflag, int vflag)
{
	int i, j, ii, jj, jnum, itype, jtype;
	double xtmp, ytmp, ztmp, delx, dely, delz;
	double rsq, imass, h, ih, ihsq;
	int *jlist;
	double wf;
	// neighbor list variables
	int inum, *ilist, *numneigh, **firstneigh;

	if (eflag || vflag)
		ev_setup(eflag, vflag);


	double **x = particle->x;
	double *rho = particle->rho;
	int *type = particle->type;
	double *mass = particle->mass;

	// check consistency of pair coefficients

	if (first) {
		for (i = 1; i <= particle->ntypes; i++) {
			for (j = 1; i <= particle->ntypes; i++) {
				if (cutrsq[i][j] > 0.0) {
					if (!setflag[i][i] || !setflag[j][j]) {
						if (parallel->procid == 0) {
							printf(
								"SPH particle types %d and %d interact, but not all of their single particle properties are set.\n",
								i, j);
						}
					}
				}
			}
		}
		first = 0;
	}

	inum = neighbor->neighlist->inum;
	ilist = neighbor->neighlist->ilist;
	numneigh = neighbor->neighlist->numneigh;
	firstneigh = neighbor->neighlist->firstneigh;

	// recompute density
	// we use a full neighborlist here

	if (nstep != 0) {
		if ((update->ntimestep % nstep) == 0) {

			// initialize density with self-contribution,
			for (ii = 0; ii < inum; ii++) {
				i = ilist[ii];
				itype = type[i];
				imass = mass[itype];

				h = cutr[itype][itype];
				if (domain->dim == 3) {
					/*
					// Lucy kernel, 3d
					wf = 2.0889086280811262819e0 / (h * h * h);
					*/

					// quadric kernel, 3d
					wf = 2.1541870227086614782 / (h * h * h);
				}
				else {
					/*
					// Lucy kernel, 2d
					wf = 1.5915494309189533576e0 / (h * h);
					*/

					// quadric kernel, 2d
					wf = 1.5915494309189533576e0 / (h * h);
				}

				rho[i] = imass * wf;
			}

			// add density at each particle via kernel function overlap
			for (ii = 0; ii < inum; ii++) {
				i = ilist[ii];
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

					if (rsq < cutrsq[itype][jtype]) {
						h = cutr[itype][jtype];
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
						}

						rho[i] += mass[jtype] * wf;
					}

				}
			}
		}
	}


}



/* ----------------------------------------------------------------------
Setting for pair_style command
------------------------------------------------------------------------- */

void PairSPH_RHOSUM::set_style(int narg, char **arg)
{
	if (narg != 2)
		error->all(FLERR, "Illegal number of setting arguments for pair_style sph/RHOSUM");
	nstep = atoi(arg[1]);

}

/* ----------------------------------------------------------------------
Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairSPH_RHOSUM::set_coeff(int narg, char **arg)
{
	if (narg != 3)
		error->all(FLERR, "Incorrect number of args for sph/rhosum coefficients");
	if (!allocated)
		allocate();

	int ilo, ihi, jlo, jhi;
	force->bounds(arg[0], particle->ntypes, ilo, ihi);
	force->bounds(arg[1], particle->ntypes, jlo, jhi);

	double cut_one = atof(arg[2]);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			//printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
			cutr[i][j] = cut_one;
			setflag[i][j] = 1;
			count++;
		}
	}

	if (count == 0)
		error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairSPH_RHOSUM::init_one(int i, int j) {

	if (setflag[i][j] == 0) {
		error->all(FLERR, "All pair sph/rhosum coeffs are not set");
	}

	cutr[j][i] = cutr[i][j];


}

/* ---------------------------------------------------------------------- */

double PairSPH_RHOSUM::single(int i, int j, int itype, int jtype,
	double rsq, double factor_coul, double factor_lj, double &fforce) {
	fforce = 0.0;
	return 0.0;
}

int PairSPH_RHOSUM::pack_forward_comm(int n, int *list, double *buf,
	int pbc_flag, int *pbc) {
	int i, j, m;
	double *rho = particle->rho;

	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		buf[m++] = rho[j];
	}
	return m;
}

/* ---------------------------------------------------------------------- */

void PairSPH_RHOSUM::unpack_forward_comm(int n, int first, double *buf) {
	int i, m, last;
	double *rho = particle->rho;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++)
		rho[i] = buf[m++];
}
