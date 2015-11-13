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
#include "pair_sph_heatconduction.h"
#include "particle.h"
#include "random_mars.h"
#include "update.h"

using namespace PDPS_NS;


#define DELTA 1
#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairSPH_HEATCONDUCTION::PairSPH_HEATCONDUCTION(PDPS *ps) : Pair(ps)
{
	 newton_pair = 1;
	 allocated = 0;
	 newton_pair = 1;
	 cut = NULL;
	 cutsq = NULL;
	 first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPH_HEATCONDUCTION::~PairSPH_HEATCONDUCTION()
{
	if (allocated) {
		memory->destroy(setflag);
		memory->destroy(cutsq);
		memory->destroy(cut);
		memory->destroy(alpha);
	}

}

/* ---------------------------------------------------------------------- */

void PairSPH_HEATCONDUCTION::allocate()
{
	 allocated = 1;
  int n = particle->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(alpha, n + 1, n + 1, "pair:alpha");

}

/* ----------------------------------------------------------------------
   Compute force for all paritcles
------------------------------------------------------------------------- */

void PairSPH_HEATCONDUCTION::compute(int eflag, int vflag)
{
	int i, j, ii, jj, inum, jnum, itype, jtype;
	double xtmp, ytmp, ztmp, delx, dely, delz;

	int *ilist, *jlist, *numneigh, **firstneigh;
	double imass, jmass, h, ih, ihsq;
	double rsq, wfd, D, deltaE;

	if (eflag || vflag)
		ev_setup(eflag, vflag);

	double **x = particle->x;
	double *e = particle->e;
	double *de = particle->de;
	double *mass = particle->mass;
	double *rho = particle->rho;
	int *type = particle->type;
	int nlocal = particle->nlocal;

	if (first) {
		for (i = 1; i <= particle->ntypes; i++) {
			for (j = 1; i <= particle->ntypes; i++) {
				if (cutsq[i][j] > 1.e-32) {
					if (!setflag[i][i] || !setflag[j][j]) {
						if (parallel->procid == 0) {
							printf(
								"SPH particle types %d and %d interact with cutoff=%g, but not all of their single particle properties are set.\n",
								i, j, sqrt(cutsq[i][j]));
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

	double wf;
	if (nstep != 0) {
		if ((update->ntimestep % nstep) == 0) {

			// initialize density with self-contribution,
			for (ii = 0; ii < inum; ii++) {
				i = ilist[ii];
				itype = type[i];
				imass = mass[itype];

				h = cut[itype][itype];
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
					//wf = 0.89 / (h * h);
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

					if (rsq < cutsq[itype][jtype]) {
						h = cut[itype][jtype];
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

						rho[i] += mass[jtype] * wf;
					}

				}
			}
		}
	}

	// loop over neighbors of my particles and do heat diffusion

	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		itype = type[i];

		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];

		jlist = firstneigh[i];
		jnum = numneigh[i];

		imass = mass[itype];

		for (jj = 0; jj < jnum; jj++) {
			j = jlist[jj];

			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			rsq = delx * delx + dely * dely + delz * delz;
			jtype = type[j];
			jmass = mass[jtype];

			if (rsq < cutsq[itype][jtype]) {
				h = cut[itype][jtype];
				ih = 1.0 / h;
				ihsq = ih * ih;

				// kernel function
				wfd = h - sqrt(rsq);
				if (domain->dim == 3) {
					// Lucy Kernel, 3d
					// Note that wfd, the derivative of the weight function with respect to r,
					// is lacking a factor of r.
					// The missing factor of r is recovered by
					// deltaE, which is missing a factor of 1/r
					wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
				}
				else {
					// Lucy Kernel, 2d
					wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
				}

				jmass = mass[jtype];
				D = alpha[itype][jtype]; // diffusion coefficient

				deltaE = 2.0 * imass * jmass / (imass + jmass);
				deltaE *= (rho[i] + rho[j]) / (rho[i] * rho[j]);
				deltaE *= D * (e[i] - e[j]) * wfd;

				de[i] += deltaE;
				if (newton_pair || j < nlocal) {
					de[j] -= deltaE;
				}

			}
		}
	}
}



/* ----------------------------------------------------------------------
   Setting for pair_style command
------------------------------------------------------------------------- */

void PairSPH_HEATCONDUCTION::set_style(int narg, char **arg)
{
	if (narg != 2)
		error->all(FLERR, "Illegal number of setting arguments for pair_style sph/heatconduction");
	nstep = atoi(arg[1]);
}

/* ----------------------------------------------------------------------
                              Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairSPH_HEATCONDUCTION::set_coeff(int narg, char **arg)
{
	if (narg != 4)
		error->all(FLERR, "Incorrect number of args for pair_style sph/heatconduction coefficients");
	if (!allocated)
		allocate();

	int ilo, ihi, jlo, jhi;
	force->bounds(arg[0], particle->ntypes, ilo, ihi);
	force->bounds(arg[1], particle->ntypes, jlo, jhi);

	double alpha_one = atof(arg[2]);
	double cut_one = atof(arg[3]);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			//printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
			cut[i][j] = cut_one;
			cutsq[i][j] = cut[i][j] * cut[i][j];
			alpha[i][j] = alpha_one;
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

void PairSPH_HEATCONDUCTION::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Not all pair sph/heatconduction coeffs are set");
  }
  force->type2pair[i][j] = pair_id;
  cut[j][i] = cut[i][j];
  cutsq[i][j] = cutsq[i][j];
  alpha[j][i] = alpha[i][j];
  setflag[j][i] = setflag[i][j];
  force->type2pair[j][i] = force->type2pair[i][j];

}

/* ---------------------------------------------------------------------- */

double PairSPH_HEATCONDUCTION::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;
  return 0.0;
}
