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
#include "pair_sph_lj.h"
#include "particle.h"
#include "random_mars.h"
#include "update.h"
#include "group.h"

using namespace PDPS_NS;


#define DELTA 1
#define EPSILON 1.0e-10
#define PI 3.1416

/* ---------------------------------------------------------------------- */

PairSPH_LJ::PairSPH_LJ(PDPS *ps) : Pair(ps)
{

	 first = 1;
	 newton_pair = 1;
	 allocated = 0;
	 cut = NULL;
	 cutsq = NULL;
}

/* ---------------------------------------------------------------------- */

PairSPH_LJ::~PairSPH_LJ()
{
	if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(viscosity);
  }

}

/* ---------------------------------------------------------------------- */

void PairSPH_LJ::allocate()
{
	 allocated = 1;
  int n = particle->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, n + 1, n + 1, "pair:viscosity");

}

/* ----------------------------------------------------------------------
   Compute force for all paritcles
------------------------------------------------------------------------- */

void PairSPH_LJ::compute(int eflag, int vflag)
{
	int i, j, ii, jj, inum, jnum, itype, jtype;
	double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

	int *ilist, *jlist, *numneigh, **firstneigh;
	double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, q;
	double rsq, rij, rij_inv, wfd, delVdotDelR, mu, deltaE, ci, cj, lrc;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
//  else
//    evflag = vflag_fdotr = 0;
 
  double **v = particle->vest;
  double **x = particle->x;
  double **f = particle->f;
  double *rho = particle->rho;
  double *mass = particle->mass;
  double *e = particle->e;
  double *de = particle->de;
  double *drho = particle->drho;
  int *type = particle->type;
  int *mask = particle->mask;
  int nlocal = particle->nlocal;
  double wf;
//  int newton_pair = force->newton_pair;
  // check consistency of pair coefficients
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

  // loop over neighbors of my particles
  for (ii = 0; ii < inum; ii++) {
	  i = ilist[ii];
	  xtmp = x[i][0];
	  ytmp = x[i][1];
	  ztmp = x[i][2];
	  vxtmp = v[i][0];
	  vytmp = v[i][1];
	  vztmp = v[i][2];
	  itype = type[i];
	  jlist = firstneigh[i];
	  jnum = numneigh[i];

	  imass = mass[itype];

	  // compute pressure of particle i with LJ EOS
	//  LJEOS2(rho[i], e[i], cv, &fi, &ci);
	//  fi /= (rho[i] * rho[i]);
	  //printf("fi = %f\n", fi);

	  for (jj = 0; jj < jnum; jj++) {
		  j = jlist[jj];

		  delx = xtmp - x[j][0];
		  dely = ytmp - x[j][1];
		  delz = ztmp - x[j][2];
		  rsq = delx * delx + dely * dely + delz * delz;
		  jtype = type[j];
		  jmass = mass[jtype];

		  if (rsq < cutsq[itype][jtype]) {

			  q = cut[itype][jtype] / sqrt(rsq);
			  rij_inv = 1.0 / sqrt(rsq);

			  // total pair force & thermal energy increment
			  fpair = K * (pow(q, 12) - pow(q, 4)) * rij_inv;

			  f[i][0] += delx * fpair * rij_inv;
			  f[i][1] += dely * fpair * rij_inv;
			  f[i][2] += delz * fpair * rij_inv;

			  if (newton_pair || j < nlocal) {
				  f[j][0] -= delx * fpair * rij_inv;
				  f[j][1] -= dely * fpair * rij_inv;
				  f[j][2] -= delz * fpair * rij_inv;

			  }

			  if (evflag)
				  ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
		  }
	  }
  }

}



/* ----------------------------------------------------------------------
   Setting for pair_style command
------------------------------------------------------------------------- */

void PairSPH_LJ::set_style(int narg, char **arg)
{
	//	if (narg != 4)
	//		error->all(FLERR, "Illegal number of setting arguments for pair_style sph/idealgas");
}

/* ----------------------------------------------------------------------
                              Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairSPH_LJ::set_coeff(int narg, char **arg)
{
	if (narg != 4)
		error->all(FLERR, "Incorrect args for pair_style sph/lj coefficients");
	if (!allocated)
		allocate();

	int ilo, ihi, jlo, jhi;
	force->bounds(arg[0], particle->ntypes, ilo, ihi);
	force->bounds(arg[1], particle->ntypes, jlo, jhi);

	K = atof(arg[2]);
	double cut_one = atof(arg[3]);
	

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			cut[i][j] = cut_one;
			cutsq[i][j] = cut[i][j] * cut[i][j];
			setflag[i][j] = 1;

			//cut[j][i] = cut[i][j];
			//viscosity[j][i] = viscosity[i][j];
			//setflag[j][i] = 1;
			count++;
		}
	}

	if (count == 0)
		error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairSPH_LJ::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Not all pair sph/lj coeffs are set");
  }
  force->type2pair[i][j] = pair_id;

  cut[j][i] = cut[i][j];
  cutsq[j][i] = cutsq[i][j];
  setflag[j][i] = setflag[i][j];
  force->type2pair[j][i] = force->type2pair[i][j];

}

