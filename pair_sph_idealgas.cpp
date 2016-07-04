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
#include "pair_sph_idealgas.h"
#include "particle.h"
#include "random_mars.h"
#include "update.h"

using namespace PDPS_NS;


#define DELTA 1
#define EPSILON 1.0e-10
#define PI 3.1416

/* ---------------------------------------------------------------------- */

PairSPH_IDEALGAS::PairSPH_IDEALGAS(PDPS *ps) : Pair(ps)
{
	cutd = cutdsq = NULL;
	 newton_pair = 1;
	 allocated = 0;
	 first = 1;
	 cubic_flag = 0;
	 quintic_flag = 0;
	 h = 0.0;
	 cut = NULL;
	 cutsq = NULL;
}

/* ---------------------------------------------------------------------- */

PairSPH_IDEALGAS::~PairSPH_IDEALGAS()
{
	if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(viscosity);
  }

}

/* ---------------------------------------------------------------------- */

void PairSPH_IDEALGAS::allocate()
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

void PairSPH_IDEALGAS::compute(int eflag, int vflag)
{
	int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, ih, q, ihsq;
  double rsq, tmp, wfd, delVdotDelR, mu, deltaE;
  double ci, cj;

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

  if (nstep != 0) {
	  if ((update->ntimestep % nstep) == 0) {

		  // initialize density with self-contribution,
		  for (i = 0; i < nlocal; i++) {
			  itype = type[i];
			  imass = mass[itype];  

			  if (domain->dim == 3) {

				  // Cubic spline kernel, 3d
				  wf = a3D;
			  }
			  else {

				  // Cubic spline kernel, 2d
				  wf = a2D;
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
					 q = sqrt(rsq) / h;

					 if (cubic_flag == 1){
						 if (q < 1)
							 wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
						 else
							 wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
					 }
					 else if (quintic_flag == 1)
						 wf = (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (2 * q + 1);

					 if (domain->dim == 3)
						 wf = wf * a3D;
					 else
						 wf = wf * a2D;

					  rho[i] += mass[jtype] * wf;
					  rho[j] += mass[itype] * wf;
				  }

			  }
   		  }
		

	  }
  }
   
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

	  fi = 0.4 * e[i] / imass / rho[i]; // ideal gas EOS; this expression is fi = pressure / rho^2
	  ci = sqrt(0.4*e[i] / imass); // speed of sound with heat capacity ratio gamma=1.4

	  for (jj = 0; jj < jnum; jj++) {
		  j = jlist[jj];
//		  j &= NEIGHMASK;

		  delx = xtmp - x[j][0];
		  dely = ytmp - x[j][1];
		  delz = ztmp - x[j][2];
		  rsq = delx * delx + dely * dely + delz * delz;
		  jtype = type[j];
		  jmass = mass[jtype];

		  if (rsq < cutsq[itype][jtype]) {
			  q = sqrt(rsq) / h;

			  if (cubic_flag == 1){
				  if (q < 1)
					  wfd = - 3 * q + 2.25 * q * q;
				  else
					  wfd = -0.75 * (2 - q) * (2 - q);
			  }
			  else if (quintic_flag == 1)
				  wfd = -2 * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (2 * q + 1) + 2 * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0);

			  if (domain->dim == 3)
				  wfd = wfd * a3D / h;
			  else
				  wfd = wfd * a2D / h;

			  fj = 0.4 * e[j] / jmass / rho[j];

			  // dot product of velocity delta and distance vector
			  delVdotDelR = delx * (vxtmp - v[j][0]) + dely * (vytmp - v[j][1])
				  + delz * (vztmp - v[j][2]);

			  // artificial viscosity (Monaghan 1992)
			  if (delVdotDelR < 0.) {
				  cj = sqrt(0.4*e[j] / jmass);
				  mu = h * delVdotDelR / (rsq + 0.01 * h * h);
				  fvisc = -viscosity[itype][jtype] * (ci + cj) * mu / (rho[i] + rho[j]);
			  }
			  else {
				  fvisc = 0.;
			  }

			  // total pair force & thermal energy increment
			  fpair = -imass * jmass * (fi + fj + fvisc) * wfd;
			  deltaE = -0.5 * fpair * delVdotDelR;

			  f[i][0] += delx * fpair;
			  f[i][1] += dely * fpair;
			  f[i][2] += delz * fpair;

			  // and change in density
			  drho[i] += jmass * delVdotDelR * wfd;

			  // change in thermal energy
			  de[i] += deltaE;

			  if (newton_pair || j < nlocal) {
				  f[j][0] -= delx * fpair;
				  f[j][1] -= dely * fpair;
				  f[j][2] -= delz * fpair;
				  de[j] += deltaE;
				  drho[j] += imass * delVdotDelR * wfd;
			  }

			  if (evflag)
				  ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);

			  if (evflag)
				  ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely,
				  delz);
		  }
	  }
  }




//  if (vflag_fdotr) virial_fdotr_compute();
}



/* ----------------------------------------------------------------------
   Setting for pair_style command
------------------------------------------------------------------------- */

void PairSPH_IDEALGAS::set_style(int narg, char **arg)
{
	if (narg != 3)
		error->all(FLERR, "Illegal number of setting arguments for pair_style sph/idealgas");
	if (strcmp(arg[1], "Cubic") == 0)
		cubic_flag = 1;
	else if (strcmp(arg[1], "Quintic") == 0)
		quintic_flag = 1;
	else
		error->all(FLERR, "Wrong Kernel function");
	nstep = atoi(arg[2]);
}

/* ----------------------------------------------------------------------
                              Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairSPH_IDEALGAS::set_coeff(int narg, char **arg)
{
	if (narg != 4)
		error->all(FLERR, "Incorrect number of args for pair_style sph/idealgas coefficients");
	if (!allocated)
		allocate();

	int ilo, ihi, jlo, jhi;
	force->bounds(arg[0], particle->ntypes, ilo, ihi);
	force->bounds(arg[1], particle->ntypes, jlo, jhi);

	double viscosity_one = atof(arg[2]);
	double cut_one = atof(arg[3]);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			viscosity[i][j] = viscosity_one;
			//printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
			cut[i][j] = cut_one * 2;
			cutsq[i][j] = cut[i][j] * cut[i][j];
			setflag[i][j] = 1;
			count++;
		}
	}

	h = cut_one;

	if (cubic_flag == 1){
		a2D = 10.0 / 7.0 / PI / h / h;
		a3D = 1.0 / PI / h / h / h;
	}
	else if (quintic_flag == 1){
		a2D = 7.0 / 4.0 / PI / h / h;
		a3D = 21.0 / 16.0 / PI / h / h / h;
	}

	if (count == 0)
		error->all(FLERR, "Incorrect args for pair sph/idealgas coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairSPH_IDEALGAS::init_one(int i, int j) {

	if (setflag[i][j] == 0) {
		error->all(FLERR, "All pair sph/idealgas coeffs are not set");
	}

	cut[j][i] = cut[i][j];
	cutsq[j][i] = cutsq[j][i];

}

/* ---------------------------------------------------------------------- */

double PairSPH_IDEALGAS::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;
  return 0.0;
}
