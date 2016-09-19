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
#include "pair_sph_couple.h"
#include "particle.h"
#include "random_mars.h"
#include "update.h"

using namespace PDPS_NS;


#define DELTA 1
#define EPSILON 1.0e-10
#define PI 3.1416

/* ---------------------------------------------------------------------- */

PairSPH_COUPLE::PairSPH_COUPLE(PDPS *ps) : Pair(ps)
{
	
	 first = 1;
	 newton_pair = 1;
	 local_kernel = 0;
	 allocated = 0;
	 sigma = 1.5;
	 cut = NULL;
	 cutsq = NULL;
	 poro_flag = 0;
	 poroset_flag = 0;
}

/* ---------------------------------------------------------------------- */

PairSPH_COUPLE::~PairSPH_COUPLE()
{
	if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(rho0);
    memory->destroy(soundspeed);
    memory->destroy(B);
    memory->destroy(viscosity);
  }

}

/* ---------------------------------------------------------------------- */

void PairSPH_COUPLE::allocate()
{
	 allocated = 1;
  int n = particle->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(rho0, n + 1, "pair:rho0");
  memory->create(soundspeed, n + 1, "pair:soundspeed");
  memory->create(B, n + 1, "pair:B");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, n + 1, n + 1, "pair:viscosity");

}

/* ----------------------------------------------------------------------
   Compute force for all paritcles
------------------------------------------------------------------------- */

void PairSPH_COUPLE::compute(int eflag, int vflag)
{
	int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc;
  double rsq, tmp, wfd, delVdotDelR, mu, deltaE;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
//  else
//    evflag = vflag_fdotr = 0;
 
  double **v = particle->vest;
  double **x = particle->x;
  double **f = particle->f;
  double *rho = particle->rho;
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


//	calcuate the coupling force on other phase
	double weight;
	double fx, fy, fz;
	for (ii = 0; ii < nlocal; ii++) {
		i = ilist[ii];

		itype = type[i];
		if (itype == phase_s){
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
				if (jtype == phase_f){
					jmass = mass[jtype];
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
				f[i][0] = fx / weight * volume[i] * rho0[phase_f] / jmass;
				f[i][1] = fy / weight * volume[i] * rho0[phase_f] / jmass;
				f[i][2] = fz / weight * volume[i] * rho0[phase_f] / jmass;

				if (poroset_flag == 1){
					f[i][0] = f[i][0] / poroset;
					f[i][1] = f[i][1] / poroset;
					f[i][2] = f[i][2] / poroset;
				}
				if (poro_flag == 1){
					f[i][0] = f[i][0] / (1 - poro[i]);
					f[i][1] = f[i][1] / (1 - poro[i]);
					f[i][2] = f[i][2] / (1 - poro[i]);

				}
				

				for (jj = 0; jj < jnum; jj++) {
					j = jlist[jj];
					//     j &= NEIGHMASK;
					jtype = type[j];
					jmass = mass[jtype];
					if (jtype == phase_f){
						delx = xtmp - x[j][0];
						dely = ytmp - x[j][1];
						delz = ztmp - x[j][2];
						rsq = delx * delx + dely * dely + delz * delz;
						if (rsq < radius[i]) {
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

		} // setflag[i][i]
	} // inum  
//  if vflag_fdotr) virial_fdotr_compute();
}



/* ----------------------------------------------------------------------
   Setting for pair_style command
------------------------------------------------------------------------- */

void PairSPH_COUPLE::set_style(int narg, char **arg)
{
	if (strcmp(arg[1], "Cubic") == 0)
		cubic_flag = 1;
	else if (strcmp(arg[1], "Quintic") == 0)
		quintic_flag = 1;
	else
		error->all(FLERR, "Wrong Kernel function");
	if (strcmp(arg[2], "poro") == 0)
		poro_flag = 1;
	else if (strcmp(arg[2], "poroset") == 0){
		poroset_flag = 1;
		poroset = atof(arg[3]);
	}


}

/* ----------------------------------------------------------------------
                              Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairSPH_COUPLE::set_coeff(int narg, char **arg)
{


  if (narg != 6)
	  error->all(FLERR, "Incorrect args for pair_style sph/taitwater coefficients");
  if (!allocated)
	  allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(arg[0], particle->ntypes, ilo, ihi);
  force->bounds(arg[1], particle->ntypes, jlo, jhi);

  // record fluid phase and solid phase
  phase_f = atoi(arg[0]);
  phase_s = atoi(arg[1]);

  double rho0_one = atof(arg[2]);
  double soundspeed_one = atof(arg[3]);
  double viscosity_one = atof(arg[4]);
  double cut_one = atof(arg[5]);
  double B_one = soundspeed_one * soundspeed_one * rho0_one / 7.0;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
	  rho0[i] = rho0_one;
	  soundspeed[i] = soundspeed_one;
	  B[i] = B_one;
	  for (int j = MAX(jlo, i); j <= jhi; j++) {
		  rho0[j] = rho0_one;
		  soundspeed[j] = soundspeed_one;
		  B[j] = B_one;
		  viscosity[i][j] = viscosity_one;
		  cut[i][j] = 2 * cut_one;
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
	  error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairSPH_COUPLE::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Not all pair sph/taitwater coeffs are set");
  }
  force->type2pair[i][j] = pair_id;

  cut[j][i] = cut[i][j];
  cutsq[j][i] = cutsq[i][j];
  setflag[j][i] = setflag[i][j];
  viscosity[j][i] = viscosity[i][j];
  force->type2pair[j][i] = force->type2pair[i][j];

}

/* ---------------------------------------------------------------------- */

double PairSPH_COUPLE::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;
  return 0.0;
}
