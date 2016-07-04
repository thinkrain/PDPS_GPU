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
	random = NULL;
	cutd = cutdsq = NULL;
	 first = 1;
	 newton_pair = 1;
	 local_kernel = 0;
	 allocated = 0;
	 sigma = 1.5;
	 cut = NULL;
	 cutsq = NULL;
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
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, h, ih, ihsq;
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
  double *volume = particle->volume;
  double *hlocal = particle->hlocal;
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


//	calcuate the coupling force on other phase
	double weight;
	double fx, fy, fz;
	for (ii = 0; ii < inum; ii++) {
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
					delx = xtmp - x[j][0];
					dely = ytmp - x[j][1];
					delz = ztmp - x[j][2];
					rsq = delx * delx + dely * dely + delz * delz;
					if (rsq < cutsq[itype][jtype]) {
						h = cut[itype][jtype];
						ih = 1.0 / h;
						ihsq = ih * ih;
						jmass = mass[jtype];
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
			}

		} // setflag[i][i]
	} // inum  
//  if (vflag_fdotr) virial_fdotr_compute();
}



/* ----------------------------------------------------------------------
   Setting for pair_style command
------------------------------------------------------------------------- */

void PairSPH_COUPLE::set_style(int narg, char **arg)
{
	if (narg != 2)
		error->all(FLERR, "Illegal number of setting arguments for pair_style sph/taitwater");
	nstep = atoi(arg[1]);

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

  // coupling force between different types of particles
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
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      viscosity[i][j] = viscosity_one;
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
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
    error->all(FLERR,"Incorrect args for pair coefficients");
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
