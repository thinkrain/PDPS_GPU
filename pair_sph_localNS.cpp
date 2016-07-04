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
#include "pair_sph_localNS.h"
#include "particle.h"
#include "random_mars.h"
#include "update.h"

using namespace PDPS_NS;


#define DELTA 1
#define EPSILON 1.0e-10
#define PI 3.1416

/* ---------------------------------------------------------------------- */

PairSPH_LOCALNS::PairSPH_LOCALNS(PDPS *ps) : Pair(ps)
{
	random = NULL;
	cutd = cutdsq = NULL;
	 first = 1;
	 newton_pair = 1;
	 local_kernel = 0;
	 couple_flag = 0;
	 allocated = 0;
	 sigma = 1.5;
	 cut = NULL;
	 cutsq = NULL;
}

/* ---------------------------------------------------------------------- */

PairSPH_LOCALNS::~PairSPH_LOCALNS()
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

void PairSPH_LOCALNS::allocate()
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

void PairSPH_LOCALNS::compute(int eflag, int vflag)
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
  // initial local kernel length
  if (local_kernel == 1){
	  for (i = 0; i < nlocal; i++){
		  itype = type[i];
		  hlocal[i] = cut[itype][itype];
	  }
		  
  }
  // loop over neighbors of my particles
  if (nstep != 0) {
	  if ((update->ntimestep % nstep) == 0) {

		  // initialize density with self-contribution,
		  for (ii = 0; ii < inum; ii++) {
			  i = ilist[ii];
			  itype = type[i];
			  if (itype == phase_f)
			  {
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
				  poro[i] = 1.0;
			  }
			  
		  }

		  // add density at each particle via kernel function overlap
		  for (ii = 0; ii < inum; ii++) {
			  i = ilist[ii];
			  itype = type[i];
			  if (itype == phase_f){
				  xtmp = x[i][0];
				  ytmp = x[i][1];
				  ztmp = x[i][2];
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
						  }		//consider porosity of liquid
						  if (itype == jtype)
							  rho[i] += mass[jtype] * wf;
						  else
							  poro[i] -= volume[j] * wf;
					  }

				  }
				  if (local_kernel == 1){
					  if (domain->dim == 3)
						  hlocal[i] = sigma * pow((mass[itype] / rho[i]), 1.0 / 3.0);
					  else
						  hlocal[i] = sigma * pow((mass[itype] / rho[i]), 1.0 / 2.0);
				  }

			  }		//  if setflag[itype][itype]
			  
   		  }		// inum
		

	  }		// update->ntimestep % nstep) == 0
  }		//	nstep != 0

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];  
//	if (itype == phase_f){
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];
		vxtmp = v[i][0];
		vytmp = v[i][1];
		vztmp = v[i][2];
		jlist = firstneigh[i];
		jnum = numneigh[i];
		imass = mass[itype];

		// compute pressure of particle i with Tait EOS
		tmp = rho[i] / rho0[itype];
		fi = tmp * tmp * tmp;
		fi = B[itype] * (fi * fi * tmp - 1.0) / (rho[i] * rho[i]);

		for (jj = 0; jj < jnum; jj++) {
			j = jlist[jj];
			//     j &= NEIGHMASK;
			jtype = type[j];

			if (setflag[itype][jtype]){
				delx = xtmp - x[j][0];
				dely = ytmp - x[j][1];
				delz = ztmp - x[j][2];
				rsq = delx * delx + dely * dely + delz * delz;

				jmass = mass[jtype];

				if (rsq < cutsq[itype][jtype]) {

					h = cut[itype][jtype];
					ih = 1.0 / h;
					ihsq = ih * ih;
					wfd = h - sqrt(rsq);
					if (domain->dim == 3) {
						// Lucy Kernel, 3d
						// Note that wfd, the derivative of the weight function with respect to r,
						// is lacking a factor of r.
						// The missing factor of r is recovered by
						// (1) using delV . delX instead of delV . (delX/r) and
						// (2) using f[i][0] += delx * fpair instead of f[i][0] += (delx/r) * fpair
						wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
					}
					else {
						// Lucy Kernel, 2d
						wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
					}

					// compute pressure  of particle j with Tait EOS
					tmp = rho[j] / rho0[jtype] / poro[j];
					fj = tmp * tmp * tmp;
					fj = B[jtype] * (fj * fj * tmp - 1.0) / (rho[j] * rho[j]);
					// dot product of velocity delta and distance vector
					delVdotDelR = delx * (vxtmp - v[j][0]) + dely * (vytmp - v[j][1])
						+ delz * (vztmp - v[j][2]);

					// artificial viscosity (Monaghan 1992)
					if (delVdotDelR < 0.) {
						mu = h * delVdotDelR / (rsq + 0.01 * h * h);
						fvisc = -viscosity[itype][jtype] * (soundspeed[itype]
							+ soundspeed[jtype]) * mu / (rho[i] + rho[j]);
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
				}		//		rsq < cutsq[itype][jtype]
			}		//  setflag[itype][jtype]
		}		//  jnum
//	}	//	setflag[itype][itype]
	


     
  }		//  inum 

//	calcuate the coupling force on other phase  
  if (couple_flag == 1)
  {
	  double weight;
	  double fx, fy, fz;
	  for (ii = 0; ii < inum; ii++) {
		  i = ilist[ii];

		  itype = type[i];
		  if (!setflag[itype][itype]){
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
				  if (setflag[jtype][jtype]){
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


  }
	

//  if (vflag_fdotr) virial_fdotr_compute();
}



/* ----------------------------------------------------------------------
   Setting for pair_style command
------------------------------------------------------------------------- */

void PairSPH_LOCALNS::set_style(int narg, char **arg)
{
	if (narg != 4)
		error->all(FLERR, "Illegal number of setting arguments for pair_style sph/taitwater");
	nstep = atoi(arg[1]);
	if (!strcmp(arg[2], "coupled")) {
		couple_flag = 1;
	}

	phase_f = atoi(arg[3]);
}

/* ----------------------------------------------------------------------
                              Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairSPH_LOCALNS::set_coeff(int narg, char **arg)
{
   if (narg != 6)
    error->all(FLERR, "Incorrect args for pair_style sph/taitwater coefficients");
  if (!allocated)
    allocate();
 // if (atoi(arg[0]) == atoi(arg[1])){
	  int ilo, ihi, jlo, jhi;
	  force->bounds(arg[0], particle->ntypes, ilo, ihi);
	  force->bounds(arg[1], particle->ntypes, jlo, jhi);

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
		  error->all(FLERR, "Incorrect args for pair coefficients");
//  }
//  else{
//	  int i = atoi(arg[0]);
//	  int j = atoi(arg[1]);
//	  setflag[i][j] = 1;
//	  cut[i][j] = atof(arg[5]);
//	  cutsq[i][j] = cut[i][j] * cut[i][j];
//  }

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairSPH_LOCALNS::init_one(int i, int j) {

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

double PairSPH_LOCALNS::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;
  return 0.0;
}
