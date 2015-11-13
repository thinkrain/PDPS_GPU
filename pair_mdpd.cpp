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

#include "error.h"
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair_mdpd.h"
#include "particle.h"
#include "phy_const.h"
#include "random_mars.h"
#include "update.h"

using namespace PDPS_NS;
using namespace PhyConst;

#define DELTA 1
#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairMDPD::PairMDPD(PDPS *ps) : Pair(ps)
{
	random = NULL;
	a0 = b0 = gamma = sigma = NULL;
	cutd = cutdsq = NULL;
	rho_local = NULL;

	allocated = 0;
	inum_old = 0;

	newton_pair = 1;
}

/* ---------------------------------------------------------------------- */

PairMDPD::~PairMDPD()
{
	delete random;
	random = NULL;

	if (allocated) {
		memory->destroy(cut);
		memory->destroy(cutsq);
		memory->destroy(cutd);
		memory->destroy(cutdsq);
		memory->destroy(a0);
		memory->destroy(gamma);
		memory->destroy(sigma);
		memory->destroy(setflag);
	}
}

/* ----------------------------------------------------------------------
                  allocate
------------------------------------------------------------------------- */

void PairMDPD::allocate()
{
	allocated = 1;

    int n = particle->ntypes;
	
	cut = memory->create(cut, n+1, n+1, "PairMDPD: cut");
	cutsq = memory->create(cutsq, n+1, n+1, "PairMDPD: cutsq");
	cutd = memory->create(cutd, n+1, n+1, "PairMDPD: cutd");
	cutdsq = memory->create(cutdsq, n+1, n+1, "PairMDPD: cutdsq");
	a0 = memory->create(a0, n+1, n+1, "PairMDPD: a0");
	b0 = memory->create(b0, n+1, n+1, "PairMDPD: b0");
	gamma = memory->create(gamma, n+1, n+1, "PairMDPD: gamma");
	sigma = memory->create(sigma, n+1, n+1, "PairMDPD: sigma");
	
	setflag = memory->create(setflag, n+1, n+1, "PairMDPD: setflag");

    // initialize setflag for pair i and j
	for (int i = 0; i <= n; i++)
	for (int j = 0; j <= n; j++) {
		setflag[i][j] = 0;
	}

}

/* ----------------------------------------------------------------------
   Compute force for all paritcles
------------------------------------------------------------------------- */

void PairMDPD::compute(int eflag, int vflag)
{
	int i, j, ii, jj, inum, jnum, itype, jtype;
	double rijsq, rij, rij_inv, randnum, evdwl, fpair, dot, wc, wd, wrho; 
	double xtmp, ytmp, ztmp, delx, dely, delz;
	double vxtmp, vytmp, vztmp, delvx, delvy, delvz;

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	int *type = particle->type;
	int numlist;
	int nlocal = particle->nlocal;
	int nall = nlocal + particle->nghost;
	
	int *ilist;
	int *jlist;
	int **firstneigh;
	int *numneigh;
	
	double dt = update->dt;
	double dtinvsqrt = 1.0/sqrt(dt);   // inverse of time step

	ilist = neighbor->neighlist->ilist;
	inum = neighbor->neighlist->inum;
	firstneigh = neighbor->neighlist->firstneigh;
	numneigh = neighbor->neighlist->numneigh;

	evdwl = 0.0;

	if (eflag || vflag) ev_setup(eflag,vflag);

	// calculate local density first
	if (inum > inum_old) {
		inum_old = inum;
		memory->destroy(rho_local);
		memory->create(rho_local, inum_old, "pair_mdpd: rho_local");
	}
	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		itype = type[i];
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];

		jlist = firstneigh[i];
		jnum = numneigh[i];
		rho_local[i] = 0.0;
		for (jj = 0; jj < jnum; jj++) {
			j = jlist[jj];
			jtype = type[j];
			if (setflag[itype][jtype]) {
				delx = xtmp - x[j][0];
				dely = ytmp - x[j][1];
				delz = ztmp - x[j][2];
				rijsq = (delx*delx + dely*dely + delz*delz);
				if (rijsq < cutdsq[itype][jtype]) {
					rij = sqrt(rijsq);
					if (rijsq < EPSILON) continue;    // rij can be 0.0 in DPD systems
					rij_inv = 1.0/cutd[itype][jtype];
					wd = (1 - rij/cutd[itype][jtype]);
					wrho = 15.0/2.0/PI*rij_inv*rij_inv*rij_inv*wd*wd;
					rho_local[i] += wrho;
				}
			}
		}
	}

	// loop for all local owned particles 
	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		itype = type[i];
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];
		vxtmp = v[i][0];
		vytmp = v[i][1];
		vztmp = v[i][2];

		fpair = 0.0;
		jlist = firstneigh[i];
		jnum = numneigh[i];
		for (jj = 0; jj < jnum; jj++) {
			j = jlist[jj];
			jtype = type[j];
			if (setflag[itype][jtype]) {
				delx = xtmp - x[j][0];
				dely = ytmp - x[j][1];
				delz = ztmp - x[j][2];
				rijsq = (delx*delx + dely*dely + delz*delz);
				// withint the cutoff range
				if (rijsq < cutsq[itype][jtype]) {
					rij = sqrt(rijsq);
					if (rijsq < EPSILON) continue;    // rij can be 0.0 in DPD systems
					rij_inv = 1.0/rij;
					delvx = vxtmp - v[j][0];
					delvy = vytmp - v[j][1];
					delvz = vztmp - v[j][2];
					dot = delx*delvx + dely*delvy + delz*delvz;
					wc = 1.0 - rij/cut[itype][jtype];
					if (rijsq < cutdsq[itype][jtype]) wd = 1.0 - rij/cutd[itype][jtype];
					else wd = 0.0;
					randnum = random->gaussian();   // need to determined later

					// conservative force = a0 * wd
					// drag force = -gamma * wd^2 * (delx dot delv) / r
					// random force = sigma * wd * rnd * dtinvsqrt;
					fpair = a0[itype][jtype]*wc + b0[itype][jtype]*(rho_local[i] + rho_local[j])*wd;
					//fpair -= gamma[itype][jtype]*wd*wd*dot*rij_inv;
					//fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt;
					fpair -= gamma[itype][jtype]*pow(wc,rank)*dot*rij_inv;
					fpair += sigma[itype][jtype]*pow(wc,rank/2.0)*randnum*dtinvsqrt;
					//fpair += sigma[itype][jtype]*pow(wd,rank/2)*dtinvsqrt;
			
					f[i][0] += delx*rij_inv*fpair;
					f[i][1] += dely*rij_inv*fpair;
					f[i][2] += delz*rij_inv*fpair;
					
					// add pair force to the other particle
					
					f[j][0] -= delx*rij_inv*fpair;
					f[j][1] -= dely*rij_inv*fpair;
					f[j][2] -= delz*rij_inv*fpair;
					
					if (eflag) {
						evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd * wd;
					}

					fpair = fpair * rij_inv;
					if (evflag) ev_tally(i,j,nlocal,1,evdwl,0.0,fpair,delx,dely,delz);

				} // if (rij < cut[itype][jtype])
			}
		} // for (jj = 0; jj < nlocal; jj++)
	} // for (i = 0; i < nlocal; i++)

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairMDPD::init_one(int i, int j)
{
	//if (setflag[i][j] == 0) error->all(FLERR,"Not all pair coeffs are not set");

	sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);
    force->type2pair[i][j] = pair_id;

	cut[j][i] = cut[i][j];
	cutsq[j][i] = cutsq[i][j];
	cutd[j][i] = cutd[i][j];
	cutdsq[j][i] = cutdsq[i][j];
	a0[j][i] = a0[i][j];
	b0[j][i] = b0[i][j];
	gamma[j][i] = gamma[i][j];
	sigma[j][i] = sigma[i][j];
	setflag[j][i] = setflag[i][j];
	force->type2pair[j][i] = force->type2pair[i][j];
}

/* ----------------------------------------------------------------------
   Setting for pair_style command
------------------------------------------------------------------------- */

void PairMDPD::set_style(int narg, char **arg)
{
	if (narg != 5) error->all(FLERR,"Illegal PairMDPD command");

	int n = strlen(arg[0]) + 1;
	style = new char[n];
	strcpy(style,arg[0]);         // store pair style's name
	
	temperature = atof(arg[1]);   // store target temperature
	cut_global = atof(arg[2]);    // store the gloabl cutoff
	seed = atoi(arg[3]);          // store seed
	rank = atoi(arg[4]);

	if (seed <= 0) error->all(FLERR,"Illegal pair_style command");
	delete random;
	random = NULL;
	random = new RanMars(ps,seed);

	// reset cutoffs that have been explicitly set

	if (allocated) {
		int i, j;
		for (i=1; i <= particle->ntypes; i++) 
		for (j = i+1; j <= particle->ntypes; j++) {
			if (setflag[i][j]) cut[i][j] = cut_global;
		}
	}
}

/* ----------------------------------------------------------------------
                              Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairMDPD::set_coeff(int narg, char **arg)
{
	int ilo,ihi,jlo,jhi;
    
	if (narg != 7) error->all(FLERR, "Incorrect args for pair_coeff of pair_style dpd"); 
	if(!allocated) allocate();

    // find the lower and upper bound of type set in the pair_coeff
    force->bounds(arg[0],particle->ntypes,ilo,ihi);
    force->bounds(arg[1],particle->ntypes,jlo,jhi);

	double cut_one = cut_global;
	cut_one = atof(arg[5]);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo,i); j <= jhi; j++) {
			a0[i][j] = atof(arg[2]);
			b0[j][i] = atof(arg[3]);
			gamma[i][j] = atof(arg[4]);
			cut[i][j] = cut_one;
			cutsq[i][j] = cut_one * cut_one;
			cutd[i][j] = atof(arg[6]);
			cutdsq[i][j] = cutd[i][j] * cutd[i][j];
			setflag[i][j] = 1;        // "1" means this pair has been setup
			count++;
		}
	}

	if (count == 0) {
		error->all(FLERR,"Incorrest args for pair coefficients");
	}
}
