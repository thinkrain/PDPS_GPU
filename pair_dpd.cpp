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
#include "pair_dpd.h"
#include "particle.h"
#include "random_mars.h"
#include "update.h"

using namespace PDPS_NS;

#define DELTA 1
#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairDPD::PairDPD(PDPS *ps) : Pair(ps)
{
	random = NULL;
	a0 = gamma = sigma = NULL;
	allocated = 0;

	newton_pair = 1;
}

/* ---------------------------------------------------------------------- */

PairDPD::~PairDPD()
{
	delete random;
	random = NULL;

	if (allocated) {
		memory->destroy(cut);
		memory->destroy(cutsq);
		memory->destroy(a0);
		memory->destroy(gamma);
		memory->destroy(sigma);
		memory->destroy(setflag);
	}
}

/* ---------------------------------------------------------------------- */

void PairDPD::allocate()
{
	allocated = 1;

    int n = particle->ntypes;
	
	cut = memory->create(cut, n+1, n+1, "PairDPD: cut");
	cutsq = memory->create(cutsq, n+1, n+1,"PairDPD: cutsq");
	a0 = memory->create(a0, n+1, n+1, "PairDPD: a0");
	gamma = memory->create(gamma, n+1, n+1,"PairDPD: gamma");
	sigma = memory->create(sigma, n+1, n+1, "PairDPD: sigma");
	setflag = memory->create(setflag, n+1, n+1, "PairDPD: setflag");

    // initialize setflag for pair i and j
	for (int i = 0; i <= n; i++)
	for (int j = 0; j <= n; j++) {
		setflag[i][j] = 0;
	}

}

/* ----------------------------------------------------------------------
   Compute force for all paritcles
------------------------------------------------------------------------- */

void PairDPD::compute(int eflag, int vflag)
{
	int i, j, ii, jj, inum, jnum, itype, jtype;
	double rijsq, rij, rij_inv, randnum, evdwl, fpair, dot, wd; 
	double xtmp, ytmp, ztmp, delx, dely, delz;
	double vxtmp, vytmp, vztmp, delvx, delvy, delvz;

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	int *type = particle->type;
	int radius_flag = particle->radius_flag;
	double *radius = particle->radius;
	int numlist;
	int nlocal = particle->nlocal;
	int nall = nlocal + particle->nghost;
	
	int *ilist;
	int *jlist;
	//int **list;
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

	// loop for all local owned particles 
	double radius_cut;
	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		itype = particle->type[i];
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
				rij = sqrt(rijsq);
				// withint the cutoff range
				if (rijsq < cutsq[itype][jtype]) {
					if (rij < EPSILON) continue;    // rij can be 0.0 in DPD systems

					if (radius_flag == 1) {
						radius_cut = radius[i] + radius[j];
						if (rij >= radius_cut) continue;
					}

					rij_inv = 1.0/rij;
					delvx = vxtmp - v[j][0];
					delvy = vytmp - v[j][1];
					delvz = vztmp - v[j][2];
					dot = delx*delvx + dely*delvy + delz*delvz;
					if (radius_flag) wd = 1.0 - rij / radius_cut;
					else wd = 1.0 - rij/cut[itype][jtype];
					randnum = random->gaussian();   // need to determined later

					// conservative force = a0 * wd
					// drag force = -gamma * wd^2 * (delx dot delv) / r
					// random force = sigma * wd * rnd * dtinvsqrt;
					fpair = a0[itype][jtype]*wd;
					//fpair -= gamma[itype][jtype]*wd*wd*dot*rij_inv;
					//fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt;
					fpair -= gamma[itype][jtype]*pow(wd,rank)*dot*rij_inv;
					fpair += sigma[itype][jtype]*pow(wd,rank/2.0)*randnum*dtinvsqrt;
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

void PairDPD::init_one(int i, int j)
{
	//if (setflag[i][j] == 0) error->all(FLERR,"Not all pair coeffs are not set");

	sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);
    force->type2pair[i][j] = pair_id;

	cut[j][i] = cut[i][j];
	cutsq[j][i] = cutsq[i][j];
	a0[j][i] = a0[i][j];
	gamma[j][i] = gamma[i][j];
	sigma[j][i] = sigma[i][j];
	setflag[j][i] = setflag[i][j];
	force->type2pair[j][i] = force->type2pair[i][j];
}

/* ----------------------------------------------------------------------
   Setting for pair_style command
------------------------------------------------------------------------- */

void PairDPD::set_style(int narg, char **arg)
{
	if (narg != 5) error->all(FLERR,"Illegal PairDPD command");

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

void PairDPD::set_coeff(int narg, char **arg)
{
	int ilo,ihi,jlo,jhi;
    
	if (narg != 5) error->all(FLERR, "Incorrect args for pair_coeff of pair_style dpd"); 
	if(!allocated) allocate();

    // find the lower and upper bound of type set in the pair_coeff
    force->bounds(arg[0],particle->ntypes,ilo,ihi);
    force->bounds(arg[1],particle->ntypes,jlo,jhi);

	double cut_one = cut_global;
	if (narg == 5) cut_one = atof(arg[4]);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo,i); j <= jhi; j++) {
			a0[i][j] = atof(arg[2]);
			gamma[i][j] = atof(arg[3]);
			cut[i][j] = cut_one;
			cutsq[i][j] = cut_one * cut_one;
			setflag[i][j] = 1;        // "1" means this pair has been setup
			count++;
		}
	}

	if (count == 0) {
		error->all(FLERR,"Incorrest args for pair coefficients");
	}
}
