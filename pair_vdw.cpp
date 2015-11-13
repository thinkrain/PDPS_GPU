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
#include "pair_vdw.h"
#include "particle.h"
#include "random_mars.h"
#include "update.h"

using namespace PDPS_NS;

#define DELTA 1
#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairVDW::PairVDW(PDPS *ps) : Pair(ps)
{
	Ha = NULL;
	allocated = 0;

	newton_pair = 1;
}

/* ---------------------------------------------------------------------- */

PairVDW::~PairVDW()
{
	if (allocated) {
		memory->destroy(cut);
		memory->destroy(cutsq);
		memory->destroy(Ha);
		memory->destroy(setflag);
	}
}

/* ----------------------------------------------------------------------
allocate
------------------------------------------------------------------------- */

void PairVDW::allocate()
{
	allocated = 1;

	int n = particle->ntypes;

	cut = memory->create(cut, n + 1, n + 1, "PairVDW: cut");
	cutsq = memory->create(cutsq, n + 1, n + 1, "PairVDW: cutsq");
	Ha = memory->create(Ha, n + 1, n + 1, "PairVDW: a0");
	setflag = memory->create(setflag, n + 1, n + 1, "PairVDW: setflag");

	// initialize setflag for pair i and j
	for (int i = 0; i <= n; i++)
	for (int j = 0; j <= n; j++) {
		setflag[i][j] = 0;
	}
}

/* ----------------------------------------------------------------------
   Compute force for all paritcles
------------------------------------------------------------------------- */

void PairVDW::compute(int eflag, int vflag)
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
	double dtinvsqrt = 1.0 / sqrt(dt);   // inverse of time step

	ilist = neighbor->neighlist->ilist;
	inum = neighbor->neighlist->inum;
	firstneigh = neighbor->neighlist->firstneigh;
	numneigh = neighbor->neighlist->numneigh;

	evdwl = 0.0;

	if (eflag || vflag) ev_setup(eflag, vflag);

	// loop for all local owned particles 
	double radius_cut;
	double sij, sijsq, Reff;
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
				if (radius_flag == 1) sij = rij - (radius[i] + radius[j]);
				else sij = rij - cut[itype][jtype];
				
				if (sij < cut[itype][jtype]) {			
					if (radius_flag == 1) {
						if (sij >= 0.5*MAX(radius[i],radius[j])) continue;
						Reff = 1.0 / (1.0/radius[i] + 1.0/radius[j]);
					}
					else Reff = cut[itype][jtype] / 4.0;
					
					if (sij < delta_min) {
						sij = delta_min;
					}
					sijsq = sij * sij;

					rij_inv = 1.0 / rij;

					fpair = -Ha[itype][jtype]*Reff/(6*sijsq) * (1 - 1.0/(1.0 + lambda/(b*sij)));

					f[i][0] += delx*rij_inv*fpair;
					f[i][1] += dely*rij_inv*fpair;
					f[i][2] += delz*rij_inv*fpair;

					// add pair force to the other particle

					f[j][0] -= delx*rij_inv*fpair;
					f[j][1] -= dely*rij_inv*fpair;
					f[j][2] -= delz*rij_inv*fpair;

				} // if (rij < cut[itype][jtype])
			}
		} // for (jj = 0; jj < nlocal; jj++)
	} // for (i = 0; i < nlocal; i++)

}

/* ----------------------------------------------------------------------
init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

void PairVDW::init_one(int i, int j)
{
	force->type2pair[i][j] = pair_id;

	cut[j][i] = cut[i][j];
	cutsq[j][i] = cutsq[i][j];
	Ha[j][i] = Ha[i][j];
	setflag[j][i] = setflag[i][j];
	force->type2pair[j][i] = force->type2pair[i][j];
}

/* ----------------------------------------------------------------------
   Setting for pair_style command
------------------------------------------------------------------------- */

void PairVDW::set_style(int narg, char **arg)
{
	if (narg != 5) error->all(FLERR, "Illegal PairVDW command");

	int n = strlen(arg[0]) + 1;
	style = new char[n];
	strcpy(style, arg[0]);         // store pair style's name

	b = atof(arg[1]);
	lambda = atof(arg[2]);
	delta_min = atof(arg[3]);      // minimum separation distance to avoid sigularity
	cut_global = atof(arg[4]);     // store the gloabl cutoff

	if (cut_global <= 0) error->all(FLERR, "Illegal pair_style command");

	// reset cutoffs that have been explicitly set

	if (allocated) {
		int i, j;
		for (i = 1; i <= particle->ntypes; i++)
		for (j = i + 1; j <= particle->ntypes; j++) {
			if (setflag[i][j]) cut[i][j] = cut_global;
		}
	}
}

/* ----------------------------------------------------------------------
   Set Coeff for pair_coeff command
------------------------------------------------------------------------- */

void PairVDW::set_coeff(int narg, char **arg)
{
	int ilo, ihi, jlo, jhi;

	if (narg != 4) error->all(FLERR, "Incorrect args for pair_coeff of pair_style dpd");
	if (!allocated) allocate();

	// find the lower and upper bound of type set in the pair_coeff
	force->bounds(arg[0], particle->ntypes, ilo, ihi);
	force->bounds(arg[1], particle->ntypes, jlo, jhi);

	double cut_one = cut_global;
	cut_one = atof(arg[3]);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			Ha[i][j] = atof(arg[2]);
			cut[i][j] = cut_one;
			cutsq[i][j] = cut_one * cut_one;
			setflag[i][j] = 1;        // "1" means this pair has been setup
			count++;
		}
	}

	if (count == 0) {
		error->all(FLERR, "Incorrest args for pair coefficients");
	}
}
