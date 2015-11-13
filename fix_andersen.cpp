/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "error.h"
#include "fix_andersen.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair.h"
#include "particle.h"
#include "random_mars.h"
#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAndersen::FixAndersen(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg < 8) error->all(FLERR,"Illegal fix anderson command");

	temperature = atof(arg[3]);
	Gamma = atof(arg[4]);
	seed = atoi(arg[5]);

	int iarg = 6;

	int ilo,ihi,jlo,jhi;
	
	int ntypes = particle->ntypes;
	setflag = new int*[ntypes+1];
	for (int i = 0; i < ntypes+1; i++) {
		setflag[i] = new int[ntypes+1];
		for (int j = 0; j < ntypes+1; j++) setflag[i][j] = 0;
	}
	while (iarg < narg) {
		// find the lower and upper bound of type set in the pair_coeff
		force->bounds(arg[iarg],particle->ntypes,ilo,ihi);
		force->bounds(arg[iarg+1],particle->ntypes,jlo,jhi);
		for (int i = ilo; i <= ihi; i++) {
			for (int j = MAX(jlo,i); j <= jhi; j++) {
				setflag[i][j] = 1;        // "1" means this pair has been setup
			}
		}
		iarg += 2;
	}
	// initialize setflag
	for (int i = 1; i <= ntypes; i++) {
		for (int j = i; j <= ntypes; j++) {
			if (setflag[i][j] == 1) {
				setflag[j][i] = setflag[i][j];
			}
		}
	}

	random = NULL;
	random = new RanMars(ps, seed);
}

/* ---------------------------------------------------------------------- */

FixAndersen::~FixAndersen()
{
	delete random;
	random = NULL;
}

/* ---------------------------------------------------------------------- */

void FixAndersen::init()
{
	kbT = force->boltz*temperature;
	probability = Gamma * (update->dt);
}

/* ----------------------------------------------------------------------
   For the algorithm, see paper 
   "2006 Advantages of a Lowe-Andersen thermostat in molecular 
   dynamics simulations"
------------------------------------------------------------------------- */

void FixAndersen::end_of_step()
{
	int i, j, ii, jj, inum, jnum, itype, jtype;
	double mi, mj;
	double rijsq, rij, rij_inv; 
	double rix, riy, riz, rijx, rijy, rijz;
	double vix, viy, viz, vijx, vijy, vijz;
	double vij0, vijx0, vijy0, vijz0;
	double nx, ny, nz;
	double deltaij, deltaijx, deltaijy, deltaijz;

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	int *type = particle->type;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	int nall = nlocal + particle->nghost;
	
	int *ilist;
	int *jlist;
	int **firstneigh;
	int *numneigh;
	
	Pair **pair = force->pair;
	int **type2pair = force->type2pair;

	ilist = neighbor->neighlist->ilist;
	inum = neighbor->neighlist->inum;
	firstneigh = neighbor->neighlist->firstneigh;
	numneigh = neighbor->neighlist->numneigh;

	int pair_id;
	double randnum_uni;
	double randnum;
	double **cutsq;

	
	// loop for all local owned particles 
	for (ii = 0; ii < inum; ii++) {
		i = ilist[ii];
		itype = particle->type[i];
		rix = x[i][0];
		riy = x[i][1];
		riz = x[i][2];
		
		jlist = firstneigh[i];
		jnum = numneigh[i];
		if (particle->rmass_flag) mi = particle->rmass[i];
		else mi = particle->mass[itype];
		for (jj = 0; jj < jnum; jj++) {
			j = jlist[jj];
			jtype = type[j];
			if (setflag[itype][jtype]) {
				pair_id = type2pair[itype][jtype];
				cutsq = pair[pair_id]->cutsq;
				rijx = rix - x[j][0];
				rijy = riy - x[j][1];
				rijz = riz - x[j][2];
				rijsq = (rijx*rijx + rijy*rijy + rijz*rijz);
				rij = sqrt(rijsq);
				// withint the cutoff range
				if (rijsq < cutsq[itype][jtype]) {
					randnum_uni = random->uniform();
					if (randnum_uni > probability) {
						continue;
					}
					if (particle->rmass_flag) mj = particle->rmass[j];
					else mj = particle->mass[jtype];
					muij = mi*mj / (mi+mj);
					
					rij_inv = 1.0 / rij;
					nx = rijx * rij_inv;
					ny = rijy * rij_inv;
					nz = rijz * rij_inv;

					//vijx0 = random->gaussian()*sqrt(kbT/muij);
					//vijy0 = random->gaussian()*sqrt(kbT/muij);
					//vijz0 = random->gaussian()*sqrt(kbT/muij);
					vij0 = random->gaussian()*sqrt(kbT/muij);

					vijx = v[i][0] - v[j][0];
					vijy = v[i][1] - v[j][1];
					vijz = v[i][2] - v[j][2];

					deltaij = vij0 - (vijx*nx + vijy*ny + vijz*nz);
					//deltaij = (vijx0 - vijx)*nx + (vijy0 - vijy)*ny + (vijz0 - vijz)*nz;
					deltaijx = deltaij * nx;
					deltaijy = deltaij * ny;
					deltaijz = deltaij * nz;

					if (mask[i] & groupbit) {
						v[i][0] = v[i][0] + deltaijx*muij/mi;
						v[i][1] = v[i][1] + deltaijy*muij/mi;
						v[i][2] = v[i][2] + deltaijz*muij/mi;
					}
					if (mask[j] & groupbit) {
						v[j][0] = v[j][0] - deltaijx*muij/mj;
						v[j][1] = v[j][1] - deltaijy*muij/mj;
						v[j][2] = v[j][2] - deltaijz*muij/mj;
					}
				} // if (rij < cut[itype][jtype])
			}
		} // for (jj = 0; jj < nlocal; jj++)
	} // for (i = 0; i < nlocal; i++)
}

/* ---------------------------------------------------------------------- */

int FixAndersen::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}
