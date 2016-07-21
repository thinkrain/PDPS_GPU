/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "compute_height.h"
#include "error.h"
#include "modify.h"
#include "particle.h"
#include "update.h"

#include "parallel.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

ComputeHeight::ComputeHeight(PDPS *ps, int narg, char **arg) : Compute(ps, narg, arg)
{
	if (narg < 8) error->all(FLERR,"Illegal compute height command");

	scalar_flag = 1;
	height = 0.0;
	xtemp = atof(arg[3]);
	ytemp = atof(arg[4]);
	z0 = atof(arg[5]);
	heightcut = atof(arg[6]);
	rho_ref = atof(arg[7]);


}

/* ---------------------------------------------------------------------- */

ComputeHeight::~ComputeHeight()
{

}

/* ---------------------------------------------------------------------- */

void ComputeHeight::init()
{

}

/* ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
------------------------------------------------------------------------- */

double ComputeHeight::compute_scalar()
{
	int i, j, k;
	double temp;
	double **x = particle->x;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	height = z0;


	height = findhigh();
	MPI_Allreduce(&height, &scalar, 1, MPI_DOUBLE, MPI_MAX, mworld);

	return scalar;
}


double ComputeHeight::findhigh(){
	double **x = particle->x;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	double *rho = particle->rho;
	double highest = 0.0;
	double cutsq = heightcut * heightcut;
	for (int i = 0; i < nlocal; i++){
		if (mask[i] & groupbit) {
			double rij = (x[i][0] - xtemp) * (x[i][0] - xtemp) + (x[i][1] - ytemp) * (x[i][1] - ytemp);
			if (rij < cutsq && x[i][2] > highest && rho[i] > rho_ref){
				highest = x[i][2];
			}
		}
	}
	return highest;
}
