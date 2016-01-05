/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "compute_level.h"
#include "error.h"
#include "modify.h"
#include "particle.h"
#include "update.h"

#include "parallel.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

ComputeLevel::ComputeLevel(PDPS *ps, int narg, char **arg) : Compute(ps, narg, arg)
{
	if (narg < 5) error->all(FLERR,"Illegal compute level command");

	scalar_flag = 1;
	frequency = atoi(arg[3]);
	level = level_pre = 0.0; 
	levelgap = atof(arg[4]);


}

/* ---------------------------------------------------------------------- */

ComputeLevel::~ComputeLevel()
{

}

/* ---------------------------------------------------------------------- */

void ComputeLevel::init()
{

}

/* ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
------------------------------------------------------------------------- */

double ComputeLevel::compute_scalar()
{
	int i, j, k;
	double temp;
	double **x = particle->x;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;

	int level_compute = 0;
	int level_computeall = 0;

	if (update->ntimestep % frequency == 0){

		double highest = findhigh();
		if (highest > level_pre - levelgap){
			level_compute = 1;
			level = createtop();
		}
		MPI_Allreduce(&level_compute, &level_computeall, 1, MPI_INT, MPI_SUM, mworld);
		MPI_Allreduce(&level, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
		scalar = scalar / level_computeall;
		level_pre = scalar;
	}
	return scalar;
}


double ComputeLevel::findhigh(){
	double **x = particle->x;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	double highest = 0.0;
	for (int i = 0; i < nlocal; i++){
		if (mask[i] & groupbit) {
			if (x[i][2] > highest){
				highest = x[i][2];
			}
		}
	}
	return highest;
}

double ComputeLevel::createtop(){

	int i, j, k;
	double **x = particle->x;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	double top[100];
	for (i = 0; i < 100; i++)
		top[i] = 0.0;
	for (i = 0; i < nlocal; i++){
		if (mask[i] & groupbit) {
			// find the top 100 highest particles
			if (x[i][2] > top[99]){
				for (j = 99; j >= 0; j--){
					if (x[i][2] < top[j])
						break;
				}
				for (k = 99; k > j + 1; k--){
					top[k] = top[k - 1];
				}
				top[j + 1] = x[i][2];
			}
		}
	}

	return top[99];

}