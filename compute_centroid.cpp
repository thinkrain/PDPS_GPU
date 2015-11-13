/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "string.h"

#include "compute_centroid.h"
#include "error.h"
#include "modify.h"
#include "particle.h"
#include "update.h"

#include "parallel.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCentroid::ComputeCentroid(PDPS *ps, int narg, char **arg) : Compute(ps, narg, arg)
{
	if (narg < 3) error->all(FLERR,"Illegal compute centroid command");

	vector_flag = 1;

	size_vector = 3;

	vector = new double[3];
}

/* ---------------------------------------------------------------------- */

ComputeCentroid::~ComputeCentroid()
{

}

/* ---------------------------------------------------------------------- */

void ComputeCentroid::init()
{

}

/* ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
------------------------------------------------------------------------- */

void ComputeCentroid::compute_vector()
{
	int i = 0;
	double **x = particle->x;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	invoked_vector = update->ntimestep;

	for (i = 0; i < 3; i++) {
		centroid[i] = 0.0;
		vector[i] = 0.0;
	}

	int num = 0;
	for (i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			centroid[0] += x[i][0];
			centroid[1] += x[i][1];
			centroid[2] += x[i][2];
			num++;
		}
	}

	int num_all;

    MPI_Allreduce(&num, &num_all, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&centroid[0], &vector[0], 3, MPI_DOUBLE, MPI_SUM, mworld);

	for (i = 0; i < 3; i++) {
		if (num_all > 0) vector[i] /= num_all;
	}
}

