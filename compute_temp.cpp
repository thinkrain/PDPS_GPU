/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdio.h"

#include "compute_temp.h"
#include "domain.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "parallel.h"
#include "particle.h"
#include "update.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTemp::ComputeTemp(PDPS *ps, int narg, char **arg) : Compute(ps, narg, arg)
{
	fix_dof = 0;                          // calculated in the fix (future work)
	vector = new double[6];

	scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeTemp::~ComputeTemp()
{
	delete [] vector;
	vector = NULL;
}

/* ----------------------------------------------------------------------
   Compute Scalar Temperature
   ---------------------------------------------------------------------- */

double ComputeTemp::compute_scalar()
{
	invoked_scalar = update->ntimestep;
	double t = 0.0;
	
	int *type = particle->type;
    double **v = particle->v;
	double *mass = particle->mass;
	double *rmass = particle->rmass;
	int rmass_flag = particle->rmass_flag;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;

	for (int i = 0; i < nlocal; i++){
		if(mask[i] & groupbit) {
			if (rmass_flag) {
				t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * rmass[i];
			}
			else t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * mass[type[i]];
		}
	}
	
	MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,mworld);
	scalar *= tfactor;
	return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTemp::dof_compute()
{
	int nparticles = group->gparticles[gid];
	dof = domain->dim * nparticles;
	dof -= fix_dof;
	if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
	else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

void ComputeTemp::init()
{
  fix_dof = 0;

  dof_compute();
}
