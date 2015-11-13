/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "string.h"

#include "compute_pressure.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "pair.h"
#include "particle.h" 
#include "update.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

ComputePressure::ComputePressure(PDPS *ps, int narg, char **arg) : Compute(ps, narg, arg)
{
	if (narg < 4) error->all(FLERR,"Illegal compute pressure command");
	if (gid) error->all(FLERR,"Compute pressure must use group all");

	id_temp == NULL;
	temperature = NULL;

	// store temperature ID used by pressure computation
	// insure it is valid for temperature computation
	int n = strlen(arg[3]) + 1;
	id_temp = new char[n];
	strcpy(id_temp,arg[3]);	

	for (int i = 0; i < 6; i++) {
		virial[i] = 0.0;
	}

	scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputePressure::~ComputePressure()
{
	delete[] id_temp;
	id_temp = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputePressure::init()
{
	boltz = force->boltz;
	nktv2p = force->nktv2p;
	dim = domain->dim;

	// set temperature compute, must be done in init()
	// fixes could have changed or compute_modify could have changed it

	cid = modify->find_compute(id_temp);
	if (cid == -1) {
		char str[128];
		sprintf(str,"Cannot find group id: %s",id_temp);
		error->all(FLERR,str);
	}
	temperature = modify->compute[cid];
}

/* ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
------------------------------------------------------------------------- */

double ComputePressure::compute_scalar()
{
	invoked_scalar = update->ntimestep;

	double t;

	if (temperature->invoked_scalar != update->ntimestep)
		t = temperature->compute_scalar();
	else {
		t = temperature->scalar;
	}

	if (dim == 3) {
		volume = (domain->boxle[0] * domain->boxle[1] * domain->boxle[2]);
		virial_compute(3);
		scalar = (temperature->dof * boltz * t + virial[0] + virial[1] + virial[2]) / 3.0 / volume * nktv2p;
	}
	else {
		volume = domain->boxle[0] * domain->boxle[1];
		virial_compute(2);
		scalar = (temperature->dof * boltz * t + virial[0] + virial[1]) / 2.0 / volume * nktv2p;
	}

	return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputePressure::virial_compute(int n)
{
	int i, j, k;
	int counter = 0;
	double **x = particle->x;
	double **f = particle->f;
	double v[6];
	
	for (k = 0; k < n; k++) {
		v[k] = 0.0;
		virial[k] = 0.0;
	}

	for (i = 0; i < force->npairs; i++)
	for (k = 0; k < n; k++) {
		v[k] += force->pair[i]->virial[k];
	}
	
	MPI_Allreduce(v,virial,n,MPI_DOUBLE,MPI_SUM,mworld);
}
