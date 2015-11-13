/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "compute.h"
#include "error.h"
#include "group.h"
#include "memory.h"

using namespace PDPS_NS;

#define DELTA 4

/* ---------------------------------------------------------------------- */

Compute::Compute(PDPS *ps, int narg, char **arg) : Pointers(ps)
{
	comm_forward = comm_reverse = 0;

	scalar_flag = vector_flag = array_flag = 0;

	invoked_scalar = invoked_vector = invoked_array = -1;

	name = NULL;
	style = NULL;

	int n = strlen(arg[0]) + 1;
	name = new char[n];
	strcpy(name,arg[0]);

	n = strlen(arg[2]) + 1;
	style = new char[n];
	strcpy(style,arg[2]);
	
	gid = group->find_group(arg[1]);
	if (gid == -1) {
		char str[128];
		sprintf(str,"Cannot find group id: %s",arg[1]);
		error->all(FLERR,str);
	}
	groupbit = group->bitmask[gid];
}

/* ---------------------------------------------------------------------- */

Compute::~Compute()
{
	delete[] name;
	name = NULL;

	delete[] style;
	style = NULL;
}
