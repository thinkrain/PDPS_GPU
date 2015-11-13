/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "string.h"

#include "error.h"
#include "fix.h"
#include "group.h"
#include "parallel.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

Fix::Fix(PDPS *ps, int narg, char **arg) : Pointers(ps)
{
	comm_forward = comm_reverse = 0;

	name = NULL;
	style = NULL;

	file = NULL;
	fname = NULL;

	nevery = 1;

	int n = strlen(arg[0]) + 1;
	name = new char[n];
	strcpy(name,arg[0]); 

	gid = group->find_group(arg[1]);
	if (gid == -1) {
		char str[128];
		sprintf(str,"Cannot find group id: %s",arg[1]);
		error->all(FLERR,str);
	}
	groupbit = group->bitmask[gid];

	n = strlen(arg[2]) + 1;
	style = new char[n];
	strcpy(style,arg[2]);

	procid = parallel->procid;
}

/* ---------------------------------------------------------------------- */

Fix::~Fix()
{
	delete[] name;
	delete[] style;
	name = NULL; 
	style = NULL;
}
