/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "analyze.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "parallel.h"
#include "particle.h"
#include "update.h"

using namespace PDPS_NS;

#define DELTA 4

/* ---------------------------------------------------------------------- */

Analyze::Analyze(PDPS *ps, int narg, char **arg) : Pointers(ps)
{
	procid = parallel->procid;
	nprocs = parallel->nprocs;

	comm_forward = comm_reverse = 0;

	scalar_flag = 0;
	vector_flag = 0;
	array_flag = 0;

	flush_flag = 0;

	name = NULL;
	style = NULL;
	
	// array related
	array = NULL;
	array_total = NULL;
	array_ave = NULL;

	// fields related
	field_data_type = NULL;
	field_format = NULL;
	field_index = NULL;
	field_name = NULL;
	field_nrows = NULL;
	field_ncols = NULL;
	field_type = NULL;
	field_index = NULL;

	// option related
	fname = NULL;
	file = NULL;

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

Analyze::~Analyze()
{
	delete[] name;
	name = NULL;

	delete[] style;
	style = NULL;

	delete[] fname;
	fname = NULL;

	if (file && procid == 0) {
		fclose(file);
		file = NULL;
	}
}



