/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "domain.h"
#include "dump.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "output.h"
#include "parallel.h"
#include "particle.h"
#include "style_dump.h"
#include "thermo.h"
#include "update.h"
#include "timer.h"
using namespace PDPS_NS;

#define DELTA 2

/* ----------------------------------------------------------------------
   initialize all output
------------------------------------------------------------------------- */

Output::Output(PDPS *ps) : Pointers(ps)
{
	maxdump = 0;
	ndumps = 0;

	dump = NULL;
	thermo = NULL;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

Output::~Output()
{
	for (int i = 0; i < ndumps; i++) {
		delete dump[i];
		dump[i] = NULL;
	}
	memory->sfree(dump);

	delete thermo;
	thermo = NULL;
}

/* ---------------------------------------------------------------------- */

void Output::init()
{
	if (thermo != NULL) thermo->init();

	for (int i = 0; i < ndumps; i++) {
		dump[i]->init();
	}
}

/* ---------------------------------------------------------------------- */

void Output::setup()
{
	for (int i = 0; i < ndumps; i++) {
		dump[i]->write();
	}

	if (thermo) {
		thermo->header();
		thermo->compute();
	}
}

/* ----------------------------------------------------------------------
   Add a dump
------------------------------------------------------------------------- */

void Output::add_dump(int narg, char **arg)
{
	if (ndumps == maxdump) {
		maxdump += DELTA;
		dump = (Dump **)
      memory->srealloc(dump,maxdump*sizeof(Dump *),"output:dump"); 
	}

	if (0) return;

#define DUMP_CLASS
#define DumpStyle(key,Class) \
	else if (strcmp(arg[2],#key) == 0) dump[ndumps] = new Class(ps,narg,arg);
#include "style_dump.h"
#undef DUMP_CLASS

	else error->all(FLERR, "Illegal dump style");

	ndumps++;
}

/* ----------------------------------------------------------------------
   Modify dump parameters
------------------------------------------------------------------------- */

void Output::modify_dump(int narg, char **arg)
{
	if (narg < 1) error->all(FLERR,"Illegal dump_modify command");

	int idump;
	for (idump = 0; idump < ndumps; idump++) {
		if (strcmp(arg[0],dump[idump]->name) == 0) break;
	}
	if (idump == ndumps) error->all(FLERR,"Cound not find dump_modify ID");

	dump[idump]->modify_params(narg-1,&arg[1]);
}

/* ----------------------------------------------------------------------
   Delete a dump 
------------------------------------------------------------------------- */

void Output::delete_dump(char *name)
{
	int idump;
	for (idump = 0; idump < ndumps; idump++) {
		if (strcmp(name,dump[idump]->name) == 0) break;
	}
	if (idump == ndumps) error->all(FLERR,"Cound not find undump ID");

	delete dump[idump];

	// move other dumps down in list one slot

	for (int i = idump+1; i < ndumps; i++) {
		dump[i-1] = dump[i];
	}
	ndumps--;
}

/* ----------------------------------------------------------------------
   Create thermo style
------------------------------------------------------------------------- */

void Output::create_thermo(int narg, char **arg)
{
	if (narg < 2) error->all(FLERR, "Illegal thermo_style command");
		
	if (thermo != NULL) {
		delete thermo;
		thermo = NULL;
	}
	thermo = new Thermo(ps, narg, arg);
}

/* ----------------------------------------------------------------------
   Set thermo command
------------------------------------------------------------------------- */

void Output::set_thermo(int narg, char **arg)
{
	if (narg != 1) error->all(FLERR, "Illegal thermo command");
	thermo_every = atoi(arg[0]);
}

/* ----------------------------------------------------------------------
						Write Dump to File
------------------------------------------------------------------------- */

void Output::write()
{
	// dump
	int ntimestep = update->ntimestep;
	
	bool transfered = false;
	for (int i = 0; i < ndumps; i++) {
		if (dump[i]->nevery != 0) {
			if ((ntimestep % dump[i]->nevery == 0)){				// check frequency
				if (!transfered){
					particle->TransferG2C();
					transfered = true;
				}
				dump[i]->write();
			}
		}
	}

	// thermo: compute all computes
	if (thermo) {
		if (ntimestep == 0) {
			thermo->header();
		}

		if (ntimestep%output->thermo_every==0) {
			thermo->compute();
		}
	}
}

/* ----------------------------------------------------------------------
   Write to screen and logfile
------------------------------------------------------------------------- */

void Output::print(const char *str)
{
	if (parallel->procid == 0) {
		if (screen)
			fprintf(screen,"%s",str);
		if (logfile)
			fprintf(logfile,"%s",str);
	}
}
