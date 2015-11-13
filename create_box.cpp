/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "create_box.h"
#include "domain.h"
#include "error.h"
#include "parallel.h"
#include "particle.h"
#include "region.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

CreateBox::CreateBox(PDPS *ps) : Pointers(ps) {}

/* ---------------------------------------------------------------------- */

void CreateBox::command(int narg, char **arg)
{
	int rid;

	if (domain->box_exist)
		error->all(FLERR,"Box has already been defined");
	if (narg != 2) 
		error->all(FLERR,"Illegal create_box command");

	if (domain->dim == 2 && domain->periodicity[2] == 0) 
		error->all(FLERR,"2D Simulation cannot have a nonperiodicity along Z coordinate");

	domain->box_exist = 1;

	rid = domain->find_region(arg[1]);         // box id in the region
	if (rid == -1) {
		char str[128];
		sprintf(str,"Create_box: Illegal region's name %s",arg[1]);
		error->all(FLERR,str);
	}

	particle->ntypes = atoi(arg[0]);   // store number of types
	 
	// Define box domain
	domain->boxlo[0] = domain->regions[rid]->extent_xlo;
	domain->boxhi[0] = domain->regions[rid]->extent_xhi;
	domain->boxlo[1] = domain->regions[rid]->extent_ylo;
	domain->boxhi[1] = domain->regions[rid]->extent_yhi;
	domain->boxlo[2] = domain->regions[rid]->extent_zlo;
	domain->boxhi[2] = domain->regions[rid]->extent_zhi;
	
	if (domain->boxlo[0] >= domain->boxhi[0] || \
		domain->boxlo[1] >= domain->boxhi[1] || \
		domain->boxlo[2] >= domain->boxhi[2])
		error->all(FLERR,"Simulation box bounds are invalid");

	domain->print_box("created");
	domain->set_initial_box();
	domain->set_global_box();
	parallel->set_proc_grid();
	domain->set_local_box();
}
