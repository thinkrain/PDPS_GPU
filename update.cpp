/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "domain.h"
#include "error.h"
#include "force.h"
#include "integrate.h"
#include "neighbor.h"
#include "region.h"
#include "update.h"
#include "v_verlet.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

Update::Update(PDPS *ps) : Pointers(ps)
{
	ntimestep = 0;

	integrate_style = NULL;
	integrate = NULL;
	//integrate = new Integrate(ps);

	units_flag = 0;
}

/* ---------------------------------------------------------------------- */

Update::~Update()
{
	delete integrate_style;
	delete integrate;
	integrate = NULL;
	integrate_style = NULL;
}

/* ---------------------------------------------------------------------- */

void Update::init()
{
	integrate->init();
}

/* ---------------------------------------------------------------------- */

void Update::dynamic_check()
{
	for (int i = 0; i < domain->nregions; i++) {
		if (domain->regions[i]->dynamic_flag) {
			domain->regions[i]->dynamic_check();
		}
	}
}

/* ----------------------------------------------------------------------
						  Set units
------------------------------------------------------------------------- */

void Update::set_units(const char *style)
{
	// physical constants from:
	// http://physics.nist.gov/cuu/Constants/Table/allascii.txt
	// using thermochemical calorie = 4.184 J

	if (strcmp(style,"lj") == 0) {
		force->boltz = 1.0;
		//force->boltz = 8.300365358821798e-06;
		force->hplanck = 0.18292026;  // using LJ parameters for argon
		force->mvv2e = 1.0;
		force->ftm2v = 1.0;
		force->mv2d = 1.0;
		force->nktv2p = 1.0;
		force->qqr2e = 1.0;
		force->qe2f = 1.0;
		force->vxmu2f = 1.0;
		force->xxt2kmu = 1.0;
		force->e_mass = 0.0;    // not yet set
		force->hhmrr2e = 0.0;
		force->mvh2r = 0.0;
		force->angstrom = 1.0;
		force->femtosecond = 1.0;
		force->qelectron = 1.0;

		dt = 0.005;
		neighbor->rskin = 0.3;
		units_flag = 1;
	}

	if (units_flag == 0) {
		char str[128];
		sprintf(str, "units style %s is invalid", style);
		error->all(FLERR,str);
	}
}

/* ----------------------------------------------------------------------
						 create integrate class
------------------------------------------------------------------------- */

void Update::create_integrate(int narg, char **arg)
{
	if (narg < 1) error->all(FLERR,"Illegal run_style command");

	delete [] integrate_style;
	integrate_style = NULL;
	delete integrate;
	integrate = NULL;
	
	new_integrate(arg[0], narg, arg);
}

/* ----------------------------------------------------------------------
						create new integrate class
------------------------------------------------------------------------- */

void Update::new_integrate(char *style, int narg, char **arg)
{
	int success = 0;
	
#define INTEGRATE_CLASS
#define IntegrateStyle(key,Class) \
    if (strcmp(style,#key) == 0) integrate = new Class(ps,narg,arg);
#include "style_integrate.h"
#undef IntegrateStyle
#undef INTEGRATE_CLASS
	
}

/* ---------------------------------------------------------------------- */

void Update::reset_timestep(int narg, char **arg)
{
	if (narg != 1) error->all(FLERR,"Illegal reset_timestep command");
	int newstep = atoi(arg[0]);
	if (newstep < 0) error->all(FLERR,"Timestep must be >= 0");
	// if new step is too big
	ntimestep = newstep;
}
