/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "string.h"

#include "integrate.h"
#include "v_verlet.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

Integrate::Integrate(PDPS *ps, int narg, char **arg) : Pointers(ps)
{
	style = NULL;

	int n = strlen(arg[0]) + 1;
	style = new char[n];
	strcpy(style,arg[0]);
}

/* ---------------------------------------------------------------------- */

Integrate::~Integrate()
{

}

/* ---------------------------------------------------------------------- */

void Integrate::init()
{
  // allow pair and Kspace compute() to be turned off via modify flags

}

/* ---------------------------------------------------------------------- */

void Integrate::ev_setup()
{
  // allow pair and Kspace compute() to be turned off via modify flags

}

/* ----------------------------------------------------------------------
   set eflag,vflag for current iteration
   invoke matchstep() on all timestep-dependent computes to clear their arrays
   eflag/vflag based on computes that need info on this ntimestep
   eflag = 0 = no energy computation
   eflag = 1 = global energy only
   eflag = 2 = per-atom energy only
   eflag = 3 = both global and per-atom energy
   vflag = 0 = no virial computation (pressure)
   vflag = 1 = global virial with pair portion via sum of pairwise interactions
   vflag = 2 = global virial with pair portion via F dot r including ghosts
   vflag = 4 = per-atom virial only
   vflag = 5 or 6 = both global and per-atom virial
------------------------------------------------------------------------- */

void Integrate::ev_set(bigint ntimestep)
{
	int i, flag;

	eflag = 1;
	vflag = 1;
}
