/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"

#include "compute.h"
#include "force.h"
#include "integrate.h"
#include "modify.h"
#include "neighbor.h"
#include "parallel.h"
#include "post_processor.h"
#include "run.h"
#include "timer.h"
#include "update.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

Run::Run(PDPS *ps) : Pointers(ps) {}

/* ---------------------------------------------------------------------- */

void Run::command(int narg, char **arg)
{
	int nsteps;           // number of steps to be run
	int runflag;          // previous run flag: 0 = first run 1 = not first run

	nsteps = atoi(arg[0]);  // store number of steps to be run

	update->firststep = update->ntimestep;
	update->laststep = update->ntimestep + nsteps;

	if (update->ntimestep == 0) runflag = 0;
	else runflag = 1;

	ps->init();
	
	update->integrate->setup();     // setup before run

	timer->stamp_start(TIME_LOOP);
    update->integrate->run(nsteps);
	timer->stamp_end(TIME_LOOP);

	postprocessor->finalize();
}
