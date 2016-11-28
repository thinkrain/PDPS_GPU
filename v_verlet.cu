/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "error.h"
#include "domain.h"
#include "dump.h"
#include "force.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "parallel.h"
#include "particle.h"
#include "pair.h"
#include "timer.h"
#include "update.h"
#include "v_verlet.h"

#include "pdps_cuda.h"
#include "cuda_engine.h"
#include "device_launch_parameters.h"
#include "device_functions.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

V_Verlet::V_Verlet(PDPS *ps, int narg, char **arg) : Integrate(ps, narg, arg) 
{
	procid = parallel->procid;
}

/* ----------------------------------------------------------------------
   Initialization before run
------------------------------------------------------------------------- */

void V_Verlet::init()
{
	//triclinic = 0;
	//eflag = vflag = 0;

	ev_setup();
}


/* ----------------------------------------------------------------------
   Setup before run
------------------------------------------------------------------------- */

void V_Verlet::setup()
{
	particle->TransferC2G();
	// Copy inputs to device
	output->print("PDPS is setting up...\n");
	// setup domain and neighbor list
	domain->pbc();
	domain->reset_box();
	
	parallel->setup();
	neighbor->setup_cells();

	parallel->exchange();
	particle->lost_check();
	parallel->borders();
	neighbor->build();
	//neighbor->debug();
	neighbor->nflag = 1;
	
	// compute all forces
	ev_set(update->ntimestep);
	force->setup();
	force->clear();
	if (modify->n_pre_force) modify->pre_force();
	force->compute(eflag,vflag);
	parallel->reverse_comm();
	// output initial structure
	
	modify->setup();

	output->setup();


}

/* ----------------------------------------------------------------------
   Run for N steps
------------------------------------------------------------------------- */

void V_Verlet::run(int n)
{
	int nflag;           // neighbor build flag
    bigint ntimestep;

	int n_pre_integrate = modify->n_pre_integrate;
	int n_pre_force = modify->n_pre_force;
	int n_post_force = modify->n_post_force;
	int n_post_integrate = modify->n_post_integrate;
	int n_end_of_step = modify->n_end_of_step;
	int nanalyzes = modify->nanalyzes;

	// run
	ntimestep = update->ntimestep;

	// Output initial structure 
	for (int i = 0; i < n; i++) {
		ntimestep = ++update->ntimestep;
		ev_set(ntimestep);


		timer->stamp();
		// group or region may be dynamic
		update->dynamic_check();
		timer->stamp(TIME_SETUP);

		if (n_pre_integrate) modify->pre_integrate();
		// first integration of Verlet algorithm
		modify->initial_integrate();
		if (n_post_integrate) modify->post_integrate();
		timer->stamp(TIME_UPDATE);

		particle->TransferG2C();
		// build neighbor
		nflag = neighbor->decide();
		if (nflag == 1) {
		
			domain->pbc();
			if (domain->box_change) {
				timer->stamp();
				domain->reset_box();
				parallel->setup();
				neighbor->setup_cells();   // setup cells for creating linked list
				timer->stamp(TIME_SETUP);
			}
			timer->stamp();
			parallel->exchange();
			particle->lost_check();
			parallel->borders();
			timer->stamp(TIME_COMM);
			neighbor->build();
			//neighbor->debug();
			timer->stamp(TIME_NEIGHBOR);
		}
		else {
			timer->stamp();
			parallel->forward_comm();
			timer->stamp(TIME_COMM);
		}
		// force computation
		force->clear();
		if (n_pre_force) modify->pre_force();

		force->compute(eflag, vflag);
		timer->stamp();
		parallel->reverse_comm();
		timer->stamp(TIME_COMM);

		// force modifications
		if (n_post_force) modify->post_force();
		

		// second integration of Verlet algorithm
		modify->final_integrate();
		if (n_end_of_step) modify->end_of_step();
		timer->stamp(TIME_UPDATE);
		if (nanalyzes) modify->check_analyze();
		timer->stamp(TIME_ANALYZE);

		// output 
		
		timer->stamp();
		output->write();
		timer->stamp(TIME_OUTPUT);
	}
}
