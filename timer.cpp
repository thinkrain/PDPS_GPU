/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"

#include "memory.h"
#include "timer.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

Timer::Timer(PDPS *ps) : Pointers(ps)
{
	start_time = 0.0;
	//elapsed_time = 0.0;

	time = NULL;
	time_all = NULL;

	time = memory->create(time,TIME_N,"Timer: time");
	time_all= memory->create(time_all,TIME_N,"Timer: time_all");

	for (int i = 0; i < TIME_N; i++) {
		time[i] = 0.0;
		time_all[i] = 0.0;
	}
}

/* ---------------------------------------------------------------------- */

Timer::~Timer()
{
	memory->destroy(time);
}

/* ---------------------------------------------------------------------- */

void Timer::init()
{
	for (int i = 0; i < TIME_N; i++) time[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

void Timer::stamp()
{
	// uncomment if want synchronized timing
	// MPI_Barrier(world);
	previous_time = MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

void Timer::stamp(int which)
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
	double current_time = MPI_Wtime();
	time[which] += current_time - previous_time;
	time_all[which] += current_time - previous_time;
	previous_time = current_time;
}

/* ---------------------------------------------------------------------- */

void Timer::stamp_start(int which)
{
	MPI_Barrier(mworld);
	start_time = MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

void Timer::stamp_end(int which)
{
	MPI_Barrier(mworld);
	double end_time = MPI_Wtime();
	time[which] += end_time - start_time;
	time_all[which] += end_time - start_time;
}
