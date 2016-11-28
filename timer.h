/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_TIMER_H
#define PS_TIMER_H

#include "pointers.h"

enum{TIME_LOOP,TIME_PAIR,TIME_NEIGHBOR,
	TIME_COMM, TIME_OUTPUT, TIME_START, TIME_END, TIME_SETUP, TIME_UPDATE, TIME_ANALYZE, TIME_COMPUTE, TIME_FIX1, TIME_FIX2, TIME_FIX3, TIME_FIX4, TIME_FIX5, TIME_FIX6,TIME_FIX7,TIME_FIX8, TIME_FIX9, TIME_FIX10,
	TIME_PAIR1, TIME_PAIR2, TIME_PAIR3, TIME_PAIR4, TIME_PAIR5, TIME_PREFIX1, TIME_PREFIX2, TIME_PREFIX3, TIME_PREFIX4, TIME_PREFIX5, TIME_PREFIX6, TIME_PREFIX7, TIME_PREFIX8, TIME_PREFIX9, TIME_PREFIX10, TIME_N,
};

namespace PDPS_NS {

class Timer : protected Pointers {
public:
	double *time;                   // time elapsed of single run     
	double *time_all;               // accumulated time elapsed of all runs

	Timer(class PDPS *);
	~Timer();

	void init();
	void stamp();
	void stamp(int);
	void stamp_start(int);
	void stamp_end(int);

private:
	double previous_time;
	double start_time;
	
};

}

#endif
