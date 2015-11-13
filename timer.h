/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_TIMER_H
#define PS_TIMER_H

#include "pointers.h"

enum{TIME_LOOP,TIME_PAIR,TIME_BOND,TIME_KSPACE,TIME_NEIGHBOR,
     TIME_COMM,TIME_OUTPUT,TIME_START,TIME_END,TIME_N};

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
