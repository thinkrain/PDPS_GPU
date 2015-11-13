/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef ANALYZE_CLASS

AnalyzeStyle(ave/time,AnalyzeAveTime)

#else

#ifndef PS_ANALYZE_AVE_TIME_H
#define PS_ANALYZE_AVE_TIME_H

#include "analyze.h"

namespace PDPS_NS {

class AnalyzeAveTime : public Analyze {
public:
	AnalyzeAveTime(class PDPS *, int, char **);
	virtual ~AnalyzeAveTime();
	
	void init();

	void invoke_analyze();

private:
	int nevery, nrepeat, nfreq;
	int nstart;         // smallest timestep to do average 
	int cid;            // compute id

	int icol;           // global icol

	typedef void (AnalyzeAveTime::*FnPtr)();
	void allocate();
	void parse_fields(int, char **);
	void addfield(const char *, FnPtr, int, int);
	FnPtr *field_func;                // list of ptrs to functions

	// computes
	class Compute **compute;

	void write_array();
	void clean_array();

	// compute different kinds of properties

	void compute_compute();

	void compute_step();
	void compute_x();
	void compute_y();
	void compute_z();
	void compute_vx();
	void compute_vy();
	void compute_vz();
	void compute_fx();
	void compute_fy();
	void compute_fz();
	void compute_wx();
	void compute_wy();
	void compute_wz();
	void compute_max_vx();
	void compute_max_vy();
	void compute_max_vz();
	void compute_min_vx();
	void compute_min_vy();
	void compute_min_vz();
};

}

#endif
#endif
