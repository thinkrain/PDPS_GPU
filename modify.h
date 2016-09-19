/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_MODIFY_H
#define PS_MODIFY_H

#include "stdio.h"
//#include "pointers.h"

namespace PDPS_NS {

class Modify : protected Pointers {
public:
	// analysis
	int nanalyzes, maxanalyze;
	class Analyze **analyze;  // 

	// fix
	int nfixes, maxfix;
	class Fix **fix;           // list of fixes
	int *fmask;       // bit mask for when each fix is applied

	// compute
	int ncomputes, maxcompute;   // list of computes
	class Compute **compute;

	// integrate
	int n_initial_integrate, n_final_integrate;
	int n_pre_integrate, n_post_integrate;
	int *list_initial_integrate, *list_final_integrate;
	int *list_pre_integrate, *list_post_integrate;

	// force
	int n_pre_force;
	int *list_pre_force;
	int n_post_force;
	int *list_post_force;

	// postprocessor
	int n_end_of_step;
	int *list_end_of_step;
	
	Modify(class PDPS *);
	virtual ~Modify();
	virtual void init();
	virtual void setup();
	
	// analyze
	void add_analyze(int, char **);
	void delete_analyze(const char *);
	int find_analyze(const char *);
	void check_analyze();
	
	// fix
	void add_fix(int, char **);
	void delete_fix(const char *);
	int find_fix(const char *);

	// compute
	void add_compute(int, char **);
	int find_compute(const char *);
	int find_compute_style(const char *);

	// integrate
	virtual void pre_integrate();
	virtual void initial_integrate();
	virtual void final_integrate();
	virtual void post_integrate();

	// force
	virtual void pre_force();
	virtual void post_force();

	// post processor
	virtual void end_of_step();

protected:
	void list_init(int, int &, int *&);
};

}

#endif
