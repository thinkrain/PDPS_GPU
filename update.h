/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_UPDATE_H
#define PS_UPDATE_H

#include "pointers.h"

namespace PDPS_NS {

class Update : protected Pointers {
public:
	double dt;                      // timestep
	double etol,ftol;               // minimizer tolerances on energy/force
	bigint ntimestep;                  // current step (dynamics or min iterations)
	bigint nsteps;                     // # of steps to run (dynamics or min iter)
	int whichflag;                  // 0 for unset, 1 for dynamics, 2 for min
	bigint firststep,laststep;         // 1st & last step of this run
	bigint beginstep,endstep;          // 1st and last step of multiple runs
	int first_update;               // 0 before initial update, 1 after

	int units_flag;                 // 1 or 0: units set or not
  
	class Integrate *integrate;
	char *integrate_style;          

	Update(class PDPS *);
	~Update();
	void init();
	void set_units(const char *);
	void create_integrate(int, char **);
	void create_minimize(int, char **);
	void reset_timestep(int, char **);

	void dynamic_check();
  
private:
	void new_integrate(char *, int, char **);

};

}

#endif
