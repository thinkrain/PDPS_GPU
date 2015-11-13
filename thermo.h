/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_THERMO_H
#define PS_THERMO_H

#include "pointers.h"

namespace PDPS_NS {

class Thermo : protected Pointers {

public:
	int modified;
	char *style;

	Thermo(class PDPS *, int, char **);
	~Thermo();

	void header();
	void compute();
	void init();

private:
	int procid,nprocs;

    char *line;
    char **keyword;
	char **format;
	int *vtype;
    int nfield;

	int ifield;
	int nfield_initial;
	int *field2index;              // which compute,fix,variable calcs this field
	int *argindex1;
	int *argindex2;

	int bivalue;
	int ivalue;
	double dvalue;

	int ncomputes;                // # of Compute objects called by thermo

	class Compute **computes;     // list of ptrs to the Compute objects

	class Compute *temperature, *pressure, *pe;

	typedef void (Thermo::*FnPtr)();
	void addfield(const char *, FnPtr, int);
	FnPtr *vfunc;                // list of ptrs to functions

	void allocate();
	void compute_compute();      // functions that compute a single value
	void compute_step();
	void compute_temp();
	void compute_press();
	void parse_fields(char *);

	void lost_check();

	int int_between_brackets(char *&);
};
}
#endif
