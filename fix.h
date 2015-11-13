/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_FIX_H
#define PS_FIX_H

#include "pointers.h"

namespace PDPS_NS {

class Fix : protected Pointers {
public:
	char *name,*style;
	int gid, groupbit;

	int invoked_flag;
	double vector[6];

	int comm_forward;                     // size of forward communication (0 if none)
	int comm_reverse;                     // size of reverse communication (0 if none)

	Fix(class PDPS *, int, char **);
	virtual ~Fix();

	virtual void init() {}
	virtual void setup() {}
	virtual void pre_integrate() {}
	virtual void initial_integrate() {}
	virtual void final_integrate() {}
	virtual void post_integrate() {}
	virtual void pre_force() {}
	virtual void post_force() {}
	virtual void end_of_step() {}
	virtual int setmask() = 0;

protected:
	int evflag;
	int vflag_global,vflag_atom;

	int procid;

	FILE *file;                             // output specific pair related data into the file
	char *fname;                            // file name
	int nevery;                             // frequency to output file

};

namespace FixConst {
  static const int INITIAL_INTEGRATE =       1<<0;
  static const int POST_INTEGRATE =          1<<1;
  static const int PRE_EXCHANGE =            1<<2;
  static const int PRE_NEIGHBOR =            1<<3;
  static const int PRE_FORCE =               1<<4;
  static const int POST_FORCE =              1<<5;
  static const int FINAL_INTEGRATE =         1<<6;
  static const int END_OF_STEP =             1<<7;
  static const int THERMO_ENERGY =           1<<8;
  static const int INITIAL_INTEGRATE_RESPA = 1<<9;
  static const int POST_INTEGRATE_RESPA =    1<<10;
  static const int PRE_FORCE_RESPA =         1<<11;
  static const int POST_FORCE_RESPA =        1<<12;
  static const int FINAL_INTEGRATE_RESPA =   1<<13;
  static const int MIN_PRE_EXCHANGE =        1<<14;
  static const int MIN_PRE_FORCE =           1<<15;
  static const int MIN_POST_FORCE =          1<<16;
  static const int MIN_ENERGY =              1<<17;
  static const int POST_RUN =                1<<18;
  static const int FIX_CONST_LAST =          1<<19;
  static const int PRE_INTEGRATE =           1<<20;
}

}

#endif
