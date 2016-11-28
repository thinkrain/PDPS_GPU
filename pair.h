/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_PAIR_H
#define PS_PAIR_H

#include "pointers.h"

namespace PDPS_NS {

class Pair : protected Pointers {
public:
	char *style;                            // store pair style's name 
	int pair_id;

	double eng_vdwl,eng_coul;               // accumulated energies
	int **setflag;                          // 0/1 = whether each i,j has been set

	int evflag;                             // energy,virial settings
	double virial[6];                       // accumulated virial

	int newton_pair;                        // 1: yes;  0: no

	double cut_max;                         // max cut-off in this pairs
	double cut_min;                         // min cut-off in this pair
	double **cut, **cutsq;                  // cut-off and its square for each pair
	double *rho_local;

	int comm_forward;                       // size of forward communication (0 if none)
	int comm_reverse;                       // size of reverse communication (0 if none)

	//	Parameters for GPU
	int *devSetflag;
	double *devCutsq;
	int *hostSetflag;
	double *hostCutsq;

	class PairList *pair_list;

	Pair(class PDPS*);
	~Pair();

	virtual void compute(int, int) = 0;
	void init();
	virtual void setup();
	virtual void init_one(int, int) {}
	virtual void set_style(int, char **) = 0;
	virtual void set_coeff(int, char **) = 0;
	virtual int pack_reverse_comm(int, int, double *) { return 0; }
	virtual void unpack_reverse_comm(int, int *, double *) {}
	virtual int pack_forward_comm(int, int *, double *) { return 0; }
	virtual void unpack_forward_comm(int, int, double *) {}
	// energy (has to be public, so that pair_style can call
	void ev_tally(int, int, int, int, double, double, double,
                double, double, double);
	  

protected:
	int allocated;                // 0/1 = whether arrays are allocated
	int procid;

	FILE *file;                             // output specific pair related data into the file
	char *fname;                            // file name
	int nevery;                             // frequency to output file
	

	virtual void ev_setup(int, int);
	void ev_unset();
};

}

#endif
