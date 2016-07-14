/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph_lj,PairSPH_LJ)

#else

#ifndef PS_PAIR_SPH_LJ_H
#define PS_PAIR_SPH_LJ_H

#include "pair.h"

namespace PDPS_NS {

class PairSPH_LJ : public Pair {
public:
	PairSPH_LJ(class PDPS *);
	~PairSPH_LJ();

	void compute(int, int);
	void set_style(int, char **);
	void set_coeff(int, char **);
	void init_one(int, int);
    double single(int, int, int, int, double, double, double, double &);
	void LJEOS2(double, double, double, double *, double *);

protected:
	double K;
	double **viscosity;
	int first, nstep;
	void allocate();
};

}

#endif
#endif
