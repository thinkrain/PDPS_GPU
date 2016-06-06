/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph_couple,PairSPH_COUPLE)

#else

#ifndef PS_PAIR_SPH_COUPLE_H
#define PS_PAIR_SPH_COUPLE_H

#include "pair.h"

namespace PDPS_NS {

class PairSPH_COUPLE : public Pair {
public:
	PairSPH_COUPLE(class PDPS *);
	~PairSPH_COUPLE();

	void compute(int, int);
	void set_style(int, char **);
	void set_coeff(int, char **);
	void init_one(int, int);
    double single(int, int, int, int, double, double, double, double &);

protected:
	double cut_global, temperature;
	int seed;
	int rank;
	int local_kernel;				//	kernel length is 1(varying among particles) 0 (uniform)
	double **cutd, **cutdsq;
	double sigma;					//  parameter to control the relationship of kernel length and density
	class RanMars *random;
	double *rho0, *soundspeed, *B;
	double **viscosity;
	int first, nstep;
	int phase_f, phase_s;
	void allocate();
};

}

#endif
#endif
