/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph_taitwater,PairSPH_TAITWATER)

#else

#ifndef PS_PAIR_SPH_TAITWATER_H
#define PS_PAIR_SPH_TAITWATER_H

#include "pair.h"

namespace PDPS_NS {

class PairSPH_TAITWATER : public Pair {
public:
	PairSPH_TAITWATER(class PDPS *);
	~PairSPH_TAITWATER();

	void compute(int, int);
	void set_style(int, char **);
	void set_coeff(int, char **);
	void init_one(int, int);
    double single(int, int, int, int, double, double, double, double &);

protected:
	double cut_global, temperature;
	int seed;
	int rank;
	double **cutd, **cutdsq;
	class RanMars *random;
	double *rho0, *soundspeed, *B;
	double **viscosity;
	int first, nstep;

	void allocate();
};

}

#endif
#endif
