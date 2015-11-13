/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph_taitwater_morris,PairSPH_TAITWATER_MORRIS)

#else

#ifndef PS_PAIR_SPH_TAITWATER_MORRIS_H
#define PS_PAIR_SPH_TAITWATER_MORRIS_H

#include "pair.h"

namespace PDPS_NS {

class PairSPH_TAITWATER_MORRIS : public Pair {
public:
	PairSPH_TAITWATER_MORRIS(class PDPS *);
	~PairSPH_TAITWATER_MORRIS();

	void compute(int, int);
	void set_style(int, char **);
	void set_coeff(int, char **);
	void init_one(int, int);
    double single(int, int, int, int, double, double, double, double &);

protected:
	//	double **cutd, **cutdsq;
	//	class RanMars *random;
	double *rho0, *soundspeed, *B;
	double **viscosity;
	int first;

	void allocate();
};

}

#endif
#endif

