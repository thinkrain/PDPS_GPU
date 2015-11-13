/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(dpd,PairDPD)

#else

#ifndef PS_PAIR_DPD_H
#define PS_PAIR_DPD_H

#include "pair.h"

namespace PDPS_NS {

class PairDPD : public Pair {
public:
	PairDPD(class PDPS *);
	~PairDPD();

	void compute(int, int);
	void init_one(int, int);
	void set_style(int, char **);
	void set_coeff(int, char **);

protected:
	double cut_global, temperature;
	int seed;
	int rank;
	double **a0,**gamma;
	double **sigma;
	class RanMars *random;
	void allocate();
};

}

#endif
#endif
