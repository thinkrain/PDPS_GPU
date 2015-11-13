/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(mdpd,PairMDPD)

#else

#ifndef PS_PAIR_MDPD_H
#define PS_PAIR_MDPD_H

#include "pair.h"

namespace PDPS_NS {

class PairMDPD : public Pair {
public:
	PairMDPD(class PDPS *);
	~PairMDPD();

	void compute(int, int);
	void init_one(int, int);
	void set_style(int, char **);
	void set_coeff(int, char **);

protected:
	double cut_global, temperature;
	int seed;
	int rank;
	double **cutd, **cutdsq;
	double **a0, **b0, **gamma;
	double **sigma;
//	double *rho_local;               // local density
	int inum_old;
	class RanMars *random;
	void allocate();
};

}

#endif
#endif
