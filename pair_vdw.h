/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(vdw, PairVDW)

#else

#ifndef PS_PAIR_VDW_H
#define PS_PAIR_VDW_H

#include "pair.h"

namespace PDPS_NS {

class PairVDW : public Pair {
public:
	PairVDW(class PDPS *);
	~PairVDW();

	void compute(int, int);
	void init_one(int, int);
	void set_style(int, char **);
	void set_coeff(int, char **);

protected:
	double cut_global;
	double **Ha;                 // Ha: Hamaker constant; 
	double lambda, b;            // lambda: dipole interaction wavelength; b: constant
	double delta_min;            // delta_min: minimum separation distance to avoid sigularity
	void allocate();
};
}

#endif
#endif
