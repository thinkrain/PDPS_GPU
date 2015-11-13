/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph_heatconduction,PairSPH_HEATCONDUCTION)

#else

#ifndef PS_PAIR_SPH_HEATCONDUCTION_H
#define PS_PAIR_SPH_HEATCONDUCTION_H

#include "pair.h"

namespace PDPS_NS {

class PairSPH_HEATCONDUCTION : public Pair {
public:
	PairSPH_HEATCONDUCTION(class PDPS *);
	~PairSPH_HEATCONDUCTION();

	void compute(int, int);
	void set_style(int, char **);
	void set_coeff(int, char **);
	void init_one(int, int);
    double single(int, int, int, int, double, double, double, double &);

protected:
	double **alpha;
	void allocate();
	int first, nstep;
};

}

#endif
#endif
