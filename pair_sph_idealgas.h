/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph_idealgas,PairSPH_IDEALGAS)

#else

#ifndef PS_PAIR_SPH_IDEALGAS_H
#define PS_PAIR_SPH_IDEALGAS_H

#include "pair.h"

namespace PDPS_NS {

class PairSPH_IDEALGAS : public Pair {
public:
	PairSPH_IDEALGAS(class PDPS *);
	~PairSPH_IDEALGAS();

	void compute(int, int);
	void set_style(int, char **);
	void set_coeff(int, char **);
	void init_one(int, int);
    double single(int, int, int, int, double, double, double, double &);

protected:

	double **cutd, **cutdsq;
	double a2D, a3D, h;
	int first, nstep;
	double **viscosity;
	int cubic_flag, quintic_flag;
	void allocate();
};

}

#endif
#endif
