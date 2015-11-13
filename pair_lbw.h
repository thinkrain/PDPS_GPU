/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lb/willet, PairLBW)

#else

#ifndef PS_PAIR_LB_H
#define PS_PAIR_LB_H

#include "pair.h"

namespace PDPS_NS {

class PairLBW : public Pair {
public:
	PairLBW(class PDPS *);
	~PairLBW();

	void compute(int, int);
	void init_one(int, int);
	void set_style(int, char **);
	void set_coeff(int, char **);

protected:
	double cut_global;

	double sij_min;
	double gamma;                                                 // surface tension
	double Vpcb, ln_Vpcb, Vpcb_square_root, Vpcb_cube_root;       // volume per capillary bridge

	double phi;                                                   // contact angle:  default value = 0

	double mu;                                                    // dynamic viscosity
	int drag_flag;                                                // drag force flag

	int ncollisions, ncollisions_total;
	double liquid_volume, liquid_volume_total;                    // liquid volume
	
	void allocate();

	double ff1();
	double ff2();
	double ff3();
	double ff4();
};
}

#endif
#endif
