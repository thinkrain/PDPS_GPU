/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph_rhosumlocalAV, PairSPH_RHOSUMLOCALAV)

#else

#ifndef PS_PAIR_SPH_RHOSUMLOCALAV_H
#define PS_PAIR_SPH_RHOSUMLOCALAV_H

#include "pair.h"

namespace PDPS_NS {

	class PairSPH_RHOSUMLOCALAV : public Pair {
	public:
		PairSPH_RHOSUMLOCALAV(class PDPS *);
		~PairSPH_RHOSUMLOCALAV();

		void compute(int, int);
		void set_style(int, char **);
		void set_coeff(int, char **);
		void init_one(int, int);
		double single(int, int, int, int, double, double, double, double &);
		int pack_reverse_comm(int, int, double *);
		void unpack_reverse_comm(int, int *, double *);
		int pack_forward_comm(int, int *, double *);
		void unpack_forward_comm(int, int, double *);

	protected:
		int nstep, first;
		void allocate();
		double a2D, a3D, h;
		int sgid, sgroupbit;
		double *rho0, *soundspeed, *B;
		int cubic_flag, quintic_flag;
		int bcbit;
	};

}

#endif
#endif
