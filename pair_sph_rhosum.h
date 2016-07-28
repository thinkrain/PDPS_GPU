/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sph_rhosum, PairSPH_RHOSUM)

#else

#ifndef PS_PAIR_SPH_RHOSUM_H
#define PS_PAIR_SPH_RHOSUM_H

#include "pair.h"

namespace PDPS_NS {

	class PairSPH_RHOSUM : public Pair {
	public:
		PairSPH_RHOSUM(class PDPS *);
		~PairSPH_RHOSUM();

		void compute(int, int);
		void set_style(int, char **);
		void set_coeff(int, char **);
		void init_one(int, int);
		double single(int, int, int, int, double, double, double, double &);

		int pack_reverse_comm(int, int, double *);			//  reverse rho during computing so that particles gets information from neighbor processors
		void unpack_reverse_comm(int, int *, double *);
		int pack_forward_comm(int, int *, double *);		//	forward rho after reverse so that neighbor processor get the new correct rho just computed
		void unpack_forward_comm(int, int, double *);

	protected:
		int nstep, first;
		void allocate();
		double a2D, a3D, h;									//	parameters in SPH governing equation

		double *rho0, *soundspeed, *B;

		int cubic_flag, quintic_flag;						//	flag to use which kernel function

	};

}

#endif
#endif
