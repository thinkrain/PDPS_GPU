/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_force_H
#define PS_force_H

#include "pointers.h"

namespace PDPS_NS {

class Force : protected Pointers {
public:
	double boltz;                      // Boltzmann constant (eng/degree-K)
	double hplanck;                    // Planck's constant (energy-time)
	double mvv2e;                      // conversion of mv^2 to energy
	double ftm2v;                      // conversion of ft/m to velocity
	double mv2d;                       // conversion of mass/volume to density
	double nktv2p;                     // conversion of NkT/V to pressure
	double qqr2e;                      // conversion of q^2/r to energy
	double qe2f;                       // conversion of qE to force
	double vxmu2f;                     // conversion of vx dynamic-visc to force
	double xxt2kmu;                    // conversion of xx/t to kinematic-visc
	double dielectric;                 // dielectric constant
	double qqrd2e;                     // q^2/r to energy w/ dielectric constant
	double e_mass;                     // electron mass
	double hhmrr2e;                    // conversion of (hbar)^2/(mr^2) to energy
	double mvh2r;                      // conversion of mv/hbar to distance
                                       // hbar = h/(2*pi)
	double angstrom;                   // 1 angstrom in native units
	double femtosecond;                // 1 femtosecond in native units
	double qelectron;                  // 1 electron charge abs() in native units

	// pair style
	class Pair **pair;                 // multiple pairs
	int npairs;                        // number of pairs
	int **type2pair;                   // itype and jtype to force pair id

	double cut_max_global;             // max cut-off among all pairs 
	double cut_min_global;             // min cut-off among all pairs

	Force(class PDPS *);
	~Force();
	void create_pair(int, char **);
	
	void init();
	void setup();
	void bounds(char *, int, int &, int &);
	void clear();
	void compute(int, int);

	int find_pair_style(const char *);

private:
	int maxpairs;	
};

}

#endif
