/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "force.h"
#include "memory.h"
#include "pair.h"
#include "parallel.h"
#include "particle.h"

#include "pair_dpd.h"

using namespace PDPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

Pair::Pair(PDPS *ps) : Pointers(ps) 
{
	style = NULL;
	setflag = NULL;
	cut = NULL;
	cutsq = NULL;

	file = NULL;
	fname = NULL;
	
	allocated = 0;        // all array are not alloacted

	nevery = 1;

	comm_forward = comm_reverse = 0;

	pair_id = force->npairs;
	procid = parallel->procid;
}

/* ---------------------------------------------------------------------- */

Pair::~Pair()
{
	delete[] style;
	style = NULL;
}

/* ----------------------------------------------------------------------
   Setup for energy, virial computation
   see integrate::ev_set() for values of eflag (0-3) and vflag (0-6)
------------------------------------------------------------------------- */

void Pair::ev_setup(int eflag, int vfalg)
{
	int i, n;

	evflag = 1;
}

/* ---------------------------------------------------------------------- */

void Pair::init()
{
	int i,j;

	cut_max = 0.0;
	cut_min = BIG;
	// find cut_max and cut_min
	for (i = 1; i <= particle->ntypes; i++) {
		for (j = i; j <= particle->ntypes; j++) {
			if (setflag[i][j] == 1) {
				init_one(i,j);
				cut_max = MAX(cut_max,cut[i][j]);
				cut_min = MIN(cut_min,cut[i][j]);
			}
		}
	}
	
	for (i = 0; i < 6; i++) {
		virial[i] = 0.0;
	}
}

/* ---------------------------------------------------------------------- */

void Pair::setup()
{
	
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   need i < nlocal test since called by bond_quartic and dihedral_charmm
------------------------------------------------------------------------- */

void Pair::ev_tally(int i, int j, int nlocal, int newton_pair,
                    double evdwl, double ecoul, double fpair,
                    double delx, double dely, double delz)
{
	double evdwlhalf,ecoulhalf,epairhalf,v[6];

	/*if (j < nlocal) {
		v[0] = delx*delx*fpair;
		v[1] = dely*dely*fpair;
		v[2] = delz*delz*fpair;
		v[3] = delx*dely*fpair;
		v[4] = delx*delz*fpair;
		v[5] = dely*delz*fpair;
	}
	else {
		v[0] = 0.5*delx*delx*fpair;
		v[1] = 0.5*dely*dely*fpair;
		v[2] = 0.5*delz*delz*fpair;
		v[3] = 0.5*delx*dely*fpair;
		v[4] = 0.5*delx*delz*fpair;
		v[5] = 0.5*dely*delz*fpair;
	}*/

	v[0] = delx*delx*fpair;
	v[1] = dely*dely*fpair;
	v[2] = delz*delz*fpair;
	v[3] = delx*dely*fpair;
	v[4] = delx*delz*fpair;
	v[5] = dely*delz*fpair;

	virial[0] += v[0];
    virial[1] += v[1];
    virial[2] += v[2];
    virial[3] += v[3];
    virial[4] += v[4];
    virial[5] += v[5];
}
