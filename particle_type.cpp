/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "domain.h"
#include "particle_type.h"

using namespace PDPS_NS;

#define DELTA 10000
#define EPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

ParticleType::ParticleType(PDPS *ps, int narg, char **arg) : Pointers(ps)
{
	nmax = 0;
}

/* ---------------------------------------------------------------------- */

ParticleType::~ParticleType()
{
	
}

/* ---------------------------------------------------------------------- */

void ParticleType::init()
{
	deform_vremap = domain->deform_vremap;
	deform_groupbit = domain->deform_groupbit;
	h_rate = domain->h_rate;
}
