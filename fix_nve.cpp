/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdio.h"

#include "error.h"
#include "fix_nve.h"
#include "force.h"
#include "particle.h"
#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVE::FixNVE(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg != 3) error->all(FLERR,"Illegal fix nve command");
}

/* ---------------------------------------------------------------------- */

void FixNVE::init()
{
	dtv = update->dt;
	dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixNVE::initial_integrate()
{
	double dtfm;

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	double *mass = particle->mass;
	double *rmass = particle->rmass;
	int *type = particle->type;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			if (particle->rmass_flag) dtfm = dtf / rmass[i];
			else dtfm = dtf / mass[type[i]];
			// update velocity at t + 0.5*dt
			v[i][0] += dtfm * f[i][0];
			v[i][1] += dtfm * f[i][1];
			v[i][2] += dtfm * f[i][2];
			// update position at time t + dt
			x[i][0] += dtv * v[i][0];
			x[i][1] += dtv * v[i][1];
			x[i][2] += dtv * v[i][2];
		}
	}
	

}

/* ---------------------------------------------------------------------- */

void FixNVE::final_integrate()
{
	double dtfm;

	// update v of particles in group

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	double *mass = particle->mass;
	double *rmass = particle->rmass;
	int *type = particle->type;
    int *mask = particle->mask;
	int nlocal = particle->nlocal;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			if (particle->rmass_flag) dtfm = dtf / rmass[i];
			else dtfm = dtf / mass[type[i]];
			// update velocity at time t + dt
			v[i][0] += dtfm * f[i][0];
			v[i][1] += dtfm * f[i][1];
			v[i][2] += dtfm * f[i][2];
		}
	}
}

/* ---------------------------------------------------------------------- */

int FixNVE::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}
