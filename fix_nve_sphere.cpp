/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdio.h"

#include "error.h"
#include "fix_nve_sphere.h"
#include "force.h"
#include "particle.h"
#include "update.h"

#define INERTIA 0.4          // moment of inertia prefactor for sphere

using namespace PDPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVESphere::FixNVESphere(PDPS *ps, int narg, char **arg) : FixNVE(ps, narg, arg)
{
	if (narg != 3) error->all(FLERR, "Illegal fix nve command");

	if (particle->sphere_flag == 0) error->all(FLERR, "Illegal particle type");
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::init()
{
	FixNVE::init();

	// check that all particles are finite-size spheres
	// no point particles allowed

	double *radius = particle->radius;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			if (radius[i] == 0.0) {
				error->all(FLERR, "Fix nve/sphere requires extended particles");
			}
		}
	}
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::initial_integrate()
{
	double dtfm, dtirotate;

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	double *mass = particle->mass;
	double *rmass = particle->rmass; 
	double **omega = particle->omega;
	double *radius = particle->radius;
	double **torque = particle->torque;
	int *type = particle->type;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;

	double dtfrotate = dtf / INERTIA;
	
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			dtfm = dtf / rmass[i];
			// update velocity at t + 0.5*dt
			v[i][0] += dtfm * f[i][0];
			v[i][1] += dtfm * f[i][1];
			v[i][2] += dtfm * f[i][2];
			// update position at time t + dt
			x[i][0] += dtv * v[i][0];
			x[i][1] += dtv * v[i][1];
			x[i][2] += dtv * v[i][2];

			dtirotate = dtfrotate / (radius[i] * radius[i] * rmass[i]);
			omega[i][0] += dtirotate * torque[i][0];
			omega[i][1] += dtirotate * torque[i][1];
			omega[i][2] += dtirotate * torque[i][2];
		}
	}
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::final_integrate()
{
	double dtfm, dtirotate;

	// update v of particles in group

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	double *mass = particle->mass;
	double **omega = particle->omega;
	double *radius = particle->radius;
	double *rmass = particle->rmass;
	double **torque = particle->torque;
	int *type = particle->type;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;

	double dtfrotate = dtf / INERTIA;
	
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			dtfm = dtf / rmass[i];
			// update velocity at time t + dt
			v[i][0] += dtfm * f[i][0];
			v[i][1] += dtfm * f[i][1];
			v[i][2] += dtfm * f[i][2];

			dtirotate = dtfrotate / (radius[i] * radius[i] * rmass[i]);
			omega[i][0] += dtirotate * torque[i][0];
			omega[i][1] += dtirotate * torque[i][1];
			omega[i][2] += dtirotate * torque[i][2];
		}
	}
}

/* ---------------------------------------------------------------------- */

int FixNVESphere::setmask()
{
	int mask = 0;
	mask |= INITIAL_INTEGRATE;
	mask |= FINAL_INTEGRATE;
	return mask;
}
