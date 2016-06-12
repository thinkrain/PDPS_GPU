/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdio.h"

#include "error.h"
#include "fix_sph_stationary.h"
#include "force.h"
#include "particle.h"
#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSPH_STATIONARY::FixSPH_STATIONARY(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if ((particle->ee_flag != 1) || (particle->rho_flag != 1))
		error->all(FLERR, "fix meso/stationary command requires particle_style with both energy and density, e.g. meso");

	if (narg != 3)
		error->all(FLERR, "Illegal number of arguments for fix sph stationary command");

}

/* ---------------------------------------------------------------------- */

void FixSPH_STATIONARY::init()
{
	dtv = update->dt;
	dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixSPH_STATIONARY::initial_integrate()
{
	double *rho = particle->rho;
	double *drho = particle->drho;
	double *e = particle->e;
	double *de = particle->de;

	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	int i;

	for (i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			e[i] += dtf * de[i]; // half-step update of particle internal energy
//			rho[i] += dtf * drho[i]; // ... and density
		}
	}


}

/* ---------------------------------------------------------------------- */

void FixSPH_STATIONARY::final_integrate()
{
	double *e = particle->e;
	double *de = particle->de;
	double *rho = particle->rho;
	double *drho = particle->drho;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			e[i] += dtf * de[i];
//			rho[i] += dtf * drho[i];
		}
	}
}

/* ---------------------------------------------------------------------- */

int FixSPH_STATIONARY::setmask()
{
	int mask = 0;
	mask |= INITIAL_INTEGRATE;
	mask |= FINAL_INTEGRATE;
	return mask;
}

/* ---------------------------------------------------------------------- */
