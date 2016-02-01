/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdio.h"

#include "error.h"
#include "fix_sph.h"
#include "force.h"
#include "particle.h"
#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSPH::FixSPH(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	if (narg != 3) error->all(FLERR,"Illegal fix sph command");
	 if ((particle->ee_flag != 1) || (particle->rho_flag != 1))
    error->all(FLERR,"fix sph command requires particle_style with both energy and density");

}

/* ---------------------------------------------------------------------- */

void FixSPH::init()
{
	dtv = update->dt;
	dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixSPH::initial_integrate()
{
	double dtfm;

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	double **vest = particle->vest;
	double *rho = particle->rho;
	double *drho = particle->drho;
	double *e = particle->e;
	double *de = particle->de;
	double *mass = particle->mass;
	double *rmass = particle->rmass;
	int *type = particle->type;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;


	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			if (particle->rmass_flag) dtfm = dtf / rmass[i];
			else dtfm = dtf / mass[type[i]];

			e[i] += dtf * de[i];			// half-step update of particle internal energy
		//	rho[i] += dtf * drho[i];		// ... and density

			vest[i][0] = v[i][0] + 2.0 * dtfm * f[i][0];
			vest[i][1] = v[i][1] + 2.0 * dtfm * f[i][1];
			vest[i][2] = v[i][2] + 2.0 * dtfm * f[i][2];

			// update velocity at t + 0.5*dt
			v[i][0] += dtfm * f[i][0];
			v[i][1] += dtfm * f[i][1];
			v[i][2] += dtfm * f[i][2];
			// update position at time t + dt
			x[i][0] += dtv * v[i][0];
			x[i][1] += dtv * v[i][1];
			x[i][2] += dtv * v[i][2];
			if (groupbit == 4){
				f[i][2] = f[i][2];
			}
		}
	}
	

}

/* ---------------------------------------------------------------------- */

void FixSPH::final_integrate()
{
	double dtfm;

	// update v of particles in group

	double **x = particle->x;
	double **v = particle->v;
	double **f = particle->f;
	double **vest = particle->vest;
	double *rho = particle->rho;
	double *drho = particle->drho;
	double *e = particle->e;
	double *de = particle->de;
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
			if (groupbit == 4){
				f[i][2] = f[i][2];
			}

			e[i] += dtf * de[i];
		//	rho[i] += dtf * drho[i];
		}
	}
}

/* ---------------------------------------------------------------------- */

int FixSPH::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPH::setup_pre_force(int vflag)
{
  // set vest equal to v 
  double **v = particle->v;
  double **vest = particle->vest;
  int *mask = particle->mask;
  int nlocal = particle->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vest[i][0] = v[i][0];
      vest[i][1] = v[i][1];
      vest[i][2] = v[i][2];
    }
  }
}