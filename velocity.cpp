/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "compute.h"
#include "compute_temp.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "parallel.h"
#include "particle.h"
#include "random_park.h"
#include "velocity.h"

using namespace PDPS_NS;

#define DELTA 4

/* ---------------------------------------------------------------------- */

Velocity::Velocity(PDPS *ps) : Pointers(ps) {}

/* ---------------------------------------------------------------------- */

void Velocity::command(int narg, char** arg)
{
	int seed;
	double t_desired;

	temperature = NULL;
	
	dist_flag = 0;                      // default
	gid = group->find_group(arg[0]);    // find group id
	groupbit = group->bitmask[gid];     // store group's bitmask
	t_desired = atof(arg[2]);           // temperature target
	seed = atoi(arg[3]);                // seed to generate random number
    //v_style = arg[4];                 // how to create velocity profile

	if(!strcmp(arg[1],"create")) {
		if(!strcmp(arg[4],"uniform")) dist_flag = 0;
		else if(!strcmp(arg[4],"gaussian")) dist_flag = 1;
		create(t_desired, seed);
	}
	else if(!strcmp(arg[1],"set")) {
		set(narg, arg);
	}
	
}

/* ----------------------------------------------------------------------
					  Create Velocity Profile
------------------------------------------------------------------------- */

void Velocity::create(double t_desired, int seed)
{
    if (seed <= 0) error->all(FLERR, "the seed has to be a positive interger");

	double **v = particle->v;
	double *mass = particle->mass;
	double *rmass = particle->rmass;
	int rmass_flag = particle->rmass_flag;
	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	int dim = domain->dim;
	double vx, vy, vz, factor;
	RanPark *random;
	int m;

    int tflag = 0;
	if (temperature == NULL) {
		char **arg = new char*[3];
		arg[0] = (char *) "velocity_temp";
		arg[1] = group->name[gid];
		arg[2] = (char *) "temp";
		temperature = new ComputeTemp(ps, 3, arg);
		tflag = 1;
		delete [] arg;
	}
	temperature->init();
	
	// create an particle map if one doesn't exist already

	int mapflag = 0;
	if (particle->map_style == 0) {
		mapflag = 1;
		// it should be 1 
		particle->map_style = 2;
		particle->nghost = 0;
		particle->map_init();
		particle->map_set();
		m = particle->map(2);
	}

	if (particle->tag_enable == 0) {
		error->all(FLERR,"Cannot use velocity create loop all unless atoms have IDs");
	}
	
	random = new RanPark(ps, seed);
	
	int nparticles = static_cast<int> (particle->nparticles);
	for (int i = 1; i <= nparticles; i++) {
		if(dist_flag == 0) {
			vx = random->uniform();
			vy = random->uniform();
			vz = random->uniform();
		} else {
			vx = random->gaussian();
			vy = random->gaussian();
			vz = random->gaussian();
		}
		m = particle->map(i);
		if (m >= 0 && m < nlocal) {
			if (mask[m] & groupbit) {
				if (rmass_flag) factor = 1.0/sqrt(rmass[m]);
				else factor = 1.0/sqrt(mass[type[m]]);
				v[m][0] = vx * factor;
				v[m][1] = vy * factor;
				if(domain->dim == 3) v[m][2] = vz * factor;
				else v[m][2] = 0.0;
			}
		}
	} // for (int i = 0; i < nparticles; i++) 

	if (mapflag) {
		particle->map_delete();
		particle->map_style = 0;
		mapflag = 0;
	}

	double t = temperature->compute_scalar();
	rescale(t, t_desired);

	delete random;
	if (tflag) delete temperature;
}

/* ----------------------------------------------------------------------
   rescale velocities of group particles to t_new from t_old
------------------------------------------------------------------------- */

void Velocity::rescale(double t_old, double t_new)
{
    //if (t_old == 0.0) error->all(FLERR,FLERR,"Attempting to rescale a 0.0 temperature");
	
	double factor = sqrt(t_new/t_old);

	int *mask = particle->mask;
	double **v = particle->v;
	int nlocal = particle->nlocal;

	for(int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			v[i][0] *= factor;
			v[i][1] *= factor;
			v[i][2] *= factor;
		}
	}
}

/* ----------------------------------------------------------------------
   set velocities of group particles  
------------------------------------------------------------------------- */

void Velocity::set(int narg, char **arg)
{
	double **v = particle->v;
	double **omega = particle->omega;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	int omega_flag = particle->omega_flag;

	if ((omega_flag == 0 && narg != 5) || (omega_flag == 1 && narg != 8)) {
		error->all(FLERR, "Illegal velocity set command, check if rotation is enabled or not");
	}

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			v[i][0] = atof(arg[2]);
			v[i][1] = atof(arg[3]);
			if (domain->dim == 3) {
				v[i][2] = atof(arg[4]);
			}
			else {
				v[i][2] = 0.0;
			}
			if (omega_flag) {
				omega[i][0] = atof(arg[5]);
				omega[i][1] = atof(arg[6]);
				omega[i][2] = atof(arg[7]);
			}
		}
	}
}
