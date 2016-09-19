/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"

#include "domain.h"
#include "error.h"
#include "memory.h"
#include "particle.h"
#include "particle_type_sphere.h"
#include "phy_const.h"

using namespace PDPS_NS;
using namespace PhyConst;

#define DELTA 10000
#define BIG 1e20

/* ---------------------------------------------------------------------- */

ParticleTypeSphere::ParticleTypeSphere(PDPS *ps, int narg, char **arg) : ParticleType(ps, narg, arg)
{
	comm_x_only = 1;
	comm_f_only = 0;
	size_forward = 3;
	size_reverse = 6;
	size_border = 8;
	size_velocity = 6;
	//size_data_atom = 7;
	//size_data_vel = 7;
	//xcol_data = 5;

	particle->atomic_flag = 0;
	particle->sphere_flag = 1;
	particle->radius_flag = 1;
	particle->rmass_flag = 1;
	particle->omega_flag = 1;
	particle->torque_flag = 1;
	particle->density_flag = 1;
}

/* ---------------------------------------------------------------------- */

ParticleTypeSphere::~ParticleTypeSphere()
{
	
}

/* ----------------------------------------------------------------------
   Grow particle arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void ParticleTypeSphere::grow(int n)
{
	if (n == 0) nmax += DELTA;
	else nmax = n;
	particle->nmax = nmax;

	if (nmax < 0 || nmax > BIG) {
		error->one(FLERR,"Too many particles to be created for this system");
	}

	tag = memory->grow(particle->tag, nmax, "particle: tag");
	type = memory->grow(particle->type, nmax, "Particle: type");
	mask = memory->grow(particle->mask, nmax, "Particle: mask");
	x = memory->grow(particle->x, nmax, 3, "particle: x");
	v = memory->grow(particle->v, nmax, 3, "particle: v");
	f = memory->grow(particle->f, nmax, 3, "particle: f");

	radius = memory->grow(particle->radius, nmax, "particle: radius");
	rmass = memory->grow(particle->rmass, nmax, "particle: rmass");
	omega = memory->grow(particle->omega, nmax, 3, "particle: omega");
	torque = memory->grow(particle->torque, nmax, 3, "particle torque");
}

/* ----------------------------------------------------------------------
Create particle
------------------------------------------------------------------------- */

void ParticleTypeSphere::create_particle(int itype, double *coord, double radi)
{
	int nlocal = particle->nlocal;
	if (nlocal == nmax) grow(0);

	tag[nlocal] = 0;
	type[nlocal] = itype;
	mask[nlocal] = 1;
	for (int i = 0; i < 3; i++) {
		x[nlocal][i] = coord[i];
		v[nlocal][i] = 0.0;
		f[nlocal][i] = 0.0;
		omega[nlocal][i] = 0.0;
		torque[nlocal][i] = 0.0;
	}
	radius[nlocal] = radi;
	rmass[nlocal] = 4.0/3 * PI * radi * radi * radi * particle->density[itype];

	particle->nlocal++;
}





/* ----------------------------------------------------------------------
   unpack one line from Particles section of data file initialize other 
   atom quantities
------------------------------------------------------------------------- */

void ParticleTypeSphere::data_particle(double *coord, char **values)
{
	int nlocal = particle->nlocal;
	if (nlocal == nmax) grow(0);

	tag[nlocal] = atoi(values[0]);
	if (tag[nlocal] <= 0)
		error->one(FLERR, "Invalid atom ID in Atoms section of data file");

	type[nlocal] = atoi(values[1]);
	if (type[nlocal] <= 0 || type[nlocal] > particle->ntypes)
		error->one(FLERR, "Invalid atom type in Atoms section of data file");

	radius[nlocal] = atof(values[2]);
	if (radius[nlocal] < 0.0) error->all(FLERR, "Invalid raidus in particles section of data file");

	density[type[nlocal]] = atof(values[3]);
	if (density[type[nlocal]] <= 0.0) error->all(FLERR, "Invalid density in particles section of data file");

	rmass[nlocal] = 4.0*PI / 3.0 * radius[nlocal] * radius[nlocal] * radius[nlocal] * density[type[nlocal]];

	x[nlocal][0] = coord[0];
	x[nlocal][1] = coord[1];
	x[nlocal][2] = coord[2];

	//image[nlocal] = imagetmp;

	mask[nlocal] = 1;
	v[nlocal][0] = 0.0;
	v[nlocal][1] = 0.0;
	v[nlocal][2] = 0.0;
	omega[nlocal][0] = 0.0;
	omega[nlocal][1] = 0.0;
	omega[nlocal][2] = 0.0;

	particle->nlocal++;
}

/* ----------------------------------------------------------------------
unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void ParticleTypeSphere::data_vel(int m, char **values)
{
	double **v = particle->v;
	v[m][0] = atof(values[0]);
	v[m][1] = atof(values[1]);
	v[m][2] = atof(values[2]);
	omega[m][0] = atof(values[3]);
	omega[m][1] = atof(values[4]);
	omega[m][2] = atof(values[5]);
}

/* ---------------------------------------------------------------------- */

int ParticleTypeSphere::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
	int i, j, m;
	double dx, dy, dz;

	m = 0;
	if (pbc_flag == 0) {
		for (i = 0; i < n; i++) {
			j = list[i];
			buf[m++] = x[j][0];
			buf[m++] = x[j][1];
			buf[m++] = x[j][2];
		}
	} 
	else {
		dx = pbc[0]*domain->xle;
		dy = pbc[1]*domain->yle;
		dz = pbc[2]*domain->zle;
  
		for (i = 0; i < n; i++) {
			j = list[i];
			buf[m++] = x[j][0] + dx;
			buf[m++] = x[j][1] + dy;
			buf[m++] = x[j][2] + dz;
		}
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int ParticleTypeSphere::pack_comm_vel(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
	int i, j, m;
	double dx, dy, dz, dvx, dvy, dvz;

	m = 0;
	if (pbc_flag == 0) {
		for (i = 0; i < n; i++) {
			j = list[i];
			buf[m++] = x[j][0];
			buf[m++] = x[j][1];
			buf[m++] = x[j][2];
			buf[m++] = v[j][0];
			buf[m++] = v[j][1];
			buf[m++] = v[j][2];
			buf[m++] = omega[j][0];
			buf[m++] = omega[j][1];
			buf[m++] = omega[j][2];
		}
	} // if (pbc_flag == 0)
	else {
		dx = pbc[0]*domain->xle;
		dy = pbc[1]*domain->yle;
		dz = pbc[2]*domain->zle;
      
		if (!deform_vremap) {
			for (i = 0; i < n; i++) {
				j = list[i];
				buf[m++] = x[j][0] + dx;
				buf[m++] = x[j][1] + dy;
				buf[m++] = x[j][2] + dz;
				buf[m++] = v[j][0];
				buf[m++] = v[j][1];
				buf[m++] = v[j][2];
				buf[m++] = omega[j][0];
				buf[m++] = omega[j][1];
				buf[m++] = omega[j][2];
			}
		} 
		else {
			dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
			dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
			dvz = pbc[2]*h_rate[2];
			for (i = 0; i < n; i++) {
				j = list[i];
				buf[m++] = x[j][0] + dx;
				buf[m++] = x[j][1] + dy;
				buf[m++] = x[j][2] + dz;
				if (mask[i] & deform_groupbit) {
					buf[m++] = v[j][0] + dvx;
					buf[m++] = v[j][1] + dvy;
					buf[m++] = v[j][2] + dvz;
				} 
				else {
					buf[m++] = v[j][0];
					buf[m++] = v[j][1];
					buf[m++] = v[j][2];
				}
				buf[m++] = omega[j][0];
				buf[m++] = omega[j][1];
				buf[m++] = omega[j][2];
			}
		} // else if (deform_vremap == 1)
	} // else if (pbc_flag != 0)
  return m;
}

/* ---------------------------------------------------------------------- */

void ParticleTypeSphere::unpack_comm(int n, int first, double *buf)
{
	int i,m,last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		x[i][0] = buf[m++];
		x[i][1] = buf[m++];
		x[i][2] = buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

void ParticleTypeSphere::unpack_comm_vel(int n, int first, double *buf)
{
	int i,m,last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		x[i][0] = buf[m++];
		x[i][1] = buf[m++];
		x[i][2] = buf[m++];
		v[i][0] = buf[m++];
		v[i][1] = buf[m++];
		v[i][2] = buf[m++];
		omega[i][0] = buf[m++];
		omega[i][1] = buf[m++];
		omega[i][2] = buf[m++];
	}
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so parallel::exchange() can test on them
------------------------------------------------------------------------- */

int ParticleTypeSphere::pack_exchange(int i, double *buf)
{
	int m = 1;

	buf[m++] = x[i][0];
	buf[m++] = x[i][1];
	buf[m++] = x[i][2];
	buf[m++] = v[i][0];
	buf[m++] = v[i][1];
	buf[m++] = v[i][2];
	buf[m++] = tag[i];
	buf[m++] = type[i];
	buf[m++] = mask[i];

	buf[m++] = radius[i];
	buf[m++] = rmass[i];
	buf[m++] = omega[i][0];
	buf[m++] = omega[i][1];
	buf[m++] = omega[i][2];
	
	buf[0] = m;
	return m;

}

/* ----------------------------------------------------------------------
   unpack data for atom I for sending to another proc
   xyz must be 1st 3 values, so parallel::exchange() can test on them
------------------------------------------------------------------------- */

int ParticleTypeSphere::unpack_exchange(double *buf)
{
	int nlocal = particle->nlocal;
	if (nlocal == nmax) grow(0);

	int m = 1;
	x[nlocal][0] = buf[m++];
	x[nlocal][1] = buf[m++];
	x[nlocal][2] = buf[m++];
	v[nlocal][0] = buf[m++];
	v[nlocal][1] = buf[m++];
	v[nlocal][2] = buf[m++];
	tag[nlocal] = static_cast<int> (buf[m++]);
	type[nlocal] = static_cast<int> (buf[m++]);
	mask[nlocal] = static_cast<int> (buf[m++]);

	radius[nlocal] = buf[m++];
	rmass[nlocal] = buf[m++];
	omega[nlocal][0] = buf[m++];
	omega[nlocal][1] = buf[m++];
	omega[nlocal][2] = buf[m++];

	particle->nlocal++;
	return m;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void ParticleTypeSphere::copyI2J(int i, int j, int delflag)
{
	tag[j] = tag[i];
	type[j] = type[i];
	mask[j] = mask[i];
	x[j][0] = x[i][0];
	x[j][1] = x[i][1];
	x[j][2] = x[i][2];
	v[j][0] = v[i][0];
	v[j][1] = v[i][1];
	v[j][2] = v[i][2];
	
	radius[j] = radius[i];
	rmass[j] = rmass[i];
	omega[j][0] = omega[i][0];
	omega[j][1] = omega[i][1];
	omega[j][2] = omega[i][2];
}

/* ---------------------------------------------------------------------- */

int ParticleTypeSphere::pack_border(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
	int i,j,m;
	double dx,dy,dz;

	m = 0;
	if (pbc_flag == 0) {
		for (i = 0; i < n; i++) {
		  j = list[i];
		  buf[m++] = x[j][0];
		  buf[m++] = x[j][1];
		  buf[m++] = x[j][2];
		  buf[m++] = tag[j];
		  buf[m++] = type[j];
		  buf[m++] = mask[j];
		  buf[m++] = radius[j];
		  buf[m++] = rmass[j];
		}
	} 
	else {
		dx = pbc[0]*domain->xle;
		dy = pbc[1]*domain->yle;
		dz = pbc[2]*domain->zle;
		for (i = 0; i < n; i++) {
		  j = list[i];
		  buf[m++] = x[j][0] + dx;
		  buf[m++] = x[j][1] + dy;
		  buf[m++] = x[j][2] + dz;
		  buf[m++] = tag[j];
		  buf[m++] = type[j];
		  buf[m++] = mask[j];
		  buf[m++] = radius[j];
		  buf[m++] = rmass[j];
		}
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int ParticleTypeSphere::pack_border_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
	int i,j,m;
	double dx,dy,dz,dvx,dvy,dvz;

	m = 0;
	if (pbc_flag == 0) {
		for (i = 0; i < n; i++) {
		  j = list[i];
		  buf[m++] = x[j][0];
		  buf[m++] = x[j][1];
		  buf[m++] = x[j][2];
		  buf[m++] = tag[j];
		  buf[m++] = type[j];
		  buf[m++] = mask[j];
		  buf[m++] = radius[j];
		  buf[m++] = rmass[j];
		  buf[m++] = v[j][0];
		  buf[m++] = v[j][1];
		  buf[m++] = v[j][2];
		  buf[m++] = omega[j][0];
		  buf[m++] = omega[j][1];
		  buf[m++] = omega[j][2];
		}
	} else {
	  dx = pbc[0]*domain->xle;
	  dy = pbc[1]*domain->yle;
	  dz = pbc[2]*domain->zle;
		if (!deform_vremap) {
			for (i = 0; i < n; i++) {
				j = list[i];
				buf[m++] = x[j][0] + dx;
				buf[m++] = x[j][1] + dy;
				buf[m++] = x[j][2] + dz;
				buf[m++] = tag[j];
				buf[m++] = type[j];
				buf[m++] = mask[j];
				buf[m++] = radius[j];
				buf[m++] = rmass[j];
				buf[m++] = v[j][0];
				buf[m++] = v[j][1];
				buf[m++] = v[j][2];
				buf[m++] = omega[j][0];
				buf[m++] = omega[j][1];
				buf[m++] = omega[j][2];
			}
		} 
		// the following part is still under investigation
		else {
			dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
			dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
			dvz = pbc[2]*h_rate[2];
			for (i = 0; i < n; i++) {
				j = list[i];
				buf[m++] = x[j][0] + dx;
				buf[m++] = x[j][1] + dy;
				buf[m++] = x[j][2] + dz;
				buf[m++] = tag[j];
				buf[m++] = type[j];
				buf[m++] = mask[j];
				buf[m++] = rmass[j];
				if (mask[i] & deform_groupbit) {
					buf[m++] = v[j][0] + dvx;
					buf[m++] = v[j][1] + dvy;
					buf[m++] = v[j][2] + dvz;
				} else {
					buf[m++] = v[j][0];
					buf[m++] = v[j][1];
					buf[m++] = v[j][2];
					buf[m++] = omega[j][0];
					buf[m++] = omega[j][1];
					buf[m++] = omega[j][2];
				}
			}
		} // if (!deform_vremap)
	} // pbc_flag == 0
	return m;
}

/* ---------------------------------------------------------------------- */

void ParticleTypeSphere::unpack_border(int n, int first, double *buf)
{
	int i,m,last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		if (i == nmax) grow(0);
		x[i][0] = buf[m++];
		x[i][1] = buf[m++];
		x[i][2] = buf[m++];
		tag[i] = static_cast<int> (buf[m++]);
		type[i] = static_cast<int> (buf[m++]);
		mask[i] = static_cast<int> (buf[m++]);
		radius[i] = buf[m++];
		rmass[i] = buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

void ParticleTypeSphere::unpack_border_vel(int n, int first, double *buf)
{
	int i,m,last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		if (i == nmax) grow(0);
		x[i][0] = buf[m++];
		x[i][1] = buf[m++];
		x[i][2] = buf[m++];
		tag[i] = static_cast<int> (buf[m++]);
		type[i] = static_cast<int> (buf[m++]);
		mask[i] = static_cast<int> (buf[m++]);
		radius[i] = buf[m++];
		rmass[i] = buf[m++];
		v[i][0] = buf[m++];
		v[i][1] = buf[m++];
		v[i][2] = buf[m++];
		omega[i][0] = buf[m++];
		omega[i][1] = buf[m++];
		omega[i][2] = buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

int ParticleTypeSphere::pack_reverse(int n, int first, double *buf)
{
	int i,m,last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		buf[m++] = f[i][0];
		buf[m++] = f[i][1];
		buf[m++] = f[i][2];
		buf[m++] = torque[i][0];
		buf[m++] = torque[i][1];
		buf[m++] = torque[i][2];
	}
	return m;
}

/* ---------------------------------------------------------------------- */

void ParticleTypeSphere::unpack_reverse(int n, int *list, double *buf)
{
	int i,j,m;

	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		f[j][0] += buf[m++];
		f[j][1] += buf[m++];
		f[j][2] += buf[m++];
		torque[j][0] += buf[m++];
		torque[j][1] += buf[m++];
		torque[j][2] += buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

bigint ParticleTypeSphere::memory_usage()
{

	return 0;
}

/* ---------------------------------------------------------------------- */

void ParticleTypeSphere::unpack_force(int n, int first, double *buf)
{
	int i, m, last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		f[i][0] = buf[m++];
		f[i][1] = buf[m++];
		f[i][2] = buf[m++];

	}
}

/* ---------------------------------------------------------------------- */

int ParticleTypeSphere::pack_force(int n, int *list, double *buf,
	int pbc_flag, int *pbc)
{
	int i, j, m;
	double dx, dy, dz;

	m = 0;

	for (i = 0; i < n; i++) {
		j = list[i];
		buf[m++] = f[j][0];
		buf[m++] = f[j][1];
		buf[m++] = f[j][2];
	}


	return m;
}