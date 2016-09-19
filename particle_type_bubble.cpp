/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "domain.h"
#include "error.h"
#include "memory.h"
#include "parallel.h"
#include "particle.h"
#include "particle_type_bubble.h"

using namespace PDPS_NS;

#define DELTA 10000  
#define EPSILON 1.0e-6
#define BIG 1e20

/* ---------------------------------------------------------------------- */

ParticleTypeBubble::ParticleTypeBubble(PDPS *ps, int narg, char **arg) : ParticleType(ps, narg, arg)
{
	 comm_x_only = 0; // we communicate not only x forward but also vest ...
	 comm_f_only = 0; // we also communicate de and drho in reverse direction
	 size_forward = 4; // 3 + rho + e + vest[3], that means we may only communicate 5 in hybrid
	 size_reverse = 3; // 3 + drho + de
	 size_border = 7; // 6 + rho + e + vest[3] + cv
	 size_velocity = 3;
 //    particle->size_data_atom = 7;
	 particle->size_data_vel = 3;
	 particle->xcol_data = 6;

	 particle->ee_flag = 1;
	 particle->rho_flag = 1;
	 particle->cv_flag = 1;
	 particle->vest_flag = 1;
}

/* ---------------------------------------------------------------------- */

ParticleTypeBubble::~ParticleTypeBubble()
{

}

/* ----------------------------------------------------------------------
   Grow particle arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void ParticleTypeBubble::grow(int n)
{
	if (n == 0) nmax += DELTA;
	else nmax = n;
	particle->nmax = nmax;

	// nmax < 0: overflow
	if (nmax < 0 || nmax > BIG) {
		error->one(FLERR,"Too many particles to be created for this system");
	}

	tag = memory->grow(particle->tag, nmax, "particle: tag");
	type = memory->grow(particle->type, nmax, "Particle: type");
	mask = memory->grow(particle->mask, nmax, "Particle: mask");
//	image = memory->grow(particle->image, nmax, "particle:image");
	x = memory->grow(particle->x, nmax, 3, "particle: x");
	v = memory->grow(particle->v, nmax, 3, "particle: v");
	f = memory->grow(particle->f, nmax, 3, "particle: f");

	rho = memory->grow(particle->rho, nmax, "particle:rho");


}

/* ----------------------------------------------------------------------
Create particle
------------------------------------------------------------------------- */

void ParticleTypeBubble::create_particle(int itype, double *coord)
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
	}
	rho[nlocal] = 0.0;
	particle->nlocal++;
}

/* ----------------------------------------------------------------------
unpack one line from Atoms section of data file
initialize other atom quantities
------------------------------------------------------------------------- */

void ParticleTypeBubble::data_particle(double *coord, char **values)
{
	int nlocal = particle->nlocal;
	if (nlocal == nmax) grow(0);

	tag[nlocal] = atoi(values[0]);
	if (tag[nlocal] <= 0)
		error->one(FLERR, "Invalid atom ID in Atoms section of data file");

	type[nlocal] = atoi(values[1]);
	if (type[nlocal] <= 0 || type[nlocal] > particle->ntypes)
		error->one(FLERR, "Invalid atom type in Atoms section of data file");

	rho[nlocal] = atof(values[2]);


	x[nlocal][0] = coord[0];
	x[nlocal][1] = coord[1];
	x[nlocal][2] = coord[2];

	//image[nlocal] = imagetmp;

	mask[nlocal] = 1;
	v[nlocal][0] = 0.0;
	v[nlocal][1] = 0.0;
	v[nlocal][2] = 0.0;


	particle->nlocal++;
}

/* ----------------------------------------------------------------------
unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void ParticleTypeBubble::data_vel(int m, char **values)
{
	double **v = particle->v;
	v[m][0] = atof(values[0]);
	v[m][1] = atof(values[1]);
	v[m][2] = atof(values[2]);
}

/* ---------------------------------------------------------------------- */

int ParticleTypeBubble::pack_comm(int n, int *list, double *buf,
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
			buf[m++] = rho[j];

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
			buf[m++] = rho[j];

		}
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int ParticleTypeBubble::pack_comm_vel(int n, int *list, double *buf,
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
			buf[m++] = rho[j];
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
				buf[m++] = rho[j];


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
				buf[m++] = rho[j];

			}
		} // else if (deform_vremap == 1)
	} // else if (pbc_flag != 0)
  return m;
}

/* ---------------------------------------------------------------------- */

void ParticleTypeBubble::unpack_comm(int n, int first, double *buf)
{
	int i,m,last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		x[i][0] = buf[m++];
		x[i][1] = buf[m++];
		x[i][2] = buf[m++];
		rho[i] = buf[m++];

	}
}

/* ---------------------------------------------------------------------- */

void ParticleTypeBubble::unpack_comm_vel(int n, int first, double *buf)
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
		rho[i] = buf[m++];

	}
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so parallel::exchange() can test on them
------------------------------------------------------------------------- */

int ParticleTypeBubble::pack_exchange(int i, double *buf)
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
//	buf[m++] = image[i];
	buf[m++] = rho[i];


	
	buf[0] = m;
	return m;

}

/* ----------------------------------------------------------------------
   unpack data for atom I for sending to another proc
   xyz must be 1st 3 values, so parallel::exchange() can test on them
------------------------------------------------------------------------- */

int ParticleTypeBubble::unpack_exchange(double *buf)
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
//	image[nlocal] =static_cast<int> (buf[m++]);
	rho[nlocal] = buf[m++];


	particle->nlocal++;
	return m;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void ParticleTypeBubble::copyI2J(int i, int j, int delflag)
{
	tag[j] = tag[i];
	type[j] = type[i];
	mask[j] = mask[i];
//	image[j] = image[i];
	x[j][0] = x[i][0];
	x[j][1] = x[i][1];
	x[j][2] = x[i][2];
	v[j][0] = v[i][0];
	v[j][1] = v[i][1];
	v[j][2] = v[i][2];
	rho[j] = rho[i];

}

/* ---------------------------------------------------------------------- */

int ParticleTypeBubble::pack_border(int n, int *list, double *buf,
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
		  buf[m++] = rho[j];

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
		  buf[m++] = rho[j];

		}
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int ParticleTypeBubble::pack_border_vel(int n, int *list, double *buf,
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
		  buf[m++] = v[j][0];
		  buf[m++] = v[j][1];
		  buf[m++] = v[j][2];
		  buf[m++] = rho[j];

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
				buf[m++] = v[j][0];
				buf[m++] = v[j][1];
				buf[m++] = v[j][2];
				buf[m++] = rho[j];
			}
		} else {
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
				if (mask[i] & deform_groupbit) {
					buf[m++] = v[j][0] + dvx;
					buf[m++] = v[j][1] + dvy;
					buf[m++] = v[j][2] + dvz;
				} else {
					buf[m++] = v[j][0];
					buf[m++] = v[j][1];
					buf[m++] = v[j][2];
				}
				buf[m++] = rho[j];

			}
		} // if (!deform_vremap)
	} // pbc_flag == 0
	return m;
}

/* ---------------------------------------------------------------------- */

void ParticleTypeBubble::unpack_border(int n, int first, double *buf)
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
		rho[i] = buf[m++];

	}
}

/* ---------------------------------------------------------------------- */

void ParticleTypeBubble::unpack_border_vel(int n, int first, double *buf)
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
		v[i][0] = buf[m++];
		v[i][1] = buf[m++];
		v[i][2] = buf[m++];
		rho[i] = buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

int ParticleTypeBubble::pack_reverse(int n, int first, double *buf)
{
	int i,m,last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		buf[m++] = f[i][0];
		buf[m++] = f[i][1];
		buf[m++] = f[i][2];

	}
	return m;
}

/* ---------------------------------------------------------------------- */

void ParticleTypeBubble::unpack_reverse(int n, int *list, double *buf)
{
	int i,j,m;

	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		f[j][0] += buf[m++];
		f[j][1] += buf[m++];
		f[j][2] += buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

bigint ParticleTypeBubble::memory_usage()
{

	return 0;
}

/* ---------------------------------------------------------------------- */

void ParticleTypeBubble::unpack_force(int n, int first, double *buf)
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

int ParticleTypeBubble::pack_force(int n, int *list, double *buf,
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