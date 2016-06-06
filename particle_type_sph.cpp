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
//#include "particle_type.h"
#include "particle_type_sph.h"

using namespace PDPS_NS;

#define DELTA 10000  
#define EPSILON 1.0e-6
#define BIG 1e20

/* ---------------------------------------------------------------------- */

ParticleTypeSph::ParticleTypeSph(PDPS *ps, int narg, char **arg) : ParticleType(ps, narg, arg)
{
	 comm_x_only = 0; // we communicate not only x forward but also vest ...
	 comm_f_only = 0; // we also communicate de and drho in reverse direction
	 size_forward = 14; // 3 + rho + e + vest[3] + density + radius + rmass + poro + volume + hlocal, that means we may only communicate 5 in hybrid
	 size_reverse = 5; // 3 + drho + de
	 size_border = 18; // 6 + rho + e + vest[3] + cv + radius + rmass + density + poro + volume + hlocal(for bubble)
	 size_velocity = 3;
  //   particle->size_data_atom = 8;
//	 particle->size_data_vel = 4;
//	 particle->xcol_data = 6;

	 particle->ee_flag = 1;
	 particle->rho_flag = 1;
	 particle->cv_flag = 1;
	 particle->vest_flag = 1;

	 particle->atomic_flag = 1; 
	 particle->sphere_flag = 1;
// 	 particle->radius_flag = 1;
//	 particle->rmass_flag = 1;
//	 particle->density_flag = 1;

}

/* ---------------------------------------------------------------------- */

ParticleTypeSph::~ParticleTypeSph()
{

}

/* ----------------------------------------------------------------------
   Grow particle arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void ParticleTypeSph::grow(int n)
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
	drho = memory->grow(particle->drho, nmax, "particle:drho");
	e = memory->grow(particle->e, nmax, "particle:e");
	de = memory->grow(particle->de, nmax, "particle:de");
	cv = memory->grow(particle->cv, nmax, "particle:cv");
	vest = memory->grow(particle->vest, nmax, 3, "particle:vest");
	density = memory->grow(particle->density, nmax, "particle: density");
	radius = memory->grow(particle->radius, nmax, "particle: radius");
	rmass = memory->grow(particle->rmass, nmax, "particle: rmass");
	poro = memory->grow(particle->poro, nmax, "particle: poro");
	volume = memory->grow(particle->volume, nmax, "particle: volume");
	hlocal = memory->grow(particle->hlocal, nmax, "particle: hlocal");

}

/* ----------------------------------------------------------------------
Create particle
------------------------------------------------------------------------- */

void ParticleTypeSph::create_particle(int itype, double *coord)
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
		vest[nlocal][i] = 0.0;
	}
	rho[nlocal] = 0.0;
	e[nlocal] = 0.0;
	cv[nlocal] = 0.0;
	de[nlocal] = 0.0;
	drho[nlocal] = 0.0;
	radius[nlocal] = 0.0;
	density[nlocal] = 0.0;
	rmass[nlocal] = 0.0;
	poro[nlocal] = 1.0;
	volume[nlocal] = 0.0;
	hlocal[nlocal] = 0.0;
	particle->nlocal++;
}

/* ----------------------------------------------------------------------
unpack one line from Atoms section of data file
initialize other atom quantities
------------------------------------------------------------------------- */

void ParticleTypeSph::data_particle(double *coord, char **values)
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
	e[nlocal] = atof(values[3]);
	cv[nlocal] = atof(values[4]);

	x[nlocal][0] = coord[0];
	x[nlocal][1] = coord[1];
	x[nlocal][2] = coord[2];

	//image[nlocal] = imagetmp;

	mask[nlocal] = 1;
	v[nlocal][0] = 0.0;
	v[nlocal][1] = 0.0;
	v[nlocal][2] = 0.0;

	vest[nlocal][0] = 0.0;
	vest[nlocal][1] = 0.0;
	vest[nlocal][2] = 0.0;

	de[nlocal] = 0.0;
	drho[nlocal] = 0.0;

	radius[nlocal] = 0.0;
	density[nlocal] = 0.0;
	rmass[nlocal] = 0.0;

	poro[nlocal] = 0.0;
	volume[nlocal] = 0.0;
	hlocal[nlocal] = 0.0;

	particle->nlocal++;
}

/* ----------------------------------------------------------------------
unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void ParticleTypeSph::data_vel(int m, char **values)
{
	double **v = particle->v;
	v[m][0] = atof(values[0]);
	v[m][1] = atof(values[1]);
	v[m][2] = atof(values[2]);
}

/* ---------------------------------------------------------------------- */

int ParticleTypeSph::pack_comm(int n, int *list, double *buf,
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
			buf[m++] = e[j];
			buf[m++] = vest[j][0];
			buf[m++] = vest[j][1];
			buf[m++] = vest[j][2];
			buf[m++] = radius[j];
			buf[m++] = density[j];
			buf[m++] = rmass[j];
			buf[m++] = poro[j];
			buf[m++] = volume[j];
			buf[m++] = hlocal[j];
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
			buf[m++] = e[j];
			buf[m++] = vest[j][0];
			buf[m++] = vest[j][1];
			buf[m++] = vest[j][2];
			buf[m++] = radius[j];
			buf[m++] = density[j];
			buf[m++] = rmass[j];
			buf[m++] = poro[j];
			buf[m++] = volume[j];
			buf[m++] = hlocal[j];
		}
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int ParticleTypeSph::pack_comm_vel(int n, int *list, double *buf,
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
			buf[m++] = e[j];
			buf[m++] = vest[j][0];
			buf[m++] = vest[j][1];
			buf[m++] = vest[j][2];
			buf[m++] = radius[j];
			buf[m++] = density[j];
			buf[m++] = rmass[j];
			buf[m++] = poro[j];
			buf[m++] = volume[j];
			buf[m++] = hlocal[j];
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
				buf[m++] = e[j];
				buf[m++] = vest[j][0];
				buf[m++] = vest[j][1];
				buf[m++] = vest[j][2];
				buf[m++] = radius[j];
				buf[m++] = density[j];
				buf[m++] = rmass[j];
				buf[m++] = poro[j];
				buf[m++] = volume[j];
				buf[m++] = hlocal[j];

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
				buf[m++] = e[j];
				buf[m++] = vest[j][0];
				buf[m++] = vest[j][1];
				buf[m++] = vest[j][2];
				buf[m++] = radius[j];
				buf[m++] = density[j];
				buf[m++] = rmass[j];
				buf[m++] = poro[j];
				buf[m++] = volume[j];
				buf[m++] = hlocal[j];
			}
		} // else if (deform_vremap == 1)
	} // else if (pbc_flag != 0)
  return m;
}

/* ---------------------------------------------------------------------- */

void ParticleTypeSph::unpack_comm(int n, int first, double *buf)
{
	int i,m,last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		x[i][0] = buf[m++];
		x[i][1] = buf[m++];
		x[i][2] = buf[m++];
		rho[i] = buf[m++];
		e[i] = buf[m++];
		vest[i][0] = buf[m++];
		vest[i][1] = buf[m++];
		vest[i][2] = buf[m++];
		radius[i] = buf[m++];
		density[i] = buf[m++];
		rmass[i] = buf[m++];
		poro[i] = buf[m++];
		volume[i] = buf[m++];
		hlocal[i] = buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

void ParticleTypeSph::unpack_comm_vel(int n, int first, double *buf)
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
		e[i] = buf[m++];
		vest[i][0] = buf[m++];
		vest[i][1] = buf[m++];
		vest[i][2] = buf[m++];
		radius[i] = buf[m++];
		density[i] = buf[m++];
		rmass[i] = buf[m++];
		poro[i] = buf[m++];
		volume[i] = buf[m++];
		hlocal[i] = buf[m++];
	
	}
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so parallel::exchange() can test on them
------------------------------------------------------------------------- */

int ParticleTypeSph::pack_exchange(int i, double *buf)
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
	buf[m++] = rho[i];
	buf[m++] = e[i];
	buf[m++] = cv[i];
	buf[m++] = vest[i][0];
	buf[m++] = vest[i][1];
	buf[m++] = vest[i][2];
	buf[m++] = density[i];
	buf[m++] = poro[i];
	buf[m++] = volume[i];
	buf[m++] = hlocal[i];

	
	buf[0] = m;
	return m;

}

/* ----------------------------------------------------------------------
   unpack data for atom I for sending to another proc
   xyz must be 1st 3 values, so parallel::exchange() can test on them
------------------------------------------------------------------------- */

int ParticleTypeSph::unpack_exchange(double *buf)
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
	rho[nlocal] = buf[m++];
	e[nlocal] = buf[m++];
	cv[nlocal] = buf[m++];
	vest[nlocal][0] = buf[m++];
	vest[nlocal][1] = buf[m++];
	vest[nlocal][2] = buf[m++];
	density[nlocal] = buf[m++];
	poro[nlocal] = buf[m++];
	volume[nlocal] = buf[m++];
	hlocal[nlocal] = buf[m++];

	particle->nlocal++;
	return m;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void ParticleTypeSph::copyI2J(int i, int j, int delflag)
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
	radius[j] = radius[i];
	rmass[j] = rmass[i];
	rho[j] = rho[i];
	drho[j] = drho[i];
//	density[j] = density[i];
	e[j] = e[i];
	de[j] = de[i];
	cv[j] = cv[i];
	vest[j][0] = vest[i][0];
	vest[j][1] = vest[i][1];
	vest[j][2] = vest[i][2];
	density[j] = density[i];
	poro[j] = poro[i];
	volume[j] = volume[i];
	hlocal[j] = hlocal[i];
}

/* ---------------------------------------------------------------------- */

int ParticleTypeSph::pack_border(int n, int *list, double *buf,
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
		  buf[m++] = rho[j];
		  buf[m++] = e[j];
          buf[m++] = cv[j];
          buf[m++] = vest[j][0];
          buf[m++] = vest[j][1];
          buf[m++] = vest[j][2];
		  buf[m++] = density[i];
		  buf[m++] = poro[j];
		  buf[m++] = volume[j];
		  buf[m++] = hlocal[i];
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
		  buf[m++] = rho[j];
		  buf[m++] = e[j];
          buf[m++] = cv[j];
          buf[m++] = vest[j][0];
          buf[m++] = vest[j][1];
          buf[m++] = vest[j][2];
		  buf[m++] = density[i];
		  buf[m++] = poro[j];
		  buf[m++] = volume[j];
		  buf[m++] = hlocal[i];

		}
	}
	return m;
}

/* ---------------------------------------------------------------------- */

int ParticleTypeSph::pack_border_vel(int n, int *list, double *buf,
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
		  buf[m++] = rho[j];
		  buf[m++] = e[j];
	 	  buf[m++] = cv[j];
          buf[m++] = vest[j][0];
		  buf[m++] = vest[j][1];
		  buf[m++] = vest[j][2];
		  buf[m++] = density[j];
		  buf[m++] = poro[j];
		  buf[m++] = volume[j];
		  buf[m++] = hlocal[i];
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
				buf[m++] = rho[j];
			    buf[m++] = e[j];
	 		    buf[m++] = cv[j];
				buf[m++] = vest[j][0];
				buf[m++] = vest[j][1];
				buf[m++] = vest[j][2];
				buf[m++] = density[j];
				buf[m++] = poro[j];
				buf[m++] = volume[j];
				buf[m++] = hlocal[i];
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
				buf[m++] = radius[j];
				buf[m++] = rmass[j];
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
			    buf[m++] = e[j];
	 		    buf[m++] = cv[j];
				buf[m++] = vest[j][0];
				buf[m++] = vest[j][1];
				buf[m++] = vest[j][2];
				buf[m++] = density[j];
				buf[m++] = poro[j];
				buf[m++] = volume[j];
				buf[m++] = hlocal[i];
			}
		} // if (!deform_vremap)
	} // pbc_flag == 0
	return m;
}

/* ---------------------------------------------------------------------- */

void ParticleTypeSph::unpack_border(int n, int first, double *buf)
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
		rho[i] = buf[m++];
		e[i] = buf[m++];
		cv[i] = buf[m++];
		vest[i][0] = buf[m++];
		vest[i][1] = buf[m++];
		vest[i][2] = buf[m++];
		density[i] = buf[m++];
		poro[i] = buf[m++];
		volume[i] = buf[m++];
		hlocal[i] = buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

void ParticleTypeSph::unpack_border_vel(int n, int first, double *buf)
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
		rho[i] = buf[m++];
		e[i] = buf[m++];
		cv[i] = buf[m++];
		vest[i][0] = buf[m++];
		vest[i][1] = buf[m++];
		vest[i][2] = buf[m++];
		density[i] = buf[m++];
		poro[i] = buf[m++];
		volume[i] = buf[m++];
		hlocal[i] = buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

int ParticleTypeSph::pack_reverse(int n, int first, double *buf)
{
	int i,m,last;

	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		buf[m++] = f[i][0];
		buf[m++] = f[i][1];
		buf[m++] = f[i][2];
		buf[m++] = drho[i];
		buf[m++] = de[i];
	}
	return m;
}

/* ---------------------------------------------------------------------- */

void ParticleTypeSph::unpack_reverse(int n, int *list, double *buf)
{
	int i,j,m;

	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		f[j][0] += buf[m++];
		f[j][1] += buf[m++];
		f[j][2] += buf[m++];
		drho[j] += buf[m++];
		de[j] += buf[m++];
	}
}

/* ---------------------------------------------------------------------- */

bigint ParticleTypeSph::memory_usage()
{

	return 0;
}
