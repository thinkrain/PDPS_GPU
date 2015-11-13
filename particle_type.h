/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_PARTICLE_TYPE_H
#define PS_PARTICLE_TYPE_H

#include "pointers.h"

namespace PDPS_NS {

class ParticleType : protected Pointers {
public:

	// parallel
	int comm_x_only;
	int comm_f_only;                  
	int size_forward;                // # of values per atom in parallel
	int size_reverse;                // # in reverse parallel
	int size_border;                 // # in border parallel
	int size_velocity;               // # of velocity based quantities
	
	ParticleType(class PDPS *, int, char **);
	~ParticleType();
	
	virtual void init();
	virtual void grow(int) = 0;
	
	virtual void create_particle(int, double *) {};
	virtual void create_particle(int, double *, double) {};

	virtual int pack_comm(int, int*, double *, int, int *) = 0;
	virtual int pack_comm_vel(int, int *, double *, int, int *) = 0;
	virtual void unpack_comm(int, int, double *) = 0;
	virtual void unpack_comm_vel(int, int, double *) = 0;

	virtual int pack_exchange(int, double *) = 0;
	virtual int unpack_exchange(double *) = 0;
	virtual void copyI2J(int, int, int) = 0;

	virtual int pack_border(int, int *, double *, int, int *) = 0;
    virtual int pack_border_vel(int, int *, double *, int, int *) = 0;
	virtual void unpack_border(int, int, double *) = 0;
    virtual void unpack_border_vel(int, int, double *) = 0;

	virtual int pack_reverse(int, int, double *) = 0;
	virtual void unpack_reverse(int, int *, double *) = 0;

	virtual void data_particle(double *, char **) = 0;
	//virtual void data_particles(int, char *) = 0;
	virtual void data_vel(int, char **) = 0;
	//virtual void data_vels(int, char *) = 0;
	

	virtual bigint memory_usage() = 0;

protected:
	int nmax;                       // local copy of particle->nmax
	int deform_vremap;              // local copy of domain properties
	int deform_groupbit;
	double *h_rate;

};
}

#endif
