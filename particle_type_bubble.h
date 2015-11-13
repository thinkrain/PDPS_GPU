/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifdef PARTICLE_CLASS

ParticleStyle(bubble, ParticleTypeBubble)

#else

#ifndef PS_PARTICLE_TYPE_BUBBLE_H
#define PS_PARTICLE_TYPE_BUBBLE_H

#include "particle_type.h"

namespace PDPS_NS {

class ParticleTypeBubble : public ParticleType {
public:

	ParticleTypeBubble(class PDPS *, int, char **);
	~ParticleTypeBubble();

	 void grow(int);
	
	void create_particle(int, double *);

	void data_particle(double *, char **);
	void data_vel(int, char **);

	int pack_comm(int, int*, double *, int, int *);
	int pack_comm_vel(int, int *, double *, int, int *);
	void unpack_comm(int, int, double *);
	void unpack_comm_vel(int, int, double *);

    int pack_exchange(int, double *);
	int unpack_exchange(double *);
	void copyI2J(int, int, int);

	int pack_border(int, int *, double *, int, int *);
    int pack_border_vel(int, int *, double *, int, int *);
	void unpack_border(int, int, double *);
    void unpack_border_vel(int, int, double *);

	int pack_reverse(int, int, double *);
	void unpack_reverse(int, int *, double *);

	bigint memory_usage();

private:
	int *tag, *type, *mask, *image;
	double **x, **v, **f;
	double *rho;


};

}

#endif
#endif
