/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_Create_Particle_H
#define PS_Create_Particle_H

#include "pointers.h"

namespace PDPS_NS {

class CreateParticle : protected Pointers {
public:
	int nparticles_previous, nparticles_now;     // number of particles created in one command
	int npx, npy, npz;                           // number of particle along each edge of the corresponding region
	double xlo, xhi, ylo, yhi, zlo, zhi;         // the lower bound of the coresponding region
	double xle, yle, zle;
	double extent_xlo, extent_xhi, extent_ylo, extent_yhi, extent_zlo, extent_zhi;
	double extent_xle, extent_yle, extent_zle;
	double dx, dy, dz;                           // lattice spacing
	double sublo[3], subhi[3];
	
	int id, localid;
	int tid, rid;                                // type index and region index

	double rho;                                  // number of density

	// create particle with radius distribution
	double rsingle, rlo, rhi, rmean, rsigma;
	int dist_style;

	CreateParticle(class PDPS *);
	void command(int, char **);
	void create_number_density(int,int);
	void create_single(double, double, double);
	void create_spacing(int, int);
	void create_random();

private:
	int n_target;
	int seed;
	double spacing;                    // the smallest radius of particle to be created
	class Lattice *lattice;
	class RanPark *random_radius;

	int random_no_overlap_flag;

	void find_factors2(int,int *);
	void find_factors3(int,int *);

	double get_uniform_radius();
	double get_gaussian_radius();
};

}

#endif
