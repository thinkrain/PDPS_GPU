/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_PARTICLE_H
#define PS_PARTICLE_H

#include "pointers.h"

namespace PDPS_NS {

class Particle : protected Pointers {
public:
	bigint nparticles;               // total number of particles
	int nlocal, nghost;              // # of owned and ghost atoms on this processor
	int nmax;                        // maxmum # of owned+ghost in arrays on this proc
	int nfirst;
	int ntypes;                      // number of types of particles
	int nmats;   

	double *mass;

    int *tag, *type, *mask;          // global ID, type of particle, mask for group
	double **x, **v, **f;
	double *radius, *rmass, *density, *poro, *volume, *hlocal;
	double **omega, **torque;

		// USER-SPH package
	double *rho,*drho,*e,*de,*cv;
	double **vest;
	int ee_flag, rho_flag, cv_flag, vest_flag;

	// flag for different particle style
	int atomic_flag, sphere_flag;
	int rmass_flag, radius_flag, omega_flag, torque_flag;
	int density_flag, e_flag;

	int tag_enable; 

	// used in read_data command
	int size_data_atom;
	int size_data_vel; 
	int xcol_data;                   // column (1-N) where x is in Atom line

	class ParticleType *ptype;
	char *particle_style;

	// map  
	int map_style;                   // map style
									 // 0 = none; 1 = arrary; 
	
	Particle(class PDPS *);
	~Particle();

	void init();
	
	void create_particle_type(const char *, int, char **);
	
	void grow(int);
	
	
	//void create_particle(int, double *);
	void add_tag();
	void set_pos(int, char **);
	void set_mass(const char *);
	void set_mass(int, char **);
	void set_density(int, char **);
	void set_radius(int, char **);
	void set_rho(int, char **);
	void set_energy(int, char **);
	void allocate_type_arrays();

	void data_particles(int, char *);
	void data_vels(int, char *);

	void delete_particle(int);
	void save_particle(int narg, char** arg);

	int count_words(const char *);


	// map
	inline int map(int global) {
		if (map_style == 1) return map_array[global];
		else return map_find_hash(global);
	};

	void map_init();
	void map_clear();
	void map_set();
	void map_one(int, int);
	void map_delete();
	int map_find_hash(int);

	void lost_check();
	int map_tag_max;

	//	variable in GPU
	int *devTag, *devType, *devMask;
	double *devMass;
	double *devCoordX, *devCoordY, *devCoordZ;
	double	*devForceX, *devForceY, *devForceZ;
	double	*devVeloX, *devVeloY, *devVeloZ;
	double	*devVestX, *devVestY, *devVestZ;
	double  *devRho, *devRadius, *devRmass, *devDensity, *devPoro, *devVolume;

	// pointer to pre-allocated device buffer
	double	*devHostCoord;
	double	*devHostVelo;
	double  *devHostVest;
	double	*devHostForce;
	double  *devHostRho;
	int *devHostTag;
	int *devHostType;
	int *devHostMask;
	double	*devHostMassType;
	double *devHostRadius, *devHostRmass, *devHostDensity, *devHostPoro, *devHostVolume;

	// pointer to host memory, e.g. atom->x
	double	*ptrHostCoord;
	double	*ptrHostVelo;
	double  *ptrHostVest;
	double	*ptrHostForce;
	double  *ptrHostRho;
	int *ptrHostTag;
	int *ptrHostType;
	int *ptrHostMask;
	double	*ptrHostMassType;
	double *ptrHostRadius, *ptrHostRmass, *ptrHostDensity, *ptrHostPoro, *ptrHostVolume;

	//	GPU pin memory
	void PinHostArray();
	void UnpinHostArray();
	void TransferC2G();
	void TransferG2C();

	// transfer data from CPU to GPU

private:
	int procid;

	
	int *map_array;

	struct HashElem {
		int global;                 // key to search on = global ID
		int local;                  // value associated with key = local index
		int next;                   // next entry in this bucket, -1 if last
	};
	int map_nhash;                  // # of entries hash table can hold
	int map_nused;                  // # of actual entries in hash table
	int map_free;                   // ptr to 1st unused entry in hash table
	int map_nbucket;                // # of hash buckets
	int *map_bucket;                // ptr to 1st entry in each bucket
	HashElem *map_hash;             // hash table
	int *primes;                    // table of prime #s for hashing
	int nprimes;                    // # of primes

	
};

}

#endif

