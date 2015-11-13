/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_NEIGHBOR_H
#define PS_NEIGHBOR_H

#include "pointers.h"

namespace PDPS_NS {

class Neighbor : protected Pointers {
public:
	//double rcut, rcutsq;                 // cutoff
	double rcut_max, rcut_min;           // rc_max/min: maximum cut-off
	double rskin;                        // skin distance
	double **rneigh, **rneighsq;         // rneigh[][] = rc_max + skin
	double *rneigh_type, *rneighsq_type; // the maximum rneigh for itype particle
	double rneigh_max, rneigh_min;       // rcut_max + rskin; rcut_min + rskin
	double rneighsq_max, rneighsq_min;
	
	//double rneighsq;                   // rneigh^2

	double rij;                          // distance between two particles
	double rijsq;                        // rij^2

	// neigh_list
	class NeighList *neighlist;          // neighbor list
	int neigh_flag;                      // 1 exist

	int pgsize;                          // size of one page
	int onesize;                         // max # of neighbor list allowed for one particle 

	// neighbor command
	int every;                           // build every this many steps
	int delay;                           // delay build for this many steps
	bigint last_build;                   // the timestep when the neig_list is built last time
	int ago;
	int build_before;                    // has the neigh list built before or not
	int nflag;                           // neighbor build flag
	int dist_check;

	int nbuilds, nbuilds_total;          // # of times the neighbor list has been built
	int ndanger, ndanger_total;          // number of dangerous builds

	int style;                           // 0,1 = single, multi

	Neighbor(class PDPS *);
	~Neighbor();

	void init();
	void settings(int, char **);
	void modify_params(int, char **);    // called by neigh_modify command

	void build();
	void allocate(int);
	void setup_cells();
	
	void create_linked_list();
	void create_neigh_list();
	int check_cell_boundary(int []);            // check cell validity based on the boundary conditions
	int check_distance();                       // check max distance moved since last build
	//int cell_upper_right(int [], int []);     // check if cell is upper right to the scanning cell
	//int particle_upper_right(int, int); // check if ghost paritcle in the neighbor list is upper right to the particle itself 

	int decide();
	double find_distance(int, int);             // return distance^2

	void debug();
	  
protected:
    int procid, nprocs;

	// other classes
	int dim;
	int nparticles;
	int nlocal_old, nmax_old;
	double *boxlo, *boxhi, *boxle;
	double *sublo, *subhi, *suble;
	double boxlo_old[3], boxhi_old[3];
	double sublo_old[3], subhi_old[3];
	double ghost_sublo[3], ghost_subhi[3];

	// linked_list cell method

	// global cells related
	double cle[3];                   // global cell length (without considering ghost particles)
	int nc[3];                       // global # of cells along one direction
	int ncxy, ncxyz;                 // global # of particle in the x-y plane and the whole simulation box

	// local cells related
	// considering local subdomain extend by the ghost particles: 
	// lower bound = sublo - rcghost; upper bound = subhi + rcghost

	int subclo[3];                    // lowerst global index in the extended subdomain      
	int subnc[3];
	int subncxy, subncxyz;
	int subncxyz_max;

	// single cell offsets related 
	int cox, coy, coz;                 // cell-offset along x, y, and z
	int co[3];                         // cell-offsets array
	int *coffsets;                     // cell offsets;
	int noffsets;                      // # of of real cell offsets
	int max_noffsets;                  // max # of cell offsets

	// multi cell offsets 
	int **coffsets_multi;
	double **distsq_multi; 
	int *noffsets_multi;
	int max_noffsets_multi;

	// linked_list related
	int *head;                         // store the header of linked_list
	int *linked_list;                  // store linked_list

	// neighbor list realated 
	double triggersq;
	double **x_old;
	int boxcheck;  

	int coord2cell(double *);                      // mapping particle coords to a cell index
	int coord2cell(double *, int &, int &, int &); // the same as above, but return the cell location also
	double cell_distance(int, int, int);           // distance between binx

	
	//void choose_build();

	// pairwise build functions
	typedef void (Neighbor::*PairPtr)(class NeighList *);
	PairPtr pair_build;

	void half_linkedlist(class NeighList *);
	void half_multi(class NeighList *);
	
	// create cell offsets
	typedef void (Neighbor::*OffsetsPtr)();
	OffsetsPtr offsets_create;

	void offsets_allocate();

	void offsets_half_single();
	void offsets_half_multi();

private:
};

}

#endif

/* ERROR/WARNING messages:

E: Neighbor delay must be 0 or multiple of every setting

The delay and every parameters set via the neigh_modify command are
inconsistent.  If the delay setting is non-zero, then it must be a
multiple of the every setting.

E: Neighbor page size must be >= 10x the one atom setting

This is required to prevent wasting too much memory.

E: Invalid atom type in neighbor exclusion list

Atom types must range from 1 to Ntypes inclusive.

E: Neighbor include group not allowed with ghost neighbors

This is a current restriction within LAMMPS.

E: Neighbor multi not yet enabled for ghost neighbors

This is a current restriction within LAMMPS.

E: Neighbor multi not yet enabled for granular

Self-explanatory.

E: Neighbor multi not yet enabled for rRESPA

Self-explanatory.

E: Neighbors of ghost atoms only allowed for full neighbor lists

This is a current restriction within LAMMPS.

E: Too many local+ghost atoms for neighbor list

The number of nlocal + nghost atoms on a processor
is limited by the size of a 32-bit integer with 2 bits
removed for masking 1-2, 1-3, 1-4 neighbors.

W: Building an occasional neighobr list when atoms may have moved too far

This can cause LAMMPS to crash when the neighbor list is built.
The solution is to check for building the regular neighbor lists
more frequently.

E: Domain too large for neighbor bins

The domain has become extremely large so that neighbor bins cannot be
used.  Most likely, one or more atoms have been blown out of the
simulation box to a great distance.

E: Cannot use neighbor bins - box size << cutoff

Too many neighbor bins will be created.  This typically happens when
the simulation box is very small in some dimension, compared to the
neighbor cutoff.  Use the "nsq" style instead of "bin" style.

E: Too many neighbor bins

This is likely due to an immense simulation box that has blown up
to a large size.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid group ID in neigh_modify command

A group ID used in the neigh_modify command does not exist.

E: Neigh_modify include group != atom_modify first group

Self-explanatory.

E: Neigh_modify exclude molecule requires atom attribute molecule

Self-explanatory.

*/
