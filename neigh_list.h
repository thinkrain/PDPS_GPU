/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_NEIGH_LIST_H
#define PS_NEIGH_LIST_H

#include "pointers.h"

namespace PDPS_NS {

class NeighList : protected Pointers {
public:
	int *ilist;            // list of local particles who need neighbor list (may not be all particles)
	int inum;

	int *numneigh;         // # of particles in the neighbor list of particle i
	int **firstneigh;      // ptr to the 1st particle in the neighbor list of particle i
	
	int pgsize;            // size of of each page
	int maxpage;           // max # of pages allocated
	int **list;            // neighbor list [page][list]

	int last_page;         // last page 
	int last_index;        // last index on last page
	
    //int *list;             // index of the particle id
	//int numlist;           // number of particles in the list
	//int nmax;              // maximum number of particles every appears in this list
	int list_exist;        // 1: exist
	                       // 0: not exist

	NeighList(class PDPS *, int);
	~NeighList();

	void grow(int);
	int **add_pages(int);


private:
	int nmax_neigh;
};

}

#endif
