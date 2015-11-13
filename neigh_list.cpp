/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"

using namespace PDPS_NS;

#define PGDELTA 1            // page increment

/* ---------------------------------------------------------------------- */

NeighList::NeighList(PDPS *ps, int size) : Pointers(ps)
{
	pgsize = size;             // size of each page 
	
	maxpage = 0;
	
	inum = 0;
	nmax_neigh = 0;

	ilist = NULL;
	list = NULL;
	firstneigh = NULL;
	numneigh = NULL;
	list_exist = 0;
}

/* ---------------------------------------------------------------------- */

NeighList::~NeighList()
{
	memory->destroy(ilist);
	memory->destroy(numneigh);
	memory->sfree(firstneigh);
	
	for (int i = 0; i < maxpage; i++) memory->destroy(list[i]);
	memory->sfree(list);
}

/* ---------------------------------------------------------------------- 
   Grow arrary due to the increase of particles
---------------------------------------------------------------------- */

void NeighList::grow(int nmax)
{
	if (nmax_neigh >= nmax) return;
		
	nmax_neigh = nmax;

	memory->destroy(ilist);
	memory->destroy(numneigh);
	memory->sfree(firstneigh);

	memory->create(ilist,nmax_neigh,"NeighList: ilist");
	memory->create(numneigh,nmax_neigh,"NeighList: numneigh");
	firstneigh = (int **) memory->smalloc(nmax_neigh*sizeof(int *),"NeighList: firstneigh");
}

/* ---------------------------------------------------------------------- 
   Add pages
---------------------------------------------------------------------- */

int **NeighList::add_pages(int n)
{
	int npages = maxpage;
	maxpage += n;

	list = (int **) memory->srealloc(list,maxpage*sizeof(int *),"NeighList: list");

	for (int i = npages; i < maxpage; i++) {
		memory->create(list[i],pgsize,"NeighList: list[i]");
	}

	return list;
}
