/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "neighbor.h"
#include "neigh_list.h"
#include "particle.h"
#include "error.h"

using namespace PDPS_NS;

/* ----------------------------------------------------------------------
   binned neighbor list construction with full Newton's 3rd law
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::half_linkedlist(NeighList *neighlist)
{
	int *tag = particle->tag;
	int *type = particle->type;
	int itype, jtype;

	int *ilist = neighlist->ilist;
	int *numneigh = neighlist->numneigh;
	int **firstneigh = neighlist->firstneigh;
	int **list = neighlist->list;

	int i, j, k;
	int inum = 0;

	int c, c1;                 // cell index and neighbor cell index
	double **x = particle->x;


	// Scan local particles

	int current_page = 0;      // current index of page
	int current_index = 0;     // current index of list on current page

	int nlocal = particle->nlocal;

	for (i = 0; i < nlocal; i++) {
		// check if we need to add pages
		if (pgsize - current_index < onesize) {
			current_page++;
			if (current_page == neighlist->maxpage) {
				list = neighlist->add_pages(1);
			}
			current_index = 0;
		}

		// Calculate a scalar cell index based on the local domain
		c = coord2cell(x[i]);

		// make sure local cells are setup correctly
		if (c < 0 || c >= subncxyz ) { 
			error->all(FLERR,"Invalid cell setup");
		}
				
		ilist[inum] = i;
		itype = type[i];
		// reset # of neighbor particles
		numneigh[inum] = 0;
		// Scan neighbor cells including itself
		
		for (int io = 0; io < noffsets; io++) {
			// find the cell index to be scan
			c1 = c + coffsets[io];

			// if c1 == c, the starting particle is linked_list[i] 
			// because the local id is monotonically increasing in the linked_list[i]
			if (c == c1) {
				j = linked_list[i];
			}
		 	else j = head[c1];
			
			// make sure local cells are setup correctly
			if (c1 < 0 || c1 >= subncxyz ) { 
				error->all(FLERR,"Invalid cell setup");
			}

			// Scan particle j in cell c1
			while (j != -1) {
				jtype = type[j];
				// scan rest of particles in this cell c1
				if (c == c1) {
					// if j is owned particle, store it, since j > i in the linked list
					// if j is ghost particle, only store when tag[i] < tag[j]
					// in LAMMPS, the author uses the rule "j coords are above and to the right of i"
					if (j >= nlocal) {
						if (tag[i] >= tag[j]) {
							j = linked_list[j];
							continue;
						}
					}
				}					
				
				rijsq = 0.0;
				rijsq = find_distance(i,j);
				// Build neighbor list for particle i
				if (rijsq < rneighsq[itype][jtype]) {
					// record the first neighbor particle's index's pointer address
					if (numneigh[inum] == 0) firstneigh[inum] = &list[current_page][current_index];
					
					list[current_page][current_index++] = j;
					numneigh[inum]++;

				}
				j = linked_list[j];
			} // while (j != -1)			
		} // for (int io = 0; io < noffsets; io++)
		
		// if neighbor list is over flow, may need to adjust onesize
		if (numneigh[inum] > onesize) {
			error->all(FLERR,"neighbor list for one particle is too large. Modify onesize by neigh_modify");
		}
		inum++;
		/* speed test purpose
		if (i % 1000 == 0) {
			fprintf(stdout, "particle = %d\n", i);
		}*/
	} // for (int i = 0; i < nlocal; i++)

	neighlist->inum = inum;
	neighlist->last_page = current_page;
	neighlist->last_index = current_index;
}
