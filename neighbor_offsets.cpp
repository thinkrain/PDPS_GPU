/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "neighbor.h"
//#include "neigh_list.h"
#include "memory.h"
#include "particle.h"
#include "error.h"

using namespace PDPS_NS;

enum{SINGLE,MULTI};     // also in neighbor.cpp

/* ----------------------------------------------------------------------
   create cell offsets
------------------------------------------------------------------------- */

void Neighbor::offsets_allocate()
{
	int i, j, k;
	double f_temp;
	int i_temp;
	int co_max;

	for (i = 0; i < 3; i++) {
		f_temp = rneigh_max / cle[i];
		i_temp = static_cast <int> (f_temp);
		if (i_temp * cle[i] < rneigh_max) co[i] = i_temp + 1;
		else co[i] = i_temp;
	}

	if (dim == 2) co[2] = 0;

	cox = co[0];
 	coy = co[1];
	coz = co[2];

	// maximum offsets
	co_max = (2 *cox + 1)*(2 *coy + 1)*(2 * coz + 1);

	if (style == SINGLE) {
		// if coffsets has been allocated before
		if (co_max > max_noffsets) {
			max_noffsets = co_max;
			memory->destroy(coffsets);
			memory->create(coffsets, max_noffsets, "Neighbor: coffsets");
		}
	}
	else if (style == MULTI) {
		int n = particle->ntypes;
		if (max_noffsets_multi == 0) {
			noffsets_multi = new int[n+1];
			coffsets_multi = new int*[n+1];
			distsq_multi = new double*[n+1];
			for (i = 1; i <= n; i++) {
				noffsets_multi[i] = 0;
				coffsets_multi[i] = NULL;
				distsq_multi[i] = NULL;
			}
		}
		if (co_max > max_noffsets_multi) {
			max_noffsets_multi = co_max;
			for (i = 1; i <= n; i++) {
				memory->destroy(coffsets_multi[i]);
				memory->destroy(distsq_multi[i]);
				memory->create(coffsets_multi[i], max_noffsets_multi, "neighbor_offsets: offsets_multi[i]");
				memory->create(distsq_multi[i], max_noffsets_multi, "neighbor_offsets: distsq_multi[i]");
			}
		}
	}
}

/* ---------------------------------------------------------------------- */

void Neighbor::offsets_half_single()
{
	int i, j, k;

	offsets_allocate();
	// lammps algorithm: upper-right
	noffsets = 0;
	for (k = 0; k <= coz; k++)
	for (j = -coy; j <= coy; j++)
	for (i = -cox; i <= cox; i++) {
		if (k > 0 || j > 0 || (j == 0 && i >= 0)) {
			if (cell_distance(i,j,k) < rneighsq_max) {
				coffsets[noffsets++] = k*subncxy + j*subnc[0] + i;
			}
		}
	}
}

/* ---------------------------------------------------------------------- */

void Neighbor::offsets_half_multi()
{
	int i, j, k, n;
	double distsq, sq_type;
	int *cm;
	double *dm;

	offsets_allocate();
	
	int ntypes = particle->ntypes;
	
	for (int itype = 1; itype <= ntypes; itype++) {
		sq_type = rneighsq_type[itype];
		cm = coffsets_multi[itype];
		dm = distsq_multi[itype];
		n = 0;
		for (k = 0; k <= coz; k++)
		for (j = -coy; j <= coy; j++)
		for (i = -cox; i <= cox; i++) {
			if (k > 0 || j > 0 || (j == 0 && i >= 0)) {
				distsq = cell_distance(i,j,k);
				if (distsq < sq_type) {
					dm[n] = distsq;
					cm[n++] = k*subncxy + j*subnc[0] + i;
				}
			}
		}
		noffsets_multi[itype] = n;
	}
}
