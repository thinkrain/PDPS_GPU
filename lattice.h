/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_LATTICE_H
#define PS_LATTICE_H

#include "pointers.h"

namespace PDPS_NS {

class Lattice : protected Pointers {
public:
	int style;
	
	int ncxy, ncxyz;
	int nc[3];
	double cle[3];
	int *cell_num;

	//double **cell_coord;

	Lattice(class PDPS *, int, char **);
	~Lattice();
	
	int coord2cell(double *, int &, int &, int &);

private:
	class Region *region;
	int rid;
	double xlo, xhi, xle, ylo, yhi, yle, zlo, zhi, zle;
	double *extent_lo, *extent_hi, *extent_le;
	
	int subnc[3];
	

	int seed;

	double rcell;
	void setup_cell();
	
};

}

#endif

