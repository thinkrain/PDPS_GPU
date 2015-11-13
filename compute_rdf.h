/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(rdf,ComputeRdf)

#else

#ifndef PS_COMPUTE_RDF_H
#define PS_COMPUTE_RDF_H

#include "compute.h"

namespace PDPS_NS {

class ComputeRdf : public Compute {
public:
	ComputeRdf(class PDPS *, int, char **);
	virtual ~ComputeRdf();
	void init();
	//double compute_scalar();
	void compute_array();

protected:
	int dim;
	int nbins;                      // # of rdf bins
	int nfields;                    // # of pair combinations

	double rcut_max;
	double rbin, rbinsinv;          // bin width and its inverse

	double **gr, **grAll;           // g(r) radial distribution function
	int *ilo, *ihi, *jlo, *jhi;
	int **type2field;               // type2field[itype][jtype] = ifield
	int *typecount;
	int *icount, *jcount;           // counter for the # of itype particles and jtype particles
	                                // participating in i,j pairs
};

}

#endif
#endif
