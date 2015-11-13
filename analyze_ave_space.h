/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef ANALYZE_CLASS

AnalyzeStyle(ave/space,AnalyzeAveSpace)

#else

#ifndef PS_ANALYZE_AVE_SPACE_H
#define PS_ANALYZE_AVE_SPACE_H

#include "analyze.h"

namespace PDPS_NS {

class AnalyzeAveSpace : public Analyze {
public:
	AnalyzeAveSpace(class PDPS *, int, char **);
	virtual ~AnalyzeAveSpace();
	
	void init();

	void invoke_analyze();

private:
	//int nprocs, procid;

	int nevery, nrepeat, nfreq;
	int nstart;
	int ndims;
	int dim[3];
	int originflag[3];
	double delta[3], inv_delta[3], origin[3];
	int cid;                                     // compute id

	int ncells, maxncells;
	int cell[3];
	double **coord_cell;
	double *area_cell, *vol_cell;
	

	int icol;                                    // global icol

	// flag
	int num_flag;                                // 1: count total number of particles in each cell
	int area_flag, vol_flag;
	int *field_ave_flag;                         // type of average operation

	typedef void (AnalyzeAveSpace::*FnPtr)();
	void allocate();
	void parse_fields(int, int, char **);
	void addfield(const char *, FnPtr, int, int);
	FnPtr *field_func;                // list of ptrs to functions

	// computes
	class Compute **compute;

	void setup_cells();
	void count_num_cell();      // count total number of particles in each cell
	int find_cell_id(int);

	void write_array();
	void clean_array();

	void compute_compute();

	void compute_x();
	void compute_y();
	void compute_z();
	void compute_density_number();
	void compute_density_mass();
	void compute_density_particle();
};

}

#endif
#endif
