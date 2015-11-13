/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_ANALYZE_H
#define PS_ANALYZE_H

#include "pointers.h"

namespace PDPS_NS {

class Analyze : protected Pointers {
public:
	char *name;
    char *style;
	int gid;
	int groupbit;

	int comm_forward;                     // size of forward communication (0 if none)
	int comm_reverse;                     // size of reverse communication (0 if none)

	// output flag
	int header_flag;                      // write header 
	int header_once_flag;                 // only write header once (good for plot purpose) 
	int flush_flag;
 
	// type flag
	int scalar_flag;
	int vector_flag;
	int array_flag;

	// ave flag
	int ave_flag;       // 1: do average 0: accumulate
	int clean_ave_flag;

	// invoke flag for compute
	int invoked_flag;                     // non-zero if invoked or accessed this step, 0 if not
	int invoked_scalar;                   // last timestep on which compute_scalar() was invoked
	int invoked_vector;
	int invoked_array;

	double scalar;
	double *vector;
	double **array, **array_total, **array_ave;  
	int nrows, ncolumns;

	int *num_cell, *numAll_cell;

	// field properties
	int ifield;
	int nfields_initial;                              // record number of arguments
	int nfields;
	int maxfields;
	char **field_format;
	char **field_name;
	int *field_nrows;                                 // # of rows for each field
	int *field_ncols;                                 // # of columns for each field
	int *field_index;
	int *field_data_type;
	int *field_type;

	Analyze(class PDPS *, int, char **);
	virtual ~Analyze();	

	virtual void init() {};
	virtual void setup() {};

	virtual void invoke_analyze() = 0;

	int iparticle;

	

protected:
	int nprocs, procid;
	
	// output file
	char *fname;
	FILE *file;
	

};
}

#endif
