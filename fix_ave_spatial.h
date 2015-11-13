/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(ave/spatial,FixAveSpatial)

#else

#ifndef PS_FIX_AVE_SPATIAL_H
#define PS_FIX_AVE_SPATIAL_H

#include "fix.h"

namespace PDPS_NS {

class FixAveSpatial : public Fix {

public:

	FixAveSpatial(class PDPS *, int, char **);
	~FixAveSpatial();
	int setmask();
	
	void init();

	//void setup();
	virtual void end_of_step();
	
protected:
	int nevery;             // number of steps to collect sample points
	int nrepeat;            // number of sample points you want to collect
	int nfreq;              // frequency to compute average
	int ago;

	int g_start;
	int g_end;
	int n_start;
	int n_end;
	
	char *fname;
	FILE *file;

	int rid; 
	double xlo, xhi, xle, ylo, yhi, yle, zlo, zhi, zle;
	char *rname;            // region's name

	int dim[3];             // dimension flag
	int nbins[3];           // number of bins
	double *hx;              // coordinate of each bin
	double *hy;
	double *hz;
	double dh[3];              // height of each bin

	char **format;		    // format to output
	char *line;				// store commands
	char **keyword;         // store keywords
	char *header;           // store variables to be output in the header

	int iter;               // number of timesteps used in one average
	int counter;            // counter of the number of averages done
	int average_flag;

	int iparticle;
	int ibin;
	int ifield;
	int nfield_initial;
	int nfield;									// number of fields to output
	int *vtype;									// Int or Double
	int *field2index;

	int nbivalue;                               // number of bivalue to store
	int bivalue;
	int **bivalues_bin;                         // field of bivalues based on bins
	int **bivalues_t;                           // field of bivalues based on time
	int ndvalue;                                // number of dvalue to store
	double dvalue;
	double **dvalues_bin;                       // field of dvalues based on bins
	double **dvalues_t;                         // field of dvalues based on time
	int *num_bins;

	void parse_field(char *);					// parse commands

	typedef void (FixAveSpatial::*FnPtr)();
	void addfield(const char *, FnPtr, int);	// add field variable
	FnPtr *vfunc;    

	void allocate();
	void average();
	int check_bins(int);
	void clear();
	void ensemble_average();
	void time_average();

	void write_header();
	void write();

	void write_vx();
	void write_vy();
	void write_vz();
};
}

#endif
#endif
