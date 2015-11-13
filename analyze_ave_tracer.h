/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef ANALYZE_CLASS

AnalyzeStyle(ave/tracer,AnalyzeAveTracer)

#else

#ifndef PS_ANALYZE_AVE_TRACER_H
#define PS_ANALYZE_AVE_TRACER_H

#include "analyze.h"

namespace PDPS_NS {

class AnalyzeAveTracer : public Analyze {
public:
	AnalyzeAveTracer(class PDPS *, int, char **);
	virtual ~AnalyzeAveTracer();

	void init();
	void invoke_analyze();

private:
	int nevery, nrepeat, nfreq;
	int nstart;					   // smallest tracerstep to do average 
	int cid;
	int rid;

	int ndims;
	int dim[3];
	double delta[3], inv_delta[3];

	int icol;

	// variables related with tracer
    int nx, ny, nz;			                              // division of the box used for the tracer 
	double txle, tyle, tzle;	                          // length of each tracer
    int ntracers, itracer;                                // total number of the tracer and the present tracer
	int *tag2tracerID_local, *tag2tracerID;               // tag2tracerID maps the global tag ID of particle into the tracer
	                                                      // this one is quite memory consuming, in the future, we may 
	                                                      // use alternative way to do the mapping
	int max_ntags;                                        // max. elements allocated for tag2tracerID

	int tracer[3];

	typedef void (AnalyzeAveTracer::*FnPtr)();
	void allocate();
	void parse_fields(int, int, char **);
	void addfield(const char *, FnPtr, int, int);
	FnPtr *field_func;             // list of ptrs to functions

	// computes
	class Compute **compute;

	void write_array();
	void clean_array();
			
	void compute_centroid();	   // compute the centroid of each tracer
    void setup_tracers();		   // setup tracers
	void assign_tracerID();        // assign tracer ID to all particles
};

}

#endif
#endif
