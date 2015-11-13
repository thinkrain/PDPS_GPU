/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef ANALYZE_CLASS

AnalyzeStyle(order,AnalyzeOrder)

#else

#ifndef PS_ANALYZE_ORDER_H
#define PS_ANALYZE_ORDER_H

#include "analyze.h"

namespace PDPS_NS {

class AnalyzeOrder : public Analyze {
public:
	AnalyzeOrder(class PDPS *, int, char **);
	virtual ~AnalyzeOrder();

	void init();
	void invoke_analyze();

private:
	int nevery, nrepeat, nfreq;
	int nstart;					   // smallest tracerstep to do average 
	int cid;
	int rid;
	int method;
	double order;

	int ndims;
	int dim[3];
	double delta[3], inv_delta[3];

	int icol;
	int orientational_flag, translational_flag;
 
	                                                      // this one is quite memory consuming, in the future, we may 
	                                                      // use alternative way to do the mapping
	int max_ntags;                                        // max. elements allocated for tag2tracerID


	typedef void (AnalyzeOrder::*FnPtr)();
	void allocate();
	void parse_fields(int, int, char **);
	void addfield(const char *, FnPtr, int, int);
	FnPtr *field_func;             // list of ptrs to functions

	// computes
	class Compute **compute;
	void orientationalorder();
	void translationalorder();
	double legendre(int m,int l,double x);
	void sphericalharmonics(double *Y, int m,int l,double thita,double fi);
	void anglecalculate(int i, int j, double *thita, double *fi);

	void write_array();
	void clean_array();

			
};

}

#endif
#endif
