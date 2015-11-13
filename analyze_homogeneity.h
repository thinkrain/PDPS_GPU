/* ----------------------------------------------------------------------
PDPS - Particle Dynamics Parallel Simulator

Copyright (2012) reserved by Lingqi Yang.
Email: ly2282@columbia.edu

See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef ANALYZE_CLASS

AnalyzeStyle(homogeneity, AnalyzeHomogeneity)

#else

#ifndef PS_ANALYZE_HOMOGENEITY_H
#define PS_ANALYZE_HOMOGENEITY_H

#include "analyze.h"

namespace PDPS_NS {

	class AnalyzeHomogeneity : public Analyze {
	public:
		AnalyzeHomogeneity(class PDPS *, int, char **);
		virtual ~AnalyzeHomogeneity();
		void init();
		//double compute_scalar();
		//void compute_array();

		void invoke_analyze();

	private:
		int nevery, nrepeat, nfreq;
		int nstart;         // smallest timestep to do average 
		int aid;            // analyze id
		int cid;            // compute id

		int icol; 

		int num_min;        // if the # of particles is less than num_min, 
		                    // the value of this cell will not be counted

		// homogeneity method
		int method;

		// local average
		int eta_s_flag;
		double *eta_s, *eta_r;        // eta_s: fully segregated; eta_r: fully mixed 
		double *eta;                  // mixing index

		// centroid

		// flag
		int scalar_flag;
		int vector_flag;
		int array_flag;

		int ave_flag;               // 1: do average 0: accumulate
		int clean_ave_flag;

		typedef void (AnalyzeHomogeneity::*FnPtr)();
		void allocate();
		void parse_fields(int, int, char **);
		void addfield(const char *, FnPtr, int, int);
		FnPtr *field_func;                // list of ptrs to functions

		// computes
		class Compute **compute;
		class Analyze **analyze;

		void write_array();
		void clean_array();

		// compute different kinds of properties

		void compute_compute();
		void compute_analyze();

		void compute_local_average(double **, int *, int, int);
		void compute_centroid(double **, int *, int, int);
	};

}

#endif
#endif
