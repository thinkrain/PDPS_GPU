/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(addforce,FixAddForce)

#else

#ifndef PS_FIX_ADD_FORCE_H
#define PS_FIX_ADD_FORCE_H

#include "fix.h"

namespace PDPS_NS {

class FixAddForce : public Fix {
public:
	FixAddForce(class PDPS *, int, char **);
	~FixAddForce();
	int setmask();
	void post_force();
	//void init();
	void setup();
	//double compute_vector(int);
	//double memory_usage();

 private:
	double fx, fy, fz;
	int coupled;
	int force_style;

	// drag_stokes
	double mu;              // viscosity
	double eta;             // volume fraction

	// drag_felice
	double cle[3];
	int ncells;
	int cell[3];
	int ndims;
	double *voidage, *vol_solid, *vol_solid_total;
	int rid;

	void add_drag_stokes();
	void add_drag_general();
	void add_drag_felice();

	void compute_voidage();

	// buoyancy
	int g_dim;
	double rho, g, rho_ref;

	void add_buoyancy();

};

}

#endif
#endif
