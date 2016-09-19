/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator/*
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(couplesph,FixCouplesph)

#else

#ifndef PS_FIX_COUPLESPH_H
#define PS_FIX_COUPLESPH_H

#include "fix.h"

namespace PDPS_NS {

class FixCouplesph : public Fix {

public:
	FixCouplesph(class PDPS *, int, char **);
	~FixCouplesph();
	int setmask();
	void init();
	void setup();
	virtual void post_force();
	
protected:
	int poro_flag, poroset_flag;
	double poroset;
	double a2D, a3D, h;
	int cubic_flag, quintic_flag;
	double rho_ref;
	int phase_f, phase_s;
	int sgid, sgroupbit, lgid, lgroupbit;					//	group bit for liquid particles and solid particles

};
}

#endif
#endif
