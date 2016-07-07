/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(custom,DumpCustom)

#else

#ifndef PS_DUMP_CUSTOM_H
#define PS_DUMP_CUSTOM_H

#include "dump.h"

namespace PDPS_NS {

class DumpCustom : public Dump {
public:
	DumpCustom(class PDPS *, int, char **);
	virtual ~DumpCustom();
	
protected:
	char *col_labels;
	
	int nselected; 
	int max_selected;
	int *selected_list;                         // list of selected particle local id

	char **vformat;								// format to output
	char *line;									// store commands
	char **keyword;
	
	int ifield;
	int nfields_initial;
	int nfields;									// number of fields to output
	int *vtype;									// Int or Double

	int iparticle;
	
	int bivalue;
	double dvalue;

	void parse_field(char *);					// parse commands

	void init_style();

	int count();
	void pack(int *);
	void write_header(bigint);
	void write_data(int, double *);

	typedef void (DumpCustom::*FnPtrHeader)(bigint);
	FnPtrHeader header_choice;                             // ptr to write header functions
	void header_item(bigint);
	void header_column_only(bigint);

	typedef void (DumpCustom::*FnPtrData)(int, double *);
	FnPtrData write_choice;
	void write_text(int, double *);

	typedef void (DumpCustom::*FnPtrPack)(int);
	FnPtrPack *pack_choice;                                 // ptr to pack functions
	void addfield(const char *,FnPtrPack,int);	// add field variable

	void pack_id(int);
	void pack_procid(int);
	void pack_fx(int);
	void pack_fy(int);
	void pack_fz(int);
	void pack_step(int);
	void pack_type(int);
	void pack_vx(int);
	void pack_vy(int);
	void pack_vz(int);
	void pack_x(int);
	void pack_y(int);
	void pack_z(int);
	void pack_radius(int);
	void pack_density(int);
	void pack_energy(int);
	void pack_rho(int);

	void allocate();    
	void write_fx();
	void write_fy();
	void write_fz();
	void write_id();
	void write_step();
	void write_type();
	void write_vx();
	void write_vy();
	void write_vz();
	void write_x();
	void write_y();
	void write_z();
};

}

#endif
#endif
