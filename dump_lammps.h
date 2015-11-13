/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(lammps,DumpLammps)

#else

#ifndef PS_DUMP_LAMMPS_H
#define PS_DUMP_LAMMPS_H

#include "dump.h"

namespace PDPS_NS {

class DumpLammps : public Dump {
public:
	DumpLammps(class PDPS *, int, char **);
	virtual ~DumpLammps();
	
protected:
	
	char *line;									// store commands
	char **keyword;

	// section related: 
	// section: atoms / masses / velocities
	// for each section, content in the columns will be handled by "nfields" and "max_nfields"

	int isection;
	int nsections_initial;
	int nsections;                              // number of sections
	int *section_mpi_flag;                      // buf gathered from 1 = all processors; 0 = only master processor
	int *section_nlines;                        // number of lines of each section

	// field related:
	// eg: in the atoms section field is like: id type x y z

	int max_nfields;       
	int *nfields;								// number of fields to output

	// section and field related

	char ***vformat;							// format to output
	int **vtype;							    // Int or Double
	
	void parse_field(char *);					// parse commands
	void allocate();

	void init_style();

	void write();

	void pack(int *);
	void write_header(bigint);
	void write_data(int, double *);

	typedef void (DumpLammps::*FnPtrHeader)(bigint);
	FnPtrHeader header_choice;                             // ptr to write header functions
	void header_item(bigint);

	typedef void (DumpLammps::*FnPtrData)(int,double *);
	FnPtrData write_choice;
	void write_text(int,double *);

	typedef void (DumpLammps::*FnPtrPack)(int);
	FnPtrPack *pack_choice;                                 // ptr to pack functions
	void addsection(const char *,FnPtrPack);	    // add field variable

	void pack_atoms(int);
	void pack_masses(int);
	void pack_velocities(int);    
};

}

#endif
#endif
