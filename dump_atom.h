/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(atom,DumpAtom)

#else

#ifndef PS_DUMP_ATOM_H
#define PS_DUMP_ATOM_H

#include "dump.h"

namespace PDPS_NS {

class DumpAtom : public Dump {
public:
	DumpAtom(PDPS *, int, char**);

private:
	char *col_labels;             // column labels

	void init_style();

	void pack(int *);
	void write_header(bigint);
	void write_data(int, double *);

	typedef void (DumpAtom::*FnPtrPack)(int *);
	FnPtrPack pack_choice;                                 // ptr to pack functions
	void pack_noscale_noimage(int *);

	typedef void (DumpAtom::*FnPtrHeader)(bigint);
	FnPtrHeader header_choice;                             // ptr to write header functions
	void header_item(bigint);

	typedef void (DumpAtom::*FnPtrData)(int, double *);
	FnPtrData write_choice;
	void write_noimage(int, double *);
};

}

#endif
#endif
