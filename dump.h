/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_DUMP_H
#define PS_DUMP_H

#include "stdio.h"
#include "pointers.h"

namespace PDPS_NS {

class Dump : protected Pointers {
public:
	char *name;                  // user-defined name of Dump
	char *style;                 // style of Dump
	FILE *fp;                    // file to be dump
	int gid, groupbit;           // group that Dump is performed on
	int nevery;

	int comm_forward;                     // size of forward communication (0 if none)
	int comm_reverse;                     // size of reverse communication (0 if none)

	char boundstr[9];            // encoding of boundary flags

	Dump(class PDPS *, int, char **);
	virtual ~Dump();
	void init();
	virtual void write();
	//virtual void write();

	void modify_params(int, char **);

protected:
	int procid, nprocs;          // proc info

	double *boxlo, *boxhi, *boxle;

	char *filename;              // user-specified file
	int compressed;              // 1 if dump file is written compressed, 0 no
	int binary;                  // 1 if dump file is written binary, 0 no
	int multifile;               // 0 = one big file, 1 = one file per timestep
	int multiproc;               // 0 = proc 0 writes for all, 1 = one file/proc

	// flag
	int column_only_flag;               // 0 = no simple header print, 1 = simple header print per run
	int flush_flag;              // 0 = flush, 1 = flush every dump
	int header_flag;             // 0 = no header, 1 = write header
	int once_flag;
	
	int nme;                     // number of particles to be output in local processor

	bigint ntotal;               // total number of particles to be output 

	int maxids;
	int *ids;                    // list of particle IDs, if sorting on IDs

	int maxbuf;                
	double *buf;

	int size_one;                // # of quantities for one particle: eg. id type x y z (size_one = 5)
							     // It can be dangerous when the output requires different # of quantities 
	                             // for each field, like "dump_lammps" command 

	char *format;                // format for output

	virtual void init_style() = 0;

	virtual int count();
	virtual void write_header(bigint) = 0;
	virtual void pack(int *) = 0;
	virtual void write_data(int, double *) = 0;
};

}

#endif
