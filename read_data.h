/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(read_data,ReadData)

#else

#ifndef PS_READ_DATA_H
#define PS_READ_DATA_H

#include "stdio.h"
#include "pointers.h"

namespace PDPS_NS {

class ReadData : protected Pointers {
public:
	ReadData(class PDPS *);
	~ReadData();
	void command(int, char **);

private:
	int procid, nprocs;
	char *line,*keyword,*buffer;
	FILE *fp;
	int narg, maxarg;
	char **arg;

	void open(char *);
	void header(int);
	void parse_keyword(int, int);

	void atoms();
	void velocities();
	void masses();
};

}

#endif
#endif
