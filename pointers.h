/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

// Pointers class contains ptrs to master copy of
//   fundamental PDPS class ptrs stored in pdps.h
// every PDPS class inherits from Pointers to access pdps.h ptrs
// these variables are auto-initialized by Pointer class constructor
// *& variables are really pointers to the pointers in pdps.h
// & enables them to be accessed directly in any class, e.g. particle->x

#ifndef PS_POINTERS_H
#define PS_POINTERS_H

//#include "pstype.h"
#include "mpi.h"
#include "pdps.h"

namespace PDPS_NS {

#define FLERR __FILE__,__LINE__

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

typedef int smallint;
typedef int bigint;
#define BIGINT_FORMAT "%d"

class Pointers {
public:
	Pointers(PDPS *ptr) :
	ps(ptr),
	memory(ptr->memory),
	error(ptr->error),
	universe(ptr->universe),
	input(ptr->input),
	particle(ptr->particle),
	update(ptr->update),
	neighbor(ptr->neighbor),
	parallel(ptr->parallel),
	domain(ptr->domain),
	force(ptr->force),
	modify(ptr->modify),
	group(ptr->group),
	output(ptr->output),
	timer(ptr->timer),
	mworld(ptr->mworld),
	infile(ptr->infile),
	screen(ptr->screen),
	//multiscale(ptr->multiscale),
	postprocessor(ptr->postprocessor),
	cudaEngine(ptr->cudaEngine),
    logfile(ptr->logfile) {}
    virtual ~Pointers() {}

protected:
	PDPS *ps;
	Memory *&memory;
	Error *&error;
	Universe *&universe;
	Input *&input;

	//CUDAParticle *&particle;
	Particle *&particle;
	Update *&update;
	Neighbor *&neighbor;
	Parallel *&parallel;
	Domain *&domain;
	Force *&force;
	Modify *&modify;
	Group *&group;
	Output *&output;
	PostProcessor *&postprocessor;
	Timer *&timer;

	CUDAEngine *&cudaEngine;

	// Multiscale test
	//MultiScale *&multiscale;

	MPI_Comm &mworld;
	FILE *&infile;
	FILE *&screen;
	FILE *&logfile;
};
}
#endif
