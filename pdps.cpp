/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Main Program to Drive PDPS
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"

#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
//#include "multiscale.h"
#include "neighbor.h"
#include "output.h"
#include "parallel.h"
#include "particle.h"
//#include "cuda_particle.h"
#include "pdps.h"
#include "post_processor.h"
#include "timer.h"
#include "update.h"
#include "cuda_engine.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

PDPS::PDPS(int narg, char** arg, MPI_Comm world) 
{
	mworld = world;

	int    device = 0;
	string profile;

	int iarg = 1;
	//	set the GPU information
	while (iarg < narg){
		if (strcmp(arg[iarg], "-device") == 0) {
			if (iarg + 1 > narg) error->all(FLERR, "Invalid command-line argument");
			device = atoi(arg[++iarg]);
			iarg++;
		}
		else if (strcmp(arg[iarg], "-profile") == 0) {
			if (iarg + 1 > narg) error->all(FLERR, "Invalid command-line argument");
			profile = arg[++iarg];
			iarg++;
		}
	}

	parallel = NULL;
	parallel = new Parallel(this);

	memory = NULL;
	memory = new Memory(this);
	error = NULL;
	error = new Error(this);
	output = NULL;
    output = new Output(this);

	cudaEngine = NULL;
	cudaEngine = new CUDAEngine(this, device, profile);
//	class CUDAEngine *cudaEngine;
	int inflag = 0;
	int screenflag = 0;
	int logflag = 0;

	infile = NULL;
	logfile = NULL;
	screen = stdout;
	
	int procid = parallel->procid;
	// only master processor open file
	if (procid == 0) {
		// infile
		//arg[1] = "liquid-water2D.in";                      // for debug only
		char *temp;
		temp = new char[20];
		strcpy(temp,"liquid.in");
		infile = fopen(temp,"r");                // open input file
		
		// logfile
		logfile = fopen("log_pdps.txt","w");

		// screen
		//char *name_screen = "screen.txt";          // for debug only
		//screen = fopen(name_screen,"w");
		char str[128];
		if (infile == NULL) {
			sprintf(str,"Cannot open input script \"%s\"\n",arg[1]);
			error->all(FLERR,str);
		}
		else {
			sprintf(str,"Open input script \"%s\" successfully\n",arg[1]);
			output->print(str);
		}
	}

	//fprintf(screen,"%d\n",sizeof(smallint));
	//fprintf(screen,"%d\n",sizeof(bigint));

	input = NULL;
	input = new Input(this,narg,arg);          // initialize Input class
	
	create();
}

/* ---------------------------------------------------------------------- */

PDPS::~PDPS()
{
	destroy();

	if (logfile) fclose(logfile);
	if (screen && screen != stdout) fclose(screen);

	delete input; 
	delete error;
	delete memory;
}

/* ---------------------------------------------------------------------- */

void PDPS::create()
{
	// Initialization of top-level class
	particle = NULL;
	particle = new Particle(this);

	domain = NULL;
	domain = new Domain(this);
	
	modify = NULL;
	modify = new Modify(this);

	group = NULL;
	group = new Group(this);

	neighbor = NULL;
    neighbor = new Neighbor(this);

	timer = NULL;

	force = NULL;
	force = new Force(this);

	update = NULL;
	update = new Update(this);

	timer = NULL;
	timer = new Timer(this);

	// Multiscale test
	//multiscale = NULL;
	//multiscale = new MultiScale(this);

	postprocessor = NULL;
	postprocessor = new PostProcessor(this);
}

/* ----------------------------------------------------------------------
   initialize top-level classes
   do not initialize Timer class, other classes like Run() do that explicitly
------------------------------------------------------------------------- */

void PDPS::init()
{
	update->init();
	force->init();
	domain->init();                 
	particle->init();               // particle must come after force and domain
    modify->init();                 // initiate modify 
	neighbor->init();
	parallel->init();
	output->init();
	timer->init();
}

/* ----------------------------------------------------------------------
   delete single instance of top-level classes
   fundamental classes are deleted in destructor
------------------------------------------------------------------------- */

void PDPS::destroy()
{
	delete update;
	delete neighbor;
	delete force;
	delete group;
	delete output;
	delete modify;               // modify must come after output, force, update
	delete domain;               // domain must come after modify
	delete particle;             // particle must come after modify, neighbor
	delete postprocessor;
	delete timer;
	delete cudaEngine;
}

