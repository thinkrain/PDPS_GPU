/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"

#include "error.h"
#include "output.h"
#include "parallel.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

Error::Error(PDPS *ps) : Pointers(ps) {}

/* ----------------------------------------------------------------------
   called by all procs in the mpi world
   close all output, screen, and log files in world
   insure all procs in world call, else will hang
   force MPI_Abort if running in multi-partition mode
------------------------------------------------------------------------- */

void Error::all(const char *file, int line, const char *str)
{
	MPI_Barrier(mworld);

	int procid = parallel->procid;
	
	if (procid == 0) {
		if (screen) fprintf(screen,"Error: %s (%s: line %d)\n",str,file,line);
		if (logfile) fprintf(logfile,"Error: %s (%s: line %d)\n",str,file,line);
	}
	
	if (output) delete output;
	if (screen && screen != stdout) fclose(screen);
	if (logfile) fclose(logfile);

	MPI_Finalize();
	exit(1);
}

/* ----------------------------------------------------------------------
   called by one proc in the mpi world
   write to world screen only if non-NULL on this proc
   always write to universe screen
   forces abort of entire world (and universe) if any proc in world calls
------------------------------------------------------------------------- */

void Error::one(const char *file, int line, const char *str)
{
	int procid;

	procid = parallel->procid;

	if (screen) fprintf(screen,"Error on Processor %d: %s (%s: line %d)\n",
					  procid,str,file,line);

	MPI_Abort(mworld,1);
	exit(1);
}

/* ----------------------------------------------------------------------
   Warning message
------------------------------------------------------------------------- */

void Error::warning(const char *file, int line, const char *str)
{
	int procid = parallel->procid;
	
	if (procid == 0) {
		if (screen) fprintf(screen,"Warning: %s (%s: line %d)\n",str,file,line);
		if (logfile) fprintf(logfile,"Warning: %s (%s: line %d)\n",str,file,line);
	}
}

/* ----------------------------------------------------------------------
   shutdown LAMMPS
   called by all procs in one world
   close all output, screen, and log files in world
   no abort, so insure all procs in world call, else will hang
------------------------------------------------------------------------- */

void Error::done()
{
  MPI_Barrier(mworld);

  if (output) delete output;
  if (screen && screen != stdout) fclose(screen);
  if (logfile) fclose(logfile);

  MPI_Finalize();
  exit(1);
}
