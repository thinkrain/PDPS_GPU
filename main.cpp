/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"

#include "input.h"
#include "pdps.h"

using namespace PDPS_NS;

int main(int argc, char** argv)
{
	MPI_Init(&argc,&argv);

	PDPS *pdps = new PDPS(argc,argv,MPI_COMM_WORLD);
    pdps->input->file();
	delete pdps;

	MPI_Finalize();
}
