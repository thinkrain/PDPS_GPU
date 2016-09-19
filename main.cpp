/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */
#include "mpi.h"

#include "input.h"
#include "pdps.h"
#include "pdps_cuda.h"

using namespace PDPS_NS;

int main(int argc, char** argv)
{

	int device = 0;
	int skip = 0;
	bool skipmode = false;
	bool specified = false;
	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "-device") == 0)
		{
			i++;
			if (argv[i][0] == '-')
			{
				skipmode = true;
				skip = abs(atoi(argv[i]));
			}
			else
			{
				skipmode = false;
				device = atoi(argv[i]);
			}
			specified = true;
		}
	}

	// determine local rank of the process
	if (!specified || skipmode)
	{
		char* var;
		int dev_count, local_rank = 0;
		if ((var = getenv("SLURM_LOCALID")) != NULL) local_rank = atoi(var);
		else if ((var = getenv("MV2_COMM_WORLD_LOCAL_RANK")) != NULL) local_rank = atoi(var);
		else if ((var = getenv("OMPI_COMM_WORLD_LOCAL_RANK")) != NULL) local_rank = atoi(var);
		cudaGetDeviceCount(&dev_count);
		if (skipmode)
		{
			device = 0;
			if (device == skip) local_rank++;
			while (local_rank-- > 0)
			{
				device = (++device) % dev_count;
				if (device == skip) local_rank++;
			}
		}
		else device = local_rank % dev_count;
	}

	// override command line arguments to make sure cudaengine get the correct one
	char **argv_new = new char*[argc + 2];
	for (int i = 0; i < argc; i++)
	{
		argv_new[i] = new char[strlen(argv[i]) + 1];
		strcpy(argv_new[i], argv[i]);
	}
	argv_new[argc] = new char[32];
	argv_new[argc + 1] = new char[32];
	strcpy(argv_new[argc], "-device");
	sprintf(argv_new[argc + 1], "%d", device);
	argc += 2;
	argv = argv_new;

	cudaSetDevice(device);

	MPI_Init(&argc,&argv);

	PDPS *pdps = new PDPS(argc,argv,MPI_COMM_WORLD);
    pdps->input->file();
	delete pdps;

	MPI_Finalize();
}
