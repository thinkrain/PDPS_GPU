/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"

#include "neighbor.h"
#include "neigh_list.h"
#include "output.h"
#include "parallel.h"
#include "particle.h"
#include "post_processor.h"
#include "timer.h"

#include "neighbor.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

PostProcessor::PostProcessor(PDPS *ps) : Pointers(ps)
{
	procid = parallel->procid;
	nprocs = parallel->nprocs;
}

/* ----------------------------------------------------------------------
   Output information regarding the whole simulation (multiple runs)
------------------------------------------------------------------------- */

void PostProcessor::finalize()
{
	char str[128];
	
	output->print("\n");

	sprintf(str,"Time Elapsed = %g\n", timer->time_all[TIME_LOOP]);
	output->print(str);
	
	sprintf(str,"Time spent on PAIR: %g\n",timer->time_all[TIME_PAIR]);
	output->print(str);
	
	sprintf(str,"TIME spent on NEIGHBOR: %g\n",timer->time_all[TIME_NEIGHBOR]);
	output->print(str);

	sprintf(str,"TIME spent on COMM: %g\n",timer->time_all[TIME_COMM]);
	output->print(str);

	sprintf(str,"TIME spent on OUTPUT: %g\n",timer->time_all[TIME_OUTPUT]);
	output->print(str);

	output->print("\n");

	sprintf(str,"Nbuilds = %d\n",neighbor->nbuilds);
	if (procid == 0) output->print(str);

	sprintf(str,"Ndanger = %d\n",neighbor->ndanger);
	if (procid == 0) output->print(str);

	output->print("\n");
	analyze_particles();
	output->print("\n");
	analyze_neighbors();
}

/* ----------------------------------------------------------------------
   Analyze particle information
------------------------------------------------------------------------- */

void PostProcessor::analyze_particles()
{
	int nlocal_min, nlocal_ave, nlocal_max;
	int nghost_min, nghost_ave, nghost_max;

	int nlocal = particle->nlocal;

	MPI_Allreduce(&nlocal,&nlocal_min,1,MPI_INT,MPI_MIN,mworld);
	MPI_Allreduce(&nlocal,&nlocal_max,1,MPI_INT,MPI_MAX,mworld);
	MPI_Allreduce(&nlocal,&nlocal_ave,1,MPI_INT,MPI_SUM,mworld);

	nlocal_ave = nlocal_ave / nprocs;

	char str[128];
	sprintf(str,"nlocal: %d ave %d max %d min\n",nlocal_ave,nlocal_max,nlocal_min);
	output->print(str);

	int nghost = particle->nghost;
	MPI_Allreduce(&nghost,&nghost_min,1,MPI_INT,MPI_MIN,mworld);
	MPI_Allreduce(&nghost,&nghost_max,1,MPI_INT,MPI_MAX,mworld);
	MPI_Allreduce(&nghost,&nghost_ave,1,MPI_INT,MPI_SUM,mworld);

	nghost_ave = nghost_ave / nprocs;

	sprintf(str,"nghost: %d ave %d max %d min\n",nghost_ave,nghost_max,nghost_min);
	output->print(str);
}

/* ----------------------------------------------------------------------
   Analyze neighbor information
------------------------------------------------------------------------- */

void PostProcessor::analyze_neighbors()
{
	int npages = neighbor->neighlist->maxpage;
	int last_index = neighbor->neighlist->last_index;
	int pgsize = neighbor->neighlist->pgsize;

	int num, num_min, num_max, num_ave;

	num = (npages - 1) * pgsize + last_index;
	
	MPI_Allreduce(&num,&num_min,1,MPI_INT,MPI_MIN,mworld);
	MPI_Allreduce(&num,&num_max,1,MPI_INT,MPI_MAX,mworld);
	MPI_Allreduce(&num,&num_ave,1,MPI_INT,MPI_SUM,mworld);

	num_ave = num_ave / nprocs; 
	
	char str[128];
	sprintf(str,"neighbors: %d ave %d max %d min\n",num_ave,num_max,num_min);
	output->print(str);
	
}

