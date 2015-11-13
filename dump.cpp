/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include "domain.h"
#include "dump.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "output.h"
#include "parallel.h"
#include "particle.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

Dump::Dump(PDPS *ps, int narg, char **arg) : Pointers(ps)
{
	procid = parallel->procid;
	nprocs = parallel->nprocs;

	flush_flag = 0;
	header_flag = 1;
	once_flag = 0;
	column_only_flag = 0;	

	comm_forward = comm_reverse = 0;

	column_only_flag = 0;

	filename = NULL;
	fp = NULL;
	name = NULL;
	style = NULL;

	int n = strlen(arg[0]) + 1;
	name = NULL;
	name = new char[n];
	strcpy(name,arg[0]);

	gid = group->find_group(arg[1]);
	groupbit = group->bitmask[gid];
	if (gid == -1) {
		char str[128];
		sprintf(str,"Cannot find group %s\n",arg[1]);
		error->all(FLERR,str);
	}

	n = strlen(arg[2]) + 1;
	style = new char[n];
	strcpy(style,arg[2]);

	nevery = atoi(arg[3]);
	if (nevery < 0) error->all(FLERR,"Illegal dump frequnecy");

	n = strlen(arg[4]) + 1;
	filename = new char[n];
	strcpy(filename,arg[4]);

	if (procid == 0) {
		fp = fopen(filename,"w");
		if (fp == NULL) {
			char str[128]; 
			sprintf(str,"Cannot open dump file %s",filename);
			error->all(FLERR,str);
		}
	}

	maxbuf = 0;
	buf = NULL;

	format = NULL;
}

/* ---------------------------------------------------------------------- */

Dump::~Dump()
{
	delete[] name;
	name = NULL;
	
	delete[] style;
	style = NULL;

	delete[] filename;
	filename = NULL;

	if (fp && procid == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

void Dump::init()
{
	init_style();

	maxids = 0;
	ids = NULL;

}

/* ---------------------------------------------------------------------- */

void Dump::write()
{
	boxlo = domain->boxlo;
	boxhi = domain->boxhi;
	boxle = domain->boxle;

	// nme = # of dump lines this proc will contribute to dump

	nme = count();
	bigint bnme = nme;

	int nmax;

	MPI_Allreduce(&bnme,&ntotal,1,MPI_INT,MPI_SUM,mworld);
    MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,mworld);

	// write timestep header
	if (header_flag) {
		write_header(ntotal);
		if (once_flag == 1) header_flag = 0;
	}
	
	if (nmax > maxbuf) {
		//if ((bigint) nmax * size_one > MAXSMALLINT)
			//error->all(FLERR,"Too much per-proc info for dump");
		maxbuf = nmax;
		memory->destroy(buf);
		memory->create(buf,maxbuf*size_one,"Dump: buf");
	}

	if (nmax > maxids) {
		maxids = nmax;
		memory->destroy(ids);
		memory->create(ids,maxids,"Dump: ids");
	}
	
	pack(ids);

	int tmp, nlines;   // tmp ?  probably used when Lammps author debug the code
	MPI_Status status;
	MPI_Request request;

	if (procid == 0) {
		for (int iproc = 0; iproc < nprocs; iproc++) {
			if (iproc > 0) {
				MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,iproc,0,mworld,&request);
				MPI_Send(&tmp,0,MPI_INT,iproc,0,mworld);
				MPI_Wait(&request,&status);
				MPI_Get_count(&status,MPI_DOUBLE,&nlines);
				nlines /= size_one;
			} else {
				nlines = nme;
			}

			write_data(nlines,buf);
		}
		fflush(fp);
	} else {
		MPI_Recv(&tmp,0,MPI_INT,0,0,mworld,&status);
		MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,0,0,mworld);
	}
}

/* ---------------------------------------------------------------------- */

int Dump::count()
{
	if (gid == 0) return particle->nlocal;

	int *mask = particle->mask;
	int nlocal = particle->nlocal;

	int m = 0;
	for (int i = 0; i < nlocal; i++)
		if (mask[i] & groupbit) m++;
	return m;
}

/* ----------------------------------------------------------------------
   Modify dump parameters
------------------------------------------------------------------------- */

void Dump::modify_params(int narg, char **arg)
{
	int iarg = 0;
	
	while (iarg < narg) {
		if (strcmp(arg[iarg],"column_only") == 0) {
			if (strcmp(arg[iarg+1],"yes") == 0) column_only_flag = 1;
			else if (strcmp(arg[iarg+1],"no") == 0) column_only_flag = 0;
			else error->all(FLERR,"Illegal dump_modify only_column command");
			iarg += 2;
		}
		else if (strcmp(arg[iarg],"flush") == 0) {
			if (strcmp(arg[iarg+1],"yes") == 0) flush_flag = 1;
			else if (strcmp(arg[iarg+1],"no") == 0) flush_flag = 0;
			else error->all(FLERR,"Illegal dump_modify flush command");
			iarg += 2;
		}
		else if (strcmp(arg[iarg],"header") == 0) {
			if (strcmp(arg[iarg+1],"yes") == 0) header_flag = 1;
			else if (strcmp(arg[iarg+1],"no") == 0) header_flag = 0;
			else error->all(FLERR,"Illegal dump_modify header command");
			iarg += 2;
		}
		else if (strcmp(arg[iarg],"once") == 0) {
			if (strcmp(arg[iarg+1],"yes") == 0) once_flag = 1;
			else if (strcmp(arg[iarg+1],"no") == 0) once_flag = 0;
			else error->all(FLERR,"Illegal dump_modify once command");
			iarg += 2;
		}
		else error->all(FLERR,"Illegal dump_modify command");
	}
}
