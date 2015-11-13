/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "domain.h"
#include "dump_lammps.h"
#include "error.h"
#include "memory.h"
#include "particle.h"
#include "update.h"

using namespace PDPS_NS;

#define MAXLINE 8192       // Maximum string length

enum{INT,DOUBLE,BIGINT};

DumpLammps::DumpLammps(PDPS *ps, int narg, char **arg) : Dump(ps, narg, arg)
{
	if (narg <= 5) error->all(FLERR,"Illegal dump custom command");

	// parse command
	keyword = NULL;
	line = NULL;

	// related to section 
	section_mpi_flag = NULL;
	section_nlines = NULL;

	// related to field
	nfields = NULL;

	// related to section and field
	vformat = NULL;
	vtype = NULL;

	// choices
	pack_choice = NULL;
	header_choice = NULL;
	write_choice = NULL;

	size_one = max_nfields = 0;

	line = new char[MAXLINE];
	line[0] = '\0';                // otherwise it will have mistakes
	for (int iarg = 5; iarg < narg; iarg++) {
	  strcat(line,arg[iarg]);
	  strcat(line," ");
	}
	line[strlen(line)-1] = '\0'; 
	
	nsections_initial = narg - 5;

	nsections = 0;

	// count number of selected particles on local processor
	nme = count();

	allocate();
	parse_field(line);	
}

/* ---------------------------------------------------------------------- */

DumpLammps::~DumpLammps()
{
	memory->destroy(vformat);
	memory->destroy(vtype);

	for (int i = 0; i < nsections; i++) {
		delete[] keyword[i];
		keyword[i] = NULL;
	}
	delete[] keyword;
	keyword = NULL;

	delete[] line;
	line = NULL;

	delete[] nfields;
	nfields = NULL;

	delete[] section_mpi_flag;
	section_mpi_flag = NULL;

	delete[] section_nlines;
	section_nlines = NULL;

	delete[] pack_choice;
	pack_choice = NULL;
}

/* ---------------------------------------------------------------------- */

void DumpLammps::allocate()
{
	int n = nsections_initial;
	
	keyword = new char*[n];

	section_nlines = new int[n];
	section_mpi_flag = new int[n];
	
	nfields = new int[n];

	pack_choice = new FnPtrPack[n];

	for (int i = 0; i < n; i++) {
		keyword[i] = NULL;
		keyword[i] = new char[32];
		nfields[i] = 0;
	}
}

/* ---------------------------------------------------------------------- */

void DumpLammps::init_style()
{
	size_one = max_nfields;       // maximum # of quantities for each particle

	header_choice = &DumpLammps::header_item;
	// pack_choice is done in add_field(const char *key, FnPtrPack func)
	write_choice = &DumpLammps::write_text;
}

/* ---------------------------------------------------------------------- */

void DumpLammps::parse_field(char *str)
{
	// customize a new keyword by adding to if statement

    char *word = strtok(str," \0");
    while (word) {
		if (strcmp(word,"atoms") == 0) {
			nfields[nsections] = 5;		
			section_nlines[nsections] = nme;
			section_mpi_flag[nsections] = 1;
			addsection("Atoms",&DumpLammps::pack_atoms);
			strcpy(vformat[nsections-1][0],"%d ");
			vtype[nsections-1][0] = INT;
			strcpy(vformat[nsections-1][1],"%d ");
			vtype[nsections-1][1] = INT;
			strcpy(vformat[nsections-1][2],"%g ");
			vtype[nsections-1][2] = DOUBLE;
			strcpy(vformat[nsections-1][3],"%g ");
			vtype[nsections-1][3] = DOUBLE;
			strcpy(vformat[nsections-1][4],"%g ");
			vtype[nsections-1][4] = DOUBLE;
		}
		else if (strcmp(word,"masses") == 0) {
			nfields[nsections] = 2;
			section_nlines[nsections] = particle->ntypes;
			section_mpi_flag[nsections] = 0;
			addsection("Masses",&DumpLammps::pack_masses);
			strcpy(vformat[nsections-1][0],"%d ");
			vtype[nsections-1][0] = INT;
			strcpy(vformat[nsections-1][1],"%g ");
			vtype[nsections-1][1] = DOUBLE;
		}
		else if (strcmp(word,"velocities") == 0) {
			nfields[nsections] = 4;
			section_nlines[nsections] = nme;
			section_mpi_flag[nsections] = 1;
			addsection("Velocities",&DumpLammps::pack_velocities);
			strcpy(vformat[nsections-1][0],"%d ");
			vtype[nsections-1][0] = INT;
			strcpy(vformat[nsections-1][1],"%g ");
			vtype[nsections-1][1] = DOUBLE;
			strcpy(vformat[nsections-1][2],"%g ");
			vtype[nsections-1][2] = DOUBLE;
			strcpy(vformat[nsections-1][3],"%g ");
			vtype[nsections-1][3] = DOUBLE;
		}
		else error->all(FLERR,"Illega dump_lammps arguments");
		word = strtok(NULL," \0");
	}
}

/* ---------------------------------------------------------------------- */

void DumpLammps::addsection(const char *key, FnPtrPack func)
{	
	if (nfields[nsections] > max_nfields) {
		max_nfields = nfields[nsections];
		memory->grow(vtype,nsections_initial,max_nfields,"DumpLammps: vtype");
		memory->grow(vformat,nsections_initial,max_nfields,32,"DumpLammps: vformat");
	}

	strcpy(keyword[nsections],key);
	pack_choice[nsections] = func;
	nsections++;
}

/* ---------------------------------------------------------------------- */

void DumpLammps::write()
{
	boxlo = domain->boxlo;
	boxhi = domain->boxhi;
	boxle = domain->boxle;

	// nme = # of selected particles (already calculated in the constructor)

	bigint bnme = nme;

	int nmax;
	int *section_nmax;                // maximum lines of each section in local processor
	section_nmax = new int[nsections];

	MPI_Allreduce(&bnme,&ntotal,1,MPI_INT,MPI_SUM,mworld); // total number of selected particles
	
	nmax = 0;
	for (int i = 0; i < nsections; i++) {
		// find the maximum lines of each field among all processors
		MPI_Allreduce(&section_nlines[i],&section_nmax[i],1,MPI_INT,MPI_MAX,mworld);
		// find the maximum lines among all fields and processors
		if (nmax < section_nmax[i])
			nmax = section_nmax[i];
	}

	if (nmax*max_nfields > maxbuf) {
		//if ((bigint) nmax * size_one > MAXSMALLINT)
			//error->all(FLERR,"Too much per-proc info for dump");
		maxbuf = nmax * max_nfields;
		memory->destroy(buf);
		memory->create(buf,maxbuf,"Dump: buf");
	}
	
	// write timestep header
	if (header_flag) {
		write_header(ntotal);
		if (once_flag == 1) header_flag = 0;
	}

	if (nmax > maxids) {
		maxids = nmax;
		memory->destroy(ids);
		memory->create(ids,maxids,"DumpLammps: ids");
	}
	
	pack(ids);

	int tmp, nlines;   // tmp ?  probably used when Lammps author debug the code
	MPI_Status status;
	MPI_Request request;
	
	for (isection = 0; isection < nsections; isection++) {
		int m = 0;
		// for field that does not require parallel output
		if (section_mpi_flag[isection] == 0) {
			if (procid == 0) {
				fprintf(fp,"%s\n",keyword[isection]);
				fprintf(fp,"\n");
				(this->*pack_choice[isection])(m);
				nlines = section_nlines[isection];
				write_data(nlines,buf);
				fprintf(fp,"\n");
				fflush(fp);
			}
			continue;
		}
		// for feild does require parallel output
		(this->*pack_choice[isection])(m);
		if (procid == 0) {
			fprintf(fp,"%s\n",keyword[isection]);
			fprintf(fp,"\n");
			for (int iproc = 0; iproc < nprocs; iproc++) {
				if (iproc > 0) {
					MPI_Irecv(buf,maxbuf,MPI_DOUBLE,iproc,0,mworld,&request);
					MPI_Send(&tmp,0,MPI_INT,iproc,0,mworld);
					MPI_Wait(&request,&status);
					MPI_Get_count(&status,MPI_DOUBLE,&nlines);
					nlines /= nfields[isection];
				} else {
					nlines = section_nlines[isection];
				}
				write_data(nlines,buf);
			}
			fprintf(fp,"\n");
			fflush(fp);
		} else {
			MPI_Recv(&tmp,0,MPI_INT,0,0,mworld,&status);
			MPI_Rsend(buf,section_nlines[isection]*nfields[isection],MPI_DOUBLE,0,0,mworld);
		} 
	} // for (int ifield = 0; ifield < nfield; ifield++)
}

/* ---------------------------------------------------------------------- */

void DumpLammps::write_header(bigint ndump)
{
	if (procid == 0) (this->*header_choice)(ndump);
}

/* ----------------------------------------------------------------------
   Choose a pack methdod
------------------------------------------------------------------------- */

void DumpLammps::pack(int *ids)
{
	int m = 0;
	int *mask = particle->mask;
	int *tag = particle->tag;
	int nlocal = particle->nlocal;

	if (ids) {
		int *tag = particle->tag;
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit)
				ids[i] = tag[i];
		}
	}
}

/* ---------------------------------------------------------------------- */

void DumpLammps::write_data(int n, double *mybuf)
{
	(this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpLammps::header_item(bigint ndump)
{
	fprintf(fp,"LAMMPS format data file (can be read by \"read_data\" command in Lammps): timestep = %d, nprocs = %d\n"
		,update->ntimestep,nprocs);
	fprintf(fp,"\n");
	fprintf(fp,BIGINT_FORMAT " atoms\n",particle->nparticles);
	fprintf(fp,"\n");
	fprintf(fp,"%d atom types\n",particle->ntypes);
	fprintf(fp,"\n");
	fprintf(fp,"%g %g xlo xhi\n",boxlo[0],boxhi[0]);
	fprintf(fp,"%g %g ylo yhi\n",boxlo[1],boxhi[1]);
	fprintf(fp,"%g %g zlo zhi\n",boxlo[2],boxhi[2]);
	fprintf(fp,"\n");
}

/* ---------------------------------------------------------------------- */

void DumpLammps::write_text(int n, double *mybuf)
{
	int m = 0;

	for (int i = 0; i < n; i++) {
		for (int ifield = 0; ifield < nfields[isection]; ifield++) {
			//fprintf(screen,"icol = %d vtype = %d vformat = %s\n",icol,vtype[ifield][icol],vformat[ifield][icol]);
			//fflush(screen);
			if (vtype[isection][ifield] == INT) fprintf(fp,vformat[isection][ifield],static_cast<int> (mybuf[m++]));
			if (vtype[isection][ifield] == DOUBLE) fprintf(fp,vformat[isection][ifield],(mybuf[m++]));
		}
		fprintf(fp,"\n");
	}
}

/* ---------------------------------------------------------------------- */

void DumpLammps::pack_atoms(int m)
{
	int *mask = particle->mask;
	int *tag = particle->tag;
	int *type = particle->type;
	double **x = particle->x;
	int nlocal = particle->nlocal;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			buf[m++] = tag[i];
			buf[m++] = type[i];
			buf[m++] = x[i][0];
			buf[m++] = x[i][1];
			buf[m++] = x[i][2];
		}
	}
}

/* ---------------------------------------------------------------------- */

void DumpLammps::pack_masses(int m)
{
	double *mass = particle->mass;
	int *type = particle->type;
	int ntypes = particle->ntypes;

	for (int i = 1; i <= ntypes; i++) {
		buf[m++] = i;
		buf[m++] = mass[i];
	}
}

/* ---------------------------------------------------------------------- */

void DumpLammps::pack_velocities(int m)
{
	int *mask = particle->mask;
	int *tag = particle->tag;
	double **v = particle->v;
	int nlocal = particle->nlocal;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			buf[m++] = tag[i];
			buf[m++] = v[i][0];
			buf[m++] = v[i][1];
			buf[m++] = v[i][2];
		}
	}
}
