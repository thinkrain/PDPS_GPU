/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

// system library
#include "mpi.h"
#include "stdlib.h"
#include "string.h"

// pdps library
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "parallel.h"
#include "particle.h"
#include "particle_type.h"
#include "read_data.h"
#include "update.h"

using namespace PDPS_NS;

#define MAXLINE 256
#define LB_FACTOR 1.1
#define CHUNK 1024
#define DELTA 4

#define NSECTIONS 23       // change when add to header::section_keywords

/* ---------------------------------------------------------------------- */

ReadData::ReadData(PDPS *ps) : Pointers(ps)
{
	procid = parallel->procid;
	nprocs = parallel->nprocs;

	line = new char[MAXLINE];
	keyword = new char[MAXLINE];
	buffer = new char[CHUNK*MAXLINE];
	narg = maxarg = 0;
	arg = NULL;
}

/* ---------------------------------------------------------------------- */

ReadData::~ReadData()
{
	delete[] line;
	line = NULL;
	
	delete[] keyword;
	keyword = NULL;
	
	delete[] buffer;
	buffer = NULL;
}

/* ---------------------------------------------------------------------- */

void ReadData::command(int narg, char **arg)
{
	if (narg < 1) error->all(FLERR,"Illegal read_data command");

	if (domain->box_exist)
		error->all(FLERR,"Cannot use read_data command after simualtion box is defined");
	if (domain->dim == 2 && domain->zperiodic == 0)
		error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");

	// read header info

	if (procid == 0) {
		if (screen) fprintf(screen,"Reading data file ... \n");
		open(arg[0]);
	}
	header(1);
	domain->box_exist = 1;

	// problem setup using info from header

	update->ntimestep = 0;

	int n;
	if (nprocs == 1) n = static_cast<int> (particle->nparticles);
	else n = static_cast<int> (LB_FACTOR * particle->nparticles / nprocs);

	particle->allocate_type_arrays();
	particle->ptype->grow(n);
	n = particle->nmax;

	domain->print_box("created");
	domain->set_initial_box();
	domain->set_global_box();
	parallel->set_proc_grid();
	domain->set_local_box();

	// customize for new sections
	// read rest of file in free format

	int atomflag = 0;

	while (strlen(keyword)) {
		if (strcmp(keyword,"Atoms") == 0) {
			atoms();
			atomflag = 1;
		} 
		else if (strcmp(keyword,"Masses") == 0) {
			masses();
		}
		else if (strcmp(keyword,"Velocities") == 0) {
			if (atomflag == 0) error->all(FLERR,"Must read Atoms before Velocities");
			velocities();
		}
		else {
			char str[128];
			sprintf(str,"Unknown identifier in data file: %s",keyword);
			error->all(FLERR,str);
		}
		parse_keyword(0,1);
	}

	// close file

	if (procid == 0) {
		fclose(fp);
	}

	// error if natoms > 0 yet no atoms were read

	if (particle->nparticles > 0 && atomflag == 0) {
		error->all(FLERR,"No atoms in data file");
	}
}

/* ----------------------------------------------------------------------
   read free-format header of data file
   if flag = 0, only called by proc 0
   if flag = 1, called by all procs so bcast lines as read them
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
   header ends with EOF or non-blank line containing no header keyword
     if EOF, line is set to blank line
     else line has first keyword line for rest of file
------------------------------------------------------------------------- */

void ReadData::header(int flag)
{
	int n;
	char *ptr;

	// customize for new sections

	const char *section_keywords[NSECTIONS] =
	{"Atoms","Velocities","Ellipsoids","Lines","Triangles",
	 "Bonds","Angles","Dihedrals","Impropers",
	 "Masses","Pair Coeffs","Bond Coeffs","Angle Coeffs",
	 "Dihedral Coeffs","Improper Coeffs",
	 "BondBond Coeffs","BondAngle Coeffs","MiddleBondTorsion Coeffs",
	 "EndBondTorsion Coeffs","AngleTorsion Coeffs",
	 "AngleAngleTorsion Coeffs","BondBond13 Coeffs","AngleAngle Coeffs"};

	// skip 1st line of file

	if (procid == 0) {
		char *eof = fgets(line,MAXLINE,fp);
		if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
	}

	// customize for new header lines

	while (1) {

		// read a line and bcast length if flag is set

		if (procid == 0) {
			if (fgets(line,MAXLINE,fp) == NULL) n = 0;
			else n = strlen(line) + 1;
		}
		if (flag) MPI_Bcast(&n,1,MPI_INT,0,mworld);

		// if n = 0 then end-of-file so return with blank line

		if (n == 0) {
			line[0] = '\0';
			return;
		}

		// bcast line if flag is set

		if (flag) MPI_Bcast(line,n,MPI_CHAR,0,mworld);

		// trim anything from '#' onward
		// if line is blank, continue

		if (ptr = strchr(line,'#')) *ptr = '\0';
		if (strspn(line," \t\n\r") == strlen(line)) continue;

		// search line for header keyword and set corresponding variable

		if (strstr(line,"atoms")) sscanf(line,BIGINT_FORMAT,&particle->nparticles);

		// check for these first
		// otherwise "triangles" will be matched as "angles"

		else if (strstr(line,"atom types")) sscanf(line,"%d",&particle->ntypes);
		else if (strstr(line,"xlo xhi"))
			sscanf(line,"%lg %lg",&domain->boxlo[0],&domain->boxhi[0]);
		else if (strstr(line,"ylo yhi"))
			sscanf(line,"%lg %lg",&domain->boxlo[1],&domain->boxhi[1]);
		else if (strstr(line,"zlo zhi"))
			sscanf(line,"%lg %lg",&domain->boxlo[2],&domain->boxhi[2]);
		else break;
	}

	// error check on total system size

	if (particle->nparticles < 0 || particle->nparticles > 2E31) {
		if (procid == 0) error->one(FLERR,"System in data file is too big");
	}

	// check that exiting string is a valid section keyword

	parse_keyword(1,flag);
	for (n = 0; n < NSECTIONS; n++)
		if (strcmp(keyword,section_keywords[n]) == 0) break;
	if (n == NSECTIONS && procid == 0) {
		char str[128];
		sprintf(str,"Unknown identifier in data file: %s",keyword);
		error->one(FLERR,str);
	}
}

/* ----------------------------------------------------------------------
   grab next keyword
   read lines until one is non-blank
   keyword is all text on line w/out leading & trailing white space
   read one additional line (assumed blank)
   if any read hits EOF, set keyword to empty
   if first = 1, line variable holds non-blank line that ended header
   if flag = 0, only proc 0 is calling so no bcast
   else flag = 1, bcast keyword line to all procs
------------------------------------------------------------------------- */

void ReadData::parse_keyword(int first, int flag)
{
	int eof = 0;

	// proc 0 reads upto non-blank line plus 1 following line
	// eof is set to 1 if any read hits end-of-file

	if (procid == 0) {
		if (!first) {
			if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
		}
		while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
			if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
		}
		if (fgets(buffer,MAXLINE,fp) == NULL) eof = 1;
	}

  // if eof, set keyword empty and return

	if (flag) MPI_Bcast(&eof,1,MPI_INT,0,mworld);
	if (eof) {
		keyword[0] = '\0';
		return;
	}

  // bcast keyword line to all procs

	if (flag) {
		int n;
		if (procid == 0) n = strlen(line) + 1;
		MPI_Bcast(&n,1,MPI_INT,0,mworld);
		MPI_Bcast(line,n,MPI_CHAR,0,mworld);
	}

  // copy non-whitespace portion of line into keyword

	int start = strspn(line," \t\n\r");
	int stop = strlen(line) - 1;
	while (line[stop] == ' ' || line[stop] == '\t'
		|| line[stop] == '\n' || line[stop] == '\r') stop--;
	line[stop+1] = '\0';
	strcpy(keyword,&line[start]);
}

/* ----------------------------------------------------------------------
   procid 0 opens data file
   Include read zip file in the future
------------------------------------------------------------------------- */

void ReadData::open(char *file)
{
	char *suffix = file + strlen(file) - 3;

	fp = fopen(file,"r");

	if (fp == NULL) {
		char str[128];
		sprintf(str,"Cannot open file %s",file);
		error->one(FLERR,str);
	}
}

/* ----------------------------------------------------------------------
   Read all atoms
------------------------------------------------------------------------- */

void ReadData::atoms()
{
	int i,m,nchunk;

	bigint nread = 0;
	bigint nparticles = particle->nparticles;

	while (nread < nparticles) {
		if (nparticles-nread > CHUNK) nchunk = CHUNK;
		else nchunk = nparticles-nread;
		if (procid == 0) {
			char *eof;
			m = 0;
			for (i = 0; i < nchunk; i++) {
				eof = fgets(&buffer[m],MAXLINE,fp);
				if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
				m += strlen(&buffer[m]);
			}
			m++;
		}
		MPI_Bcast(&m,1,MPI_INT,0,mworld);
		MPI_Bcast(buffer,m,MPI_CHAR,0,mworld);

		particle->data_particles(nchunk,buffer);
		nread += nchunk;
	}

  // check that all atoms were assigned correctly

	bigint tmp = particle->nlocal;
	MPI_Allreduce(&tmp,&nparticles,1,MPI_INT,MPI_SUM,mworld);

	if (procid == 0) {
		if (screen) fprintf(screen,"  " BIGINT_FORMAT " atoms\n",nparticles);
		if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " atoms\n",nparticles);
	}

	group->update_all(nparticles);

	if (nparticles != particle->nparticles)
		error->all(FLERR,"Did not assign all atoms correctly");

	// if any atom ID < 0, error
	// if all atom IDs = 0, tag_enable = 0
	// if any atom ID > 0, error if any atom ID == 0
	// not checking if atom IDs > natoms or are unique

	int nlocal = particle->nlocal;
	int *tag = particle->tag;

	int flag = 0;
	for (int i = 0; i < nlocal; i++)
		if (tag[i] < 0) flag = 1;
	int flag_all;
	MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,mworld);
	if (flag_all)
		error->all(FLERR,"Invalid atom ID in Atoms section of data file");

	flag = 0;
	for (int i = 0; i < nlocal; i++)
		if (tag[i] > 0) flag = 1;
	MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_MAX,mworld);
	if (flag_all == 0) particle->tag_enable = 0;

	if (particle->tag_enable) {
		flag = 0;
		for (int i = 0; i < nlocal; i++)
			if (tag[i] == 0) flag = 1;
		MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,mworld);
		if (flag_all)
			error->all(FLERR,"Invalid atom ID in Atoms section of data file");
	}

  // create global mapping

	if (particle->map_style) {
		particle->map_init();
		particle->map_set();
	}
}

/* ---------------------------------------------------------------------- */

void ReadData::masses()
{
	int i,m;
	char *buf = new char[particle->ntypes*MAXLINE];
	char *original = buf;

	if (procid == 0) {
		char *eof;
		m = 0;
		for (i = 0; i < particle->ntypes; i++) {
			eof = fgets(&buf[m],MAXLINE,fp);
			if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
			m += strlen(&buf[m]);
			buf[m-1] = '\0';
		}
	}

	MPI_Bcast(&m,1,MPI_INT,0,mworld);
	MPI_Bcast(buf,m,MPI_CHAR,0,mworld);

	for (i = 0; i < particle->ntypes; i++) {
		particle->set_mass(buf);
		buf += strlen(buf) + 1;
	}
	delete [] original;
}

/* ----------------------------------------------------------------------
   read all velocities
   to find atoms, must build atom map if not a molecular system
------------------------------------------------------------------------- */

void ReadData::velocities()
{
	int i,m,nchunk;

	int mapflag = 0;
	if (particle->map_style == 0) {	
		mapflag = 1;
		particle->map_style = 1;
		particle->map_init();
		particle->map_set();
	}

	bigint nread = 0;
	bigint nparticles = particle->nparticles;

	while (nread < nparticles) {
		if (nparticles-nread > CHUNK) nchunk = CHUNK;
		else nchunk = nparticles-nread;
		if (procid == 0) {
			char *eof;
			m = 0;
			for (i = 0; i < nchunk; i++) {
				eof = fgets(&buffer[m],MAXLINE,fp);
				if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
				m += strlen(&buffer[m]);
			}
			m++;
		}
		MPI_Bcast(&m,1,MPI_INT,0,mworld);
		MPI_Bcast(buffer,m,MPI_CHAR,0,mworld);

		particle->data_vels(nchunk,buffer);
		nread += nchunk;
	}

	if (mapflag) {
		particle->map_delete();
		particle->map_style = 0;
	}

	if (procid == 0) {
		if (screen) fprintf(screen,"  " BIGINT_FORMAT " velocities\n",nparticles);
		if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " velocities\n",nparticles);
	}
}


