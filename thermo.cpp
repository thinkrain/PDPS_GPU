/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "ctype.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "compute.h"
#include "compute_temp.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "thermo.h"
#include "output.h"
#include "parallel.h"
#include "particle.h"
#include "update.h"

using namespace PDPS_NS;

//#define ONE "step temp epair emol etotal press"
//#define MULTI "etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press"

enum{SCALAR,VECTOR,ARRAY};
enum{INT,FLOAT,BIGINT};

#define INVOKED_SCALAR 1

#define MAXLINE 8192       // Maximum string length
#define DELTA 4            // For reallocation of memory, if the size of the array need to be expand

/* ---------------------------------------------------------------------- */

Thermo::Thermo(PDPS *ps, int narg, char **arg) : Pointers(ps)
{
	procid = parallel->procid;
	nprocs = parallel->nprocs;

	style = NULL;
	int n = strlen(arg[0]) + 1;
	style = new char[n];
	strcpy(style,arg[0]);

	line = NULL;
	line = new char[MAXLINE];
	if (strcmp(style,"custom") == 0){
		  //if (narg == 1) error->all(FLERR,screen,"Illegal thermo style custom command");
		line[0] = '\0';   // otherwise it will have mistakes
		for (int iarg = 1; iarg < narg; iarg++) {
		  strcat(line,arg[iarg]);
		  strcat(line," ");
		}
		line[strlen(line)-1] = '\0'; 
	}
	else error->all(FLERR, "Illegal thermo_style command");

	temperature = NULL;
	pressure = NULL;
	pe = NULL;

	keyword = NULL;
	format = NULL;

	nfield_initial = narg - 1;
	nfield = 0;

	allocate();
	parse_fields(line);

	// compute
	ncomputes = modify->ncomputes;
}

/* ---------------------------------------------------------------------- */

Thermo::~Thermo()
{
	delete [] style;
	style = NULL;
    delete [] line;
	line = NULL;
	
	int n = nfield;
	for (int i = 0; i < n; i++){
		delete [] keyword[i];
		keyword[i] = NULL;
	}
	
	delete[] keyword;
	delete[] format;
	delete[] vtype;
	delete[] vfunc;
	format = NULL;
	vtype = NULL;
	keyword = NULL;
	vfunc = NULL;

	if (temperature) delete temperature;
	temperature = NULL;
}

/* ---------------------------------------------------------------------- */

void Thermo::init()
{
	// update compute pointer
	if (modify->ncomputes > 0) computes = modify->compute;
	
	for (int i = 0; i < nfield; i++) {
		if (vtype[i] == INT) strcpy(format[i],"%8d ");
		if (vtype[i] == FLOAT) strcpy(format[i],"%12.8f ");
	}
}

/* ---------------------------------------------------------------------- */

void Thermo::allocate()
{
	int n = nfield_initial;

	format = new char*[n];
	// store keywords
	keyword = new char*[n];
	vtype = new int[n];
	vfunc = new FnPtr[n];
	field2index = new int[n];
	argindex1 = new int[n];
	argindex2 = new int[n];

	for (int i = 0; i < n; i++) {
		keyword[i] = NULL;
		keyword[i] = new char[32];
        format[i] = NULL;
		format[i] = new char[32];
		argindex1[i] = argindex2[i] = -1;
	}
}

/* ----------------------------------------------------------------------
   parse list of thermo keywords from str
   set compute flags (temp, press, pe, etc)
------------------------------------------------------------------------- */

void Thermo::parse_fields(char *str)
{
	// customize a new keyword by adding to if statement

    char *word = strtok(str," \0");
    while (word) {
		if (strcmp(word,"step") == 0) {
			field2index[nfield] = 0;
			addfield("Step",&Thermo::compute_step,INT);
		} else if (strcmp(word,"temp") == 0) {
			if (temperature == NULL) {
				char **arg = new char*[3];
				arg[0] = (char *) "velocity_temp";
				arg[1] = (char *) "all";             // always group all
				arg[2] = (char *) "temp";
				temperature = new ComputeTemp(ps, 3, arg);
				delete [] arg;
				arg = NULL;
			}
		    temperature->init();
			field2index[nfield] = modify->find_compute_style("temp");
			addfield("Temp",&Thermo::compute_temp,FLOAT);
		}
		else if (strcmp(word,"press") == 0) {
			field2index[nfield] = modify->find_compute_style("pressure");
			if (field2index[nfield] == -1) error->all(FLERR, "Cannot find the compute pressure style");
			addfield("Press",&Thermo::compute_compute,FLOAT);
		}
		else if ((strncmp(word, "c_", 2) == 0)) {
			// find compute by name
			int n = strlen(word);
			char *id = new char[n-2];
			strcpy(id, &word[2]);
			//char *copy = new char[n-2];
			//strcpy(copy, id);

			char *ptr = strchr(id, '[');
			if (ptr != NULL) {
				*ptr = '\0';
				argindex1[nfield] = int_between_brackets(ptr);
				ptr++;
				if (*ptr == '[') {
					argindex2[nfield] = int_between_brackets(ptr);
					ptr++;
				}
				else if (*ptr != 0) error->all(FLERR, "Illegal brackets");
			}
			int cid = modify->find_compute(id);
			if (cid < 0) error->all(FLERR, "Cannot find thermo custom compute ID");
			if (argindex1[nfield] < 0 && modify->compute[cid]->scalar_flag == 0) {
				error->all(FLERR, "Thermo compute does not compute scalar");
			}
			if (argindex1[nfield] >= 0 && argindex2[nfield] < 0) {
				if (modify->compute[cid]->vector_flag == 0) {
					error->all(FLERR, "Thermo compute does not compute vector");
				}
				if (argindex1[nfield] > modify->compute[cid]->size_vector) {
					error->all(FLERR, "Thermo compute vector exceeds the upper limit");
				}
			}
			if (argindex1[nfield] >=0 && argindex2[nfield] >= 0) {
				if (modify->compute[cid]->array_flag == 0) {
					error->all(FLERR, "Thermo compute does not compute array");
				}
				if (argindex1[nfield] > modify->compute[cid]->size_array_rows || 
					argindex2[nfield] > modify->compute[cid]->size_array_columns) {
					error->all(FLERR, "Thermo compute array exceeds the upper limit");
				}
			}
			field2index[nfield] = cid;
			addfield(word, &Thermo::compute_compute, FLOAT);
		}
		else error->all(FLERR, "Illegal thermo style");
		word = strtok(NULL," \0");
	}
}

/* ----------------------------------------------------------------------
   Find the integer inside a bracket: eg. [32]
------------------------------------------------------------------------- */

int Thermo::int_between_brackets(char *&ptr)
{
	char *start = ++ptr;

	while (*ptr && *ptr != ']') {
		if (!isdigit(*ptr)) {
			error->all(FLERR,"Non digit character between brackets in variable");
		}
		ptr++;
	}

	if (*ptr != ']') error->all(FLERR, "Mismatched brackets in variable");
	if (ptr == start) error->all(FLERR, "Empty brackets in variable");

	*ptr = '\0';
	int index = atoi(start);
	*ptr = ']';

	return index;
}

/* ----------------------------------------------------------------------
   add field to list of quantities to print
------------------------------------------------------------------------- */

void Thermo::addfield(const char *key, FnPtr func, int typeflag)
{
	strcpy(keyword[nfield],key);
	vfunc[nfield] = func;
	vtype[nfield] = typeflag;
	nfield++;
}

/* ---------------------------------------------------------------------- */

void Thermo::header()
{
	int loc = 0;
	for (int i = 0; i < nfield; i++)
	{
		loc += sprintf(&line[loc],"%s ",keyword[i]);
	}
	sprintf(&line[loc],"\n");
  
	if (procid == 0) {
		if (screen) fprintf(screen,"%s",line);
		if (logfile) fprintf(logfile,"%s",line);
	}
}

/* ---------------------------------------------------------------------- */

void Thermo::compute()
{
	int i;

	//lost_check();
	// invoke Compute methods needed for thermo keywords

	for (i = 0; i < ncomputes; i++){
		if (computes[i]->scalar_flag == 1) {
			if (computes[i]->invoked_scalar != update->ntimestep) {
				computes[i]->compute_scalar();
			}
		}
		if (computes[i]->vector_flag == 1) {
			if (computes[i]->invoked_vector != update->ntimestep) {
				computes[i]->compute_vector();
			}
		}
		if (computes[i]->array_flag == 1) {
			if (computes[i]->invoked_array != update->ntimestep) {
				computes[i]->compute_array();
			}
		}
	}
    
	int loc = 0;
    for (ifield = 0; ifield < nfield; ifield++) {
		(this->*vfunc[ifield])();
		if (vtype[ifield] == INT) {
			loc += sprintf(&line[loc], format[ifield], bivalue);
		}
		if (vtype[ifield] == FLOAT) {
			loc += sprintf(&line[loc], format[ifield], dvalue);
		}
	}

	if (procid == 0) {
		if (screen) fprintf(screen,"%s\n",line);
		if (logfile) {
		  fprintf(logfile,"%s\n",line);
		}
		fflush(logfile);
	}
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
   compute/fix are normalized by atoms if returning extensive value
   variable value is not normalized (formula should normalize if desired)
------------------------------------------------------------------------- */

void Thermo::compute_compute()
{
	int m = field2index[ifield];

	Compute *compute = computes[m];

	if (compute->scalar_flag == 1) {
		dvalue = compute->scalar;
	} 
	else if (compute->vector_flag == 1) {
		dvalue = compute->vector[argindex1[ifield]];
	}
	else if (compute->array_flag == 1) {
		dvalue = compute->array[argindex1[ifield]][argindex2[ifield]];
	}
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_step()
{
	bivalue = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void Thermo::compute_temp()
{
	dvalue = temperature->compute_scalar();
}

/* ---------------------------------------------------------------------- */

void Thermo::lost_check()
{
	int ntotal;
	int nlocal = particle->nlocal;

	MPI_Allreduce(&nlocal, &ntotal, 1, MPI_INT, MPI_SUM, mworld);

	if (ntotal < particle->nparticles) {
		char str[128];
		sprintf(str, "Particle lost from %d total particles to %d", particle->nparticles, ntotal);
		error->all(FLERR, str);
	}
}
