/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "parallel.h"
#include "domain.h"
#include "dump_custom.h"
#include "error.h"
#include "memory.h"
#include "particle.h"
#include "style_dump.h"
#include "update.h"
#include "pair_mdpd.h"
#include "force.h"

using namespace PDPS_NS;

#define MAXLINE 8192       // Maximum string length

enum{INT,DOUBLE,BIGINT};

DumpCustom::DumpCustom(PDPS *ps, int narg, char **arg) : Dump(ps, narg, arg)
{
	if (narg <= 5) error->all(FLERR,"Illegal dump custom command");

	// parse line related
	line = NULL;
	keyword = NULL;
	
	// field related
	col_labels = NULL;
	vformat = NULL;

	// output choice
	pack_choice = NULL;

	// selected particles related
	selected_list = NULL;

	max_selected = 0;

	line = new char[MAXLINE];
	line[0] = '\0';                // otherwise it will have mistakes
	for (int iarg = 5; iarg < narg; iarg++) {
	  strcat(line,arg[iarg]);
	  strcat(line," ");
	}
	line[strlen(line)-1] = '\0'; 
	
	size_one = nfields_initial = narg - 5;

	nfields = 0;
	allocate();
	parse_field(line);
}

/* ---------------------------------------------------------------------- */

DumpCustom::~DumpCustom()
{
	for (int i = 0; i < nfields; i++) {
		delete[] keyword[i];
		keyword[i] = NULL;
		delete[] vformat[i];
		vformat[i] = NULL;
	}
	delete[] keyword;
	keyword = NULL;

	delete[] vformat;
	vformat = NULL;

	delete[] line;
	line = NULL;

	delete[] vtype;
	vtype = NULL;

	delete[] pack_choice;
	pack_choice = NULL;
}

/* ---------------------------------------------------------------------- */

void DumpCustom::allocate()
{
	int n = nfields_initial;

	vformat = new char*[n];
	// store keywords
	keyword = new char*[n];
	vtype = new int[n];
	pack_choice = new FnPtrPack[n];
	//field2index = new int[n];

	for (int i = 0; i < n; i++) {
		keyword[i] = NULL;
		keyword[i] = new char[32];
        vformat[i] = NULL;
		vformat[i] = new char[32];
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::init_style()
{
	for (int i = 0; i < nfields; i++) {
		if (vtype[i] == INT) 
			strcpy(vformat[i],"%d ");
		else if (vtype[i] == DOUBLE) 
			strcpy(vformat[i],"%g ");
		else if (vtype[i] == BIGINT)
			strcpy(vformat[i],BIGINT_FORMAT);
	}

	//pack_choice = &DumpCustom::pack;
	if (column_only_flag == 0) {
		header_choice = &DumpCustom::header_item;
	}
	if (column_only_flag == 1) {
		header_choice = &DumpCustom::header_column_only;
	}
	write_choice = &DumpCustom::write_text;

	// setup boundary string
	domain->boundary_string(boundstr);

	// setup column's labels string

	int sizes = 0;
	for (int ifield = 0; ifield < nfields; ifield++) {
		sizes += strlen(keyword[ifield]);
	}
	
	col_labels = new char[sizes + nfields];
	strcpy(col_labels,keyword[0]);
	for (int ifield = 1; ifield < nfields; ifield++) {
		strcat(col_labels," ");
		strcat(col_labels,keyword[ifield]);
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::parse_field(char *str)
{
	// customize a new keyword by adding to if statement

    char *word = strtok(str," \0");
    while (word) {
		if (strcmp(word,"id") == 0) {
			addfield("ID",&DumpCustom::pack_id,INT);
			//field2index[nfield] = 0;
		} 
		else if (strcmp(word, "procid") == 0) {
			addfield("Procid", &DumpCustom::pack_procid, DOUBLE);
		}
		else if (strcmp(word,"fx") == 0) {
			addfield("Fx",&DumpCustom::pack_fx,DOUBLE);
		}
		else if (strcmp(word,"fy") == 0) {
			addfield("Fy",&DumpCustom::pack_fy,DOUBLE);
		}
		else if (strcmp(word,"fz") == 0) {
			addfield("Fz",&DumpCustom::pack_fz,DOUBLE);
		}
		else if (strcmp(word,"step") == 0) {
			addfield("Step",&DumpCustom::pack_step,INT);
			//field2index[nfield-1] = modify->find_compute_style("pressure");
		}
		else if (strcmp(word,"type") == 0) {
			addfield("Type",&DumpCustom::pack_type,INT);
			//field2index[nfield-1] = modify->find_compute_style("pressure");
		}
		else if (strcmp(word,"vx") == 0) {
			addfield("Vx",&DumpCustom::pack_vx,DOUBLE);
		} 
		else if (strcmp(word,"vy") == 0) {
			addfield("Vy",&DumpCustom::pack_vy,DOUBLE);
		} 
		else if (strcmp(word,"vz") == 0) {
			addfield("Vz",&DumpCustom::pack_vz,DOUBLE);
			//field2index[nfield-1] = modify->find_compute_style("temp");
		}
		else if (strcmp(word,"x") == 0) {
			addfield("x",&DumpCustom::pack_x,DOUBLE);
		} 
		else if (strcmp(word,"y") == 0) {
			addfield("y",&DumpCustom::pack_y,DOUBLE);
		} 
		else if (strcmp(word,"z") == 0) {
			addfield("z",&DumpCustom::pack_z,DOUBLE);
			//field2index[nfield-1] = modify->find_compute_style("temp");
		}
		else if (strcmp(word, "radius") == 0) {
	//		if (particle->radius_flag == 0) error->all(FLERR, "Illegal particle_style");
			addfield("radius", &DumpCustom::pack_radius, DOUBLE);
		}
		else if (strcmp(word, "density") == 0) {
			addfield("density", &DumpCustom::pack_density, DOUBLE);
		}
		else if (strcmp(word, "energy") == 0) {
			addfield("energy", &DumpCustom::pack_energy, DOUBLE);
		}
		else if (strcmp(word, "rho") == 0) {
			addfield("rho", &DumpCustom::pack_rho, DOUBLE);
		}
		word = strtok(NULL," \0");
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::addfield(const char *key, FnPtrPack func, int typeflag)
{
	strcpy(keyword[nfields],key);
	pack_choice[nfields] = func;
	vtype[nfields] = typeflag;
	nfields++;
}

/* ---------------------------------------------------------------------- */

int DumpCustom::count()
{
	int *mask = particle->mask;
	int i ;

	int nlocal = particle->nlocal;

	if (nlocal > max_selected) {
		max_selected = particle->nmax;
		
		memory->destroy(selected_list);
		memory->create(selected_list,max_selected,"DumpCustom: selected_list");
	}

	nselected = 0;
	for (i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit)
			selected_list[nselected++] = i;
	}

	return nselected;
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_header(bigint ndump)
{
	if (procid == 0) (this->*header_choice)(ndump);
}

/* ----------------------------------------------------------------------
   Choose a pack methdod
------------------------------------------------------------------------- */

void DumpCustom::pack(int *ids)
{
	for (int ifield = 0; ifield < nfields; ifield++) 
		(this->*pack_choice[ifield])(ifield);

	if (ids) {
		int *tag = particle->tag;
		for (int i = 0; i < nselected; i++)
			ids[i] = tag[selected_list[i]];
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_data(int n, double *mybuf)
{
	(this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_item(bigint ndump)
{
	fprintf(fp,"ITEM: TIMESTEP\n");
	fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
	fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
	fprintf(fp,"%d \n",ndump);
	fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
	fprintf(fp,"%g %g\n",boxlo[0],boxhi[0]);
	fprintf(fp,"%g %g\n",boxlo[1],boxhi[1]);
	fprintf(fp,"%g %g\n",boxlo[2],boxhi[2]);
	fprintf(fp,"ITEM: ATOMS %s\n",col_labels);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::header_column_only(bigint ndump)
{
	fprintf(fp,"%s\n",col_labels);
}

/* ---------------------------------------------------------------------- */

void DumpCustom::write_text(int n, double *mybuf)
{
	int i, j;

	int m = 0;
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < nfields; j++) {
			if (vtype[j] == INT) fprintf(fp,vformat[j],static_cast<int> (mybuf[m]));
			else if (vtype[j] == DOUBLE) fprintf(fp,vformat[j],(mybuf[m]));
			m++;
		}
		fprintf(fp,"\n");
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_id(int n)
{
	int *tag = particle->tag;
	
	for (int i = 0; i < nselected; i++) {
		buf[n] = tag[selected_list[i]];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fx(int n)
{
	double **f = particle->f;

	for (int i = 0; i < nselected; i++) {
		buf[n] = f[selected_list[i]][0];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fy(int n)
{
	double **f = particle->f;

	for (int i = 0; i < nselected; i++) {
		buf[n] = f[selected_list[i]][1];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_fz(int n)
{
	double **f = particle->f;

	for (int i = 0; i < nselected; i++) {
		buf[n] = f[selected_list[i]][2];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_step(int n)
{
	for (int i = 0; i < nselected; i++) {
		buf[n] = update->ntimestep;
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_type(int n)
{
	int *type = particle->type;
	
	for (int i = 0; i < nselected; i++) {
		buf[n] = type[selected_list[i]];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vx(int n)
{
	double **v = particle->v;

	for (int i = 0; i < nselected; i++) {
		buf[n] = v[selected_list[i]][0];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vy(int n)
{
	double **v = particle->v;

	for (int i = 0; i < nselected; i++) {
		buf[n] = v[selected_list[i]][1];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_vz(int n)
{
	double **v = particle->v;

	for (int i = 0; i < nselected; i++) {
		buf[n] = v[selected_list[i]][2];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_x(int n)
{
	double **x = particle->x;

	for (int i = 0; i < nselected; i++) {
		buf[n] = x[selected_list[i]][0];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_y(int n)
{
	double **x = particle->x;

	for (int i = 0; i < nselected; i++) {
		buf[n] = x[selected_list[i]][1];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_z(int n)
{
	double **x = particle->x;

	for (int i = 0; i < nselected; i++) {
		buf[n] = x[selected_list[i]][2];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_radius(int n)
{
	double *radius = particle->radius;

	for (int i = 0; i < nselected; i++) {
		buf[n] = radius[selected_list[i]];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_density(int n)
{
	double *rho_local = force->pair[0]->rho_local;
	double *density = particle->density;

	for (int i = 0; i < nselected; i++) {
		//buf[n] = rho_local[selected_list[i]];
		buf[n] = density[selected_list[i]];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_energy(int n)
{
	double *e = particle->e;

	for (int i = 0; i < nselected; i++) {
		buf[n] = e[selected_list[i]];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_rho(int n)
{
	double *rho = particle->rho;

	for (int i = 0; i < nselected; i++) {
		buf[n] = rho[selected_list[i]];
		n += nfields;
	}
}

/* ---------------------------------------------------------------------- */

void DumpCustom::pack_procid(int n)
{
	int procid = parallel->procid;

	for (int i = 0; i < nselected; i++) {
		buf[n] = procid;
		n += nfields;
	}
}