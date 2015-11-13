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
#include "error.h"
#include "fix_ave_spatial.h"
#include "particle.h"
#include "region.h"
#include "style_fix.h"
#include "update.h"

using namespace PDPS_NS;
using namespace FixConst;

#define MAXLINE 8192       // Maximum string length

enum{INT,FLOAT,BIGINT};

/* ---------------------------------------------------------------------- */

FixAveSpatial::FixAveSpatial(PDPS *ps, int narg, char **arg) : Fix(ps, narg, arg)
{
	int n;
	if (narg < 10) error->all(FLERR,"Illegal fix ave/spatial command");

	// initilialize important variabels
	invoked_flag = 1;
	average_flag = 0;
	iter = 0;
	counter = 0;

	// initialize pointers
	bivalues_bin = NULL;
	bivalues_t = NULL;
	dvalues_bin = NULL;
	dvalues_t = NULL;
	field2index = NULL;
	file = NULL;
	fname = NULL;
	format = NULL;
	header = NULL;
	hx = NULL;
	hy = NULL;
	hz = NULL;
	keyword = NULL;
	line = NULL;
	num_bins = NULL;
	rname = NULL;
	vtype = NULL;
	vfunc = NULL;

	g_start = atoi(arg[3]);         // global start timestep
	g_end = atoi(arg[4]);           // global end timestep
	nevery = atoi(arg[5]);          
	nrepeat = atoi(arg[6]);
	nfreq = atoi(arg[7]);

	// region's name 
	if (strcmp(arg[8],"region") == 0) {
		n = strlen(arg[9]) + 1;
		rname = new char[n];
		strcpy(rname,arg[9]);
		rid = domain->find_region(rname);
		if (rid == -1) {
			char str[128];
			sprintf(str,"Cannot find the region %s", rname);
			error->all(FLERR,str);
		}
	}
	
	// initialization
	for (int i = 0; i < 3; i++) {
		dim[i] = 0;
		nbins[i] = 0;
		dh[i] =0.0;
	}

	// create bins along which direction, so far we only have one direction available
	// needs to be more general
	if (strcmp(arg[10],"x") == 0) {
		dim[0] = 1;
		nbins[0] = atoi(arg[11]);
	}
	else if (strcmp(arg[10],"y") == 0) {
		dim[1] = 1;
		nbins[1] = atoi(arg[11]);
	}
	else if (strcmp(arg[10],"z") == 0) {
		dim[2] = 1;
		nbins[2] = atoi(arg[11]);
	}
	
	// restore all keywords into "line"
	line = new char[MAXLINE];
	line[0] = '\0';                // otherwise it will have mistakes
	for (int iarg = 12; iarg < narg - 1; iarg++) {
	  strcat(line,arg[iarg]);
	  strcat(line," ");
	}
	line[strlen(line)-1] = '\0';
	n = strlen(line) + 1;
	header = new char[n];
	strcpy(header,line);

	// create output file
	n = strlen(arg[narg-1]) + 1;
	fname = new char[n];
	strcpy(fname,arg[narg-1]);
	file = fopen(fname,"w");
	if (file == NULL) {
		char str[128];
		sprintf(str,"Cannot open file \"s\"",fname);
		error->all(FLERR,str);
	}

	// number of keywords
	nfield_initial = narg - 13;
	nfield = 0;
	nbivalue = 0;
	ndvalue = 0;

	allocate();
	parse_field(line);
	
	// allocate bivalues for each bin if it exist
	if (nbivalue > 0) {
		bivalues_bin = new int*[nbivalue];
		bivalues_t = new int*[nbivalue];
		for (int i = 0; i < nbivalue; i++) {
			bivalues_bin[i] = NULL;
			bivalues_bin[i] = new int[nbins[2]];
			bivalues_t[i] = NULL;
			bivalues_t[i] = new int[nbins[2]];
		}
	}
	// allocate dvalues for each bin if it exist
	if (ndvalue > 0) {
		dvalues_bin = new double*[ndvalue];
		dvalues_t = new double*[ndvalue];
		for (int i = 0; i < ndvalue; i++) {
			dvalues_bin[i] = NULL;
			dvalues_bin[i] = new double[nbins[2]];
			dvalues_t[i] = NULL;
			dvalues_t[i] = new double[nbins[2]];
		}
	}
}

/* ---------------------------------------------------------------------- */

FixAveSpatial::~FixAveSpatial()
{
	delete fname;
	fname = NULL;
	delete [] rname;
	rname = NULL;

	if (file) fclose(file);
	file = NULL;

	if (dim[0] == 1) {
		delete [] hx;
		hx = NULL;
	}
	if (dim[1] == 1) {
		delete [] hy;
		hy = NULL;
	}
	if (dim[2] == 1) {
		delete [] hz;
		hz = NULL;
	}

	for (ifield = 0; ifield < 3; ifield++) {
		delete [] format[ifield];
		format[ifield] = NULL;
		delete [] keyword[ifield];
		keyword[ifield] = NULL;
		int id = field2index[ifield];
		if (vtype[ifield] == INT) {
			delete [] bivalues_bin[id];
			bivalues_bin[id] = NULL;
			delete [] bivalues_t[id];
			bivalues_t[id] = NULL;
		}
		if (vtype[ifield] == FLOAT) {
			delete [] dvalues_bin[id];
			dvalues_bin[id] = NULL;
			delete [] dvalues_t[id];
			dvalues_t[id] = NULL;
		}
	}
	delete [] format;
	format = NULL;
	delete [] keyword;
	keyword = NULL;
	
	if (nbivalue > 0) {
		delete [] bivalues_bin;
		bivalues_bin = NULL;
		delete [] bivalues_t;
		bivalues_t = NULL;
	}
	if (ndvalue > 0) {
		delete [] dvalues_bin;
		dvalues_bin = NULL;
		delete [] dvalues_t;
		dvalues_t = NULL;
	}

	delete [] line;
	line = NULL;
	delete [] header;
	header = NULL;
	
	delete [] num_bins;
	num_bins = NULL;

	delete [] vfunc;
	vfunc = NULL;

	delete [] vtype;
	vtype = NULL;

	delete [] field2index;
	field2index = NULL;
}

/* ---------------------------------------------------------------------- */

int FixAveSpatial::setmask()
{
	int mask = 0;
	mask |= END_OF_STEP;
	return mask;
}

/* ----------------------------------------------------------------------
   init
------------------------------------------------------------------------- */

void FixAveSpatial::init()
{
	// store region's information
	xlo = domain->regions[rid]->extent_xlo;
	xhi = domain->regions[rid]->extent_xlo;
	xle = domain->regions[rid]->extent_xle;
	ylo = domain->regions[rid]->extent_ylo;
	yhi = domain->regions[rid]->extent_yhi;
	yle = domain->regions[rid]->extent_yle;
	zlo = domain->regions[rid]->extent_zlo;
	zhi = domain->regions[rid]->extent_zhi;
	zle = domain->regions[rid]->extent_zle;

	// compute coordinate of each bin
	if (dim[0] == 1) {
		dh[0] = xle/nbins[0];
		for (int j = 0; j < nbins[0]; j++) {
			hx[j] = xlo + 0.5*dh[0] + dh[0]*j;	
		}	
	}
	if (dim[1] == 1) {
		dh[1] = yle/nbins[1];
		for (int j = 0; j < nbins[1]; j++) {
			hy[j] = ylo + 0.5*dh[1] + dh[1]*j;	
		}	
	}
	if (dim[2] == 1) {
		dh[2] = zle/nbins[2];
		for (int j = 0; j < nbins[2]; j++) {
			hz[j] = zlo + 0.5*dh[2] + dh[2]*j;	
		}	
	}

	// number of bins created along z coordinate
	// needs to be more general
	
	for (int i = 0; i < nbins[2]; i++) {
		num_bins[i] = 0;
	}

	// record format used for output
	for (int i = 0; i < nfield; i++) {
		if (vtype[i] == INT) 
			strcpy(format[i],"%8d ");
		if (vtype[i] == FLOAT) 
			strcpy(format[i],"%12.8f ");
	}
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::allocate()
{
	int n = nfield_initial;

	format = new char*[n];
	// store keywords
	keyword = new char*[n];
	vtype = new int[n];
	vfunc = new FnPtr[n];
	field2index = new int[n];

	for (int i = 0; i < n; i++) {
		keyword[i] = NULL;
		keyword[i] = new char[32];
        format[i] = NULL;
		format[i] = new char[32];
	}

	if (dim[0] == 1) {
		hx = new double[nbins[0]];
	}
	if (dim[1] == 1) {
		hy = new double[nbins[1]];
	}
	if (dim[2] == 1) {
		hz = new double[nbins[2]];
	}
	num_bins = new int[nbins[2]];
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::parse_field(char *str)
{
	// customize a new keyword by adding to if statement

    char *word = strtok(str," \0");
    while (word) {
		if (strcmp(word,"vx") == 0) {
			field2index[nfield] = ndvalue;
			addfield("Vx",&FixAveSpatial::write_vx,FLOAT);
			ndvalue++;
		} 
		else if (strcmp(word,"vy") == 0) {
			field2index[nfield] = ndvalue;
			addfield("Vy",&FixAveSpatial::write_vy,FLOAT);
			ndvalue++;
		}
		else if (strcmp(word,"vz") == 0) {
			field2index[nfield] = ndvalue;
			addfield("Vz",&FixAveSpatial::write_vz,FLOAT);
			ndvalue++;
		}
		word = strtok(NULL," \0");
	}
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::addfield(const char *key, FnPtr func, int typeflag)
{
	strcpy(keyword[nfield],key);
	vfunc[nfield] = func;
	vtype[nfield] = typeflag;
	nfield++;
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::end_of_step()
{
	int ntimestep = update->ntimestep;
	int n;

	// global start and end timestep
	if (ntimestep >= g_start && ntimestep <= g_end) {
		n = (ntimestep - g_start)/nfreq;

		// define local start and end timestep
		if (average_flag == 0) {
			n_start = g_start + (n + 1)*nfreq - (nrepeat - 1)*nevery;
			n_end = g_start + (n + 1)*nfreq;
			average_flag = 1;
			iter = 0;
			counter++;
		}
		// local end timestep cannot be larger than the global end timestep
		if (n_end > g_end) {
			n_end = g_end;
		}
		
		if (ntimestep == n_start) {
			clear();
		}
		if (ntimestep >= n_start && ntimestep <= n_end) {
			ago = ntimestep - n_start;
			if (ago % nevery == 0) {
				average();
			}
		}
	}
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::average()
{
	int ntimestep = update->ntimestep;
	int *mask = particle->mask;
	int nparticles = particle->nparticles;
	int id, ib;

	iter++;    // timesteps used in one average
	for (ibin = 0; ibin < nbins[2]; ibin++) {
		num_bins[ibin] = 0;
	}
	// Compute average
	for (iparticle = 0; iparticle < nparticles; iparticle++) {
		if (mask[iparticle] & groupbit) {
			ib = check_bins(iparticle);    // locate bin;
			if (ib < 0 || ib > (nbins[2]-1)) {
				error->all(FLERR,"Invalid ibin\n");
			}
			num_bins[ib]++;                // calculate number of particles in each bin
			for (ifield = 0; ifield < nfield; ifield++) {
				(this->*vfunc[ifield])();
				if (vtype[ifield] == FLOAT) {
					id = field2index[ifield];
					dvalues_bin[id][ib] += dvalue;   // sum all dvalues in each bin
				}
			}
		}
	} // loop for all particles
	// ensemble average: get average for dvalues_bin and sum it for dvalues_t
	ensemble_average();

	// time average
	if (ntimestep == n_end) {
		time_average();
		write_header();
		write();
		average_flag = 0;
		invoked_flag = ntimestep;
		//vector[0] = dvalues_
	}
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::write_header()
{
	char str[128];
	int ntimestep = update->ntimestep;

	fprintf(file, "ITEM: AVERAGE COUNTER\n");
	fprintf(file, "%d \n", counter);
	fprintf(file, "ITEM: TIMESTEPS SPAN OF THIS AVERAGE\n");
	fprintf(file, " %d %d\n", n_start, n_end);
	fprintf(file, "ITEM: NUMBER OF BINS\n");
	fprintf(file, "  %d \n", nbins[2]);
	sprintf(str, "bin_x bin_y bin_z %s\n", header);
	fprintf(file,str);
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::write()
{
	int loc = 0;
	int id;

	for (ibin = 0; ibin < nbins[2]; ibin++) {
		loc = 0;
		loc += sprintf(&line[loc],"%12.8f ",hz[ibin]);
		for (ifield = 0; ifield < nfield; ifield++) {
				id = field2index[ifield];
			if (vtype[ifield] == FLOAT) {
				loc += sprintf(&line[loc],format[ifield],dvalues_t[id][ibin]);
			}
		}
		fprintf(file,"%s\n",line);
	}
	fflush(file);
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::clear()
{
	int id;

	for (ifield = 0; ifield < nfield; ifield++) {
		id = field2index[ifield];
		for (ibin = 0; ibin < nbins[2]; ibin++) {	
			if (vtype[ifield] == INT) {
				bivalues_bin[id][ibin] = 0;
				bivalues_t[id][ibin] = 0;
			}
			if (vtype[ifield] == FLOAT) {
				dvalues_bin[id][ibin] = 0.0;
				dvalues_t[id][ibin] = 0.0;
			}
		}
	}
}

/* ----------------------------------------------------------------------
   locate bin id for each particle
------------------------------------------------------------------------- */

int FixAveSpatial::check_bins(int id) 
{
	double **x = particle->x;
	int icell;

	if (dim[2] == 1) {
		icell = int((x[id][2] - zlo)/dh[2]);
		if (icell < 0) icell = 0;
		if (icell >= nbins[2]) icell = nbins[2] - 1;
	}

	return icell;
}

/* ----------------------------------------------------------------------
   Ensemble average: average dvalues_bin, add it for dvalues_t 
                     and zero dvalues_bin
------------------------------------------------------------------------- */

void FixAveSpatial::ensemble_average()
{
	int id;

	for (ifield = 0; ifield < nfield; ifield++) {
		id = field2index[ifield];
		for (ibin = 0; ibin < nbins[2]; ibin++) {
			if (num_bins[ibin]!= 0) {
				if (vtype[ifield] == FLOAT) {
					dvalues_t[id][ibin] += dvalues_bin[id][ibin]/num_bins[ibin];
					dvalues_bin[id][ibin] = 0.0;
				}
			}
		}
	}
}

/* ----------------------------------------------------------------------
   Time average: average dvalues_t 
------------------------------------------------------------------------- */

void FixAveSpatial::time_average()
{
	int id;

	for (ifield = 0; ifield < nfield; ifield++) {
		id = field2index[ifield];
		for (ibin = 0; ibin < nbins[2]; ibin++) {
			if (vtype[ifield] == FLOAT) {
				dvalues_t[id][ibin] /= iter;
			}
		}
	}
}	

/* ---------------------------------------------------------------------- */

void FixAveSpatial::write_vx()
{
	double **v = particle->v;

	dvalue = v[iparticle][0];
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::write_vy()
{
	double **v = particle->v;

	dvalue = v[iparticle][1];
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::write_vz()
{
	double **v = particle->v;

	dvalue = v[iparticle][2];
}
