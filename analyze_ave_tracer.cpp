/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include <limits>
#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "analyze_ave_tracer.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "region.h"
#include "pair.h"
#include "parallel.h"
#include "particle.h"
#include "update.h"

using namespace PDPS_NS;

enum{X, V, MAX_V, F, STEP, COMPUTE};
enum{INT, FLOAT};
enum{SCALAR, VECTOR, ARRAY};

/* ---------------------------------------------------------------------- */

AnalyzeAveTracer::AnalyzeAveTracer(PDPS *ps, int narg, char **arg) : Analyze(ps, narg, arg)
{
	if (narg < 7) error->all(FLERR, "Illegal analyze ave/tracer command");

	nstart = 0;
	ave_flag = 0;
	clean_ave_flag = 0;

	scalar_flag = vector_flag = 0;
	array_flag = 1;
	nrows = ncolumns = 0;
	
	max_ntags = 0;

	header_flag = 1;
	header_once_flag = 0;

	fname = NULL;
	file = NULL;
	field_func = NULL;

	array = NULL;
	array_total = NULL;
	array_ave = NULL;

	nx = ny = nz = 0;
	ntracers = 0;
	tag2tracerID_local = NULL;
	tag2tracerID = NULL;
	num_cell = NULL;
	numAll_cell = NULL;

	nevery = atoi(arg[3]);
	nrepeat = atoi(arg[4]);
	nfreq = atoi(arg[5]);
	rid = domain->find_region(arg[6]);
	if (rid == -1) error->all(FLERR, "Cannot find the region id");

	// initialization of dim[3] and originflag[3]
	ndims = 0;
	for (int i = 0; i < 3; i++) {
		dim[i] = -1;
	}

	int iarg = 7;
	// parse spatial cell information
	while (iarg < narg) {
		if (iarg + 3 > narg) break;
		if (strcmp(arg[iarg], "x") == 0) dim[ndims] = 0;
		else if (strcmp(arg[iarg], "y") == 0) dim[ndims] = 1;
		else if (strcmp(arg[iarg], "z") == 0) dim[ndims] = 2;
		else break;

		if (dim[ndims] == 2 && domain->dim == 2) {
			error->all(FLERR, "Cannot use analyze ave/spatial z for 2D model");
		}

		delta[ndims] = atof(arg[iarg+1]);

		ndims++;
		iarg += 2;
		if (ndims == domain->dim) break;
	}

	// define other delta
	double *boxle = domain->boxle;
	// 1D
	if (ndims == 1) {
		if (dim[0] == 0) {
			delta[1] = boxle[1];
			delta[2] = boxle[2];
		}
		else if (dim[0] == 1) {
			delta[1] = boxle[0];
			delta[2] = boxle[2];
		}
		else if (dim[0] == 2) {
			delta[1] = boxle[0];
			delta[2] = boxle[1];
		}
	}
	// 2D
	else if (ndims == 2) {
		if ((dim[0] + dim[1]) == 1) delta[2] = boxle[2];
		else if ((dim[0] + dim[1]) == 2) delta[2] = boxle[1];
		else if ((dim[0] + dim[1]) == 3) delta[2] = boxle[0];
	}
	for (int i = 0; i < 3; i++)	{
		inv_delta[i] = 1.0 / delta[i];
	}

	if (ndims == 0) error->all(FLERR, "Illegal analyze ave/space command");
	if (ndims == 2 && dim[0] == dim[1]) {
		error->all(FLERR, "Same dimension twice in analyze ave/space command");
	}
	if (ndims == 3 && (dim[0] == dim[1] || dim[1] == dim[2] || dim[0] == dim[2])) {
		error->all(FLERR, "Same dimension twice in analyze ave/space");
	}

	             		
	nfields_initial = narg - 6;
	allocate();

	nfields = 0;
	compute = modify->compute;
	parse_fields(iarg, narg, arg);
	iarg += nfields;
	
	// parse options
	while (iarg < narg) {
		if (strcmp(arg[iarg], "file") == 0) {
			int n = strlen(arg[iarg + 1]) + 1;
			fname = new char[n];
			strcpy(fname, arg[iarg + 1]);
			if (procid == 0) {
				file = fopen(fname, "w");
				if (file == NULL) {
					char str[128];
					sprintf(str, "Cannot open analyze ave/space file %s", arg[iarg + 1]);
					error->one(FLERR, str);
				}
			}
			iarg += 2;
		}
		else if (!strcmp(arg[iarg], "header")) {
			if (!strcmp(arg[iarg + 1], "no")) header_flag = 0;
			else if (!strcmp(arg[iarg + 1], "yes")) header_flag = 1;
			else if (!strcmp(arg[iarg + 1], "once")) header_once_flag = 1;
			else error->all(FLERR, "Illegal analyze ave/time header option");
			iarg += 2;
		}
		else if (!strcmp(arg[iarg], "start")) {
			nstart = atoi(arg[iarg + 1]);
			iarg += 2;
		}
		else error->all(FLERR, "Illegal analyz ave/space filed name or option");
	}

	if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0) {
		error->all(FLERR, "Illegal analyze ave/space average operation parameters");
	}
	if (nfreq % nevery || (nrepeat - 1)*nevery >= nfreq) {
		error->all(FLERR, "Illegal analyze ave/space average operation parameters");
	}

	for (int i = 0; i < nfields; i++) {
		ncolumns += field_ncols[i];
	}

	setup_tracers();
	assign_tracerID();
}

/* ---------------------------------------------------------------------- */

AnalyzeAveTracer::~AnalyzeAveTracer()
{
	memory->destroy(tag2tracerID);
	memory->destroy(numAll_cell);
	
	if (array_flag == 1) {
		memory->destroy(array);
		memory->destroy(array_total);
		memory->destroy(array_ave);
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTracer::addfield(const char *key, FnPtr func, int typeflag, int data_typeflag)
{
	strcpy(field_name[nfields], key);
	field_func[nfields] = func;
	field_data_type[nfields] = data_typeflag;
	field_type[nfields] = typeflag;

	if (typeflag == SCALAR) {
		nrows = 1;
		field_ncols[nfields] = 1;
	}
	else if (typeflag == VECTOR) {
		nrows = 1;
	}

	nfields++;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTracer::allocate()
{
	int n = nfields_initial;

	field_format = new char*[n];
	// store field_names
	field_name = new char*[n];
	field_data_type = new int[n];
	field_type = new int[n];
	field_func = new FnPtr[n];
	field_ncols = new int[n];
	field_index = new int[n];

	for (int i = 0; i < n; i++) {
		field_name[i] = NULL;
		field_name[i] = new char[32];
		field_format[i] = NULL;
		field_format[i] = new char[32];
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTracer::init()
{
	// update compute pointer
	if (modify->ncomputes > 0) compute = modify->compute;

	for (int i = 0; i < nfields; i++) {
		if (field_data_type[i] == INT) strcpy(field_format[i], "%8d ");
		else if (field_data_type[i] == FLOAT) strcpy(field_format[i], "%12.8f ");
	}

	// update tag2tracerID
	/*if (max_ntags < particle->nparticles+1) {
		max_ntags = particle->nparticles + 1;
		memory->grow(tag2tracerID, max_ntags, "AnalyzeAveTracer: tag2tracerID");
	}*/
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTracer::setup_tracers()
{
	int i, j, k, idim;
	double area, vol;
	double *extent_lo, *extent_hi, *extent_le;

	extent_lo = domain->regions[rid]->extent_lo;
	extent_hi = domain->regions[rid]->extent_hi;
	extent_le = domain->regions[rid]->extent_le;

	// define cell origin
	ntracers = 1;
	for (i = 0; i < 3; i++) tracer[i] = 1;
	for (i = 0; i < ndims; i++) {
		tracer[i] = static_cast<int> (extent_le[dim[i]] * inv_delta[i]);
		if (tracer[i] == 0) tracer[i] = 1;
		ntracers *= tracer[i];
	}

	nrows = ntracers;
	
	memory->grow(array, nrows, ncolumns, "AnalyzeAveSpace: array");
	memory->grow(array_total, nrows, ncolumns, "AnalyzeAveSpace: array_total");
	memory->grow(array_ave, nrows, ncolumns, "AnalyzeAveSpace: array_ave");
	max_ntags = particle->nparticles + 1;
	memory->grow(tag2tracerID_local, max_ntags, "AnalyzeAveTracer: tag2tracerID_local");
	memory->grow(tag2tracerID, max_ntags, "AnalyzeAveTracer: tag2tracerID"); 
	memory->grow(num_cell, nrows, "AnalyzeAveSpace: num_cell");
	memory->grow(numAll_cell, nrows, "AnalyzeAveSpace: numAll_cell");

	for (int i = 0; i < nrows; i++) {
		num_cell[i] = numAll_cell[i] = 0;
	}
	for (int i = 0; i < max_ntags; i++) {
		tag2tracerID_local[i] = -1;
		tag2tracerID[i] = -1;
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTracer::assign_tracerID()
{
	int inside_flag, tracer_id;
	int *mask = particle->mask;
	int *tag = particle->tag;
	double **x = particle->x;
	double *extent_lo, *extent_hi, *extent_le;

	extent_lo = domain->regions[rid]->extent_lo;
	extent_hi = domain->regions[rid]->extent_hi;
	extent_le = domain->regions[rid]->extent_le;

	for (int i = 0; i < particle->nlocal; i++) {
		if (mask[i] & groupbit) {
			inside_flag = domain->regions[rid]->inside(x[i]);
			if (inside_flag == 0) continue;
			int c[3], idim;
			for (int j = 0; j < ndims; j++) {
				idim = dim[j];
				c[j] = static_cast<int> ((x[i][idim] - extent_lo[idim]) / delta[j]);
				// the last cell may have large size
				if (c[j] == tracer[j]) c[j] = tracer[j] - 1;
			}
			tracer_id = c[2] * tracer[1] * tracer[0] + c[1] * tracer[0] + c[0];
			tag2tracerID_local[tag[i]] = tracer_id;
			num_cell[tracer_id]++;
		}
	}

	MPI_Allreduce(&tag2tracerID_local[0], &tag2tracerID[0], max_ntags, MPI_INT, MPI_MAX, mworld);
	MPI_Allreduce(&num_cell[0], &numAll_cell[0], ntracers, MPI_INT, MPI_SUM, mworld);

	memory->destroy(tag2tracerID_local);
	memory->destroy(num_cell);
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTracer::invoke_analyze() 
{
	int ntimestep = update->ntimestep;

	int mod = ntimestep % nevery;

	int temp, start, end;

	if (ntimestep % nfreq == 0) temp = ntimestep / nfreq;
	else temp = ntimestep / nfreq + 1;

	// start and end steps
	end = temp * nfreq;
	start = end - (nrepeat - 1) * nevery;

	// example: nevery = 10 nrepeat = 5 nfreq = 1000
	// time averages are calculated on 960 970 980 990 1000 steps
	if ((end - ntimestep) % nevery == 0 && ntimestep >= start && ntimestep >= nstart) {
		if (ntimestep == start) clean_ave_flag = 1;
		else clean_ave_flag = 0;
		if (ntimestep == end) ave_flag = 1; // do average at this timestep
		else ave_flag = 0;

		clean_array();

		// compute
		icol = 0;
		for (ifield = 0; ifield < nfields; ifield++) {
			(this->*field_func[ifield])();
			icol += field_ncols[ifield];
		}

		MPI_Allreduce(&array[0][0], &array_total[0][0], nrows*ncolumns, MPI_DOUBLE, MPI_SUM, mworld);

		icol = 0;
		for (ifield = 0; ifield < nfields; ifield++) {
			for (int itracer = 0; itracer < ntracers; itracer++) {
				for (int j = 0; j < field_ncols[ifield]; j++) {
						if (numAll_cell[itracer] > 0) {
							array_total[itracer][icol + j] /= numAll_cell[itracer];
						}
					}
				}
			}

		// sum all for average
		for (int i = 0; i < nrows; i++)
		for (int j = 0; j < ncolumns; j++) {
			array_ave[i][j] += array_total[i][j];
		}

		// average
		if (ave_flag == 1) {
			for (int i = 0; i < nrows; i++) 
			for (int j = 0; j < ncolumns; j++) {
				array_ave[i][j] = array_ave[i][j] / nrepeat;
			}
			write_array();
		}
	}	
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTracer::parse_fields(int iarg, int narg, char **arg)
{
	while (iarg < narg) {
		if (strcmp(arg[iarg], "centroid") == 0) {
			field_ncols[nfields] = 3;
			addfield("centroid", &AnalyzeAveTracer::compute_centroid, ARRAY, FLOAT);
			iarg++;
		}
		else break;
	} // while (iarg < nargs)
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTracer::write_array()
{
	int i, j;

	// when output file is specified
	if (file != NULL) {
		// master processor
		if (procid == 0) {
			if (header_flag == 1) {
				fprintf(file, "TimeStep NumberOfRows\n");
				fprintf(file, "%d %d\n", update->ntimestep, nrows);
				fprintf(file, "Row ");
				for (ifield = 0; ifield < nfields; ifield++) {
					if (field_type[ifield] == SCALAR) fprintf(file, "%s ", field_name[ifield]);
					else if (field_type[ifield] == VECTOR || field_type[ifield] == ARRAY){
						for (j = 0; j < field_ncols[ifield]; j++) {
							fprintf(file, "%s[%d] ", field_name[ifield], j);
						}
					}
				}
				fprintf(file,"\n");
				if (header_once_flag == 1) header_flag = 0;
			}
			int temp = field_ncols[0];
			int curr_field = 0;
			int ivalue;
			double dvalue;
			for (i = 0; i < nrows; i++) {
				fprintf(file, "%d ", i+1); 
				for (j = 0; j < ncolumns; j++) {
					if (j == temp) {
						curr_field++;
						temp += field_ncols[curr_field];
					}
					if (field_data_type[curr_field] == INT) {
						fprintf(file, field_format[curr_field], (int)array_ave[i][j]);
					}
					else if (field_data_type[curr_field] == FLOAT) {
						fprintf(file, field_format[curr_field], array_ave[i][j]);
					}
				}
				fprintf(file, "\n");
			}
		} // if (procid == 0)
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTracer::clean_array()
{
	for (int i = 0; i < nrows; i++) 
	for (int j = 0; j < ncolumns; j++) {
		array[i][j] = 0.0;
		array_total[i][j] = 0.0;
		if (clean_ave_flag) {
			array_ave[i][j] = 0.0;
		}
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTracer::compute_centroid()
{
	int *mask = particle->mask;
	int *tag = particle->tag;
	double **x = particle->x;
	int nlocal = particle->nlocal;
	int tracer_id;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			tracer_id = tag2tracerID[tag[i]];
			if (tracer_id < 0) {
				error->all(FLERR, "Illegal tracer ID. Check the code\n");
			}
			array[tracer_id][icol] += x[i][0];
			array[tracer_id][icol+1] += x[i][1];
			array[tracer_id][icol+2] += x[i][2];
		}
	}	
}
