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

#include "analyze_ave_time.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair.h"
#include "parallel.h"
#include "particle.h"
#include "update.h"

using namespace PDPS_NS;

enum{INT, FLOAT};
enum{SCALAR, VECTOR, ARRAY};
#define EPSILON 1.0e-10
#define PI 3.1416

/* ---------------------------------------------------------------------- */

AnalyzeAveTime::AnalyzeAveTime(PDPS *ps, int narg, char **arg) : Analyze(ps, narg, arg)
{
	if (narg < 7) error->all(FLERR, "Illegal analyze ave/time command");

	nstart = 0;
	ave_flag = 0;
	clean_ave_flag = 0;

	scalar_flag = vector_flag = array_flag = 0;
	nrows = ncolumns = 0;
	
	header_flag = 1;
	header_once_flag = 0;

	field_func = NULL;

	nevery = atoi(arg[3]);
	nrepeat = atoi(arg[4]);
	nfreq = atoi(arg[5]);

	// parse field_name
	nfields_initial = narg - 6;
	allocate();

	nfields = 0;
	compute = modify->compute;
	parse_fields(narg, arg);

	int iarg = 6 + nfields;
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
		else if (!strcmp(arg[iarg], "flush")) {
			if (!strcmp(arg[iarg + 1], "yes")) flush_flag = 1;
			else if (!strcmp(arg[iarg + 1], "no")) flush_flag = 0;
			else error->all(FLERR, "Illegal analyze ave/time flush option");
			iarg += 2;
		}
		else if (!strcmp(arg[iarg], "start")) {
			nstart = atoi(arg[iarg + 1]);
			iarg += 2;
		}
		else error->all(FLERR, "Illegal analyz ave/space filed name or option");
	}

	if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0) {
		error->all(FLERR,"Illegal analyze ave/time command");
	}
	if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq) {
		error->all(FLERR,"Illegal analyze ave/time command");
	}

	if ((scalar_flag || vector_flag == 1) && (array_flag == 1)) {
		error->all(FLERR, "analyze ave/time field_name is not consisent with its style");
	}

	for (int i = 0; i < nfields; i++) {
		ncolumns += field_ncols[i];
	}
	
	memory->create(array, nrows, ncolumns, "AnalyzeAveTime: array");
	memory->create(array_total, nrows, ncolumns, "AnalyzeAveTime: array_total");
	memory->create(array_ave, nrows, ncolumns, "AnalyzeAveTime: array_ave");
	clean_array();		
}

/* ---------------------------------------------------------------------- */

AnalyzeAveTime::~AnalyzeAveTime()
{
	delete[] fname;
	fname = NULL;

	if (array_flag == 1) {
		memory->destroy(array);
		memory->destroy(array_total);
		memory->destroy(array_ave);
	}

	if (file && procid == 0) {
		fclose(file);
		file = NULL;
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::init()
{
	// update compute pointer
	if (modify->ncomputes > 0) compute = modify->compute;

	for (int i = 0; i < nfields; i++) {
		if (field_data_type[i] == INT) strcpy(field_format[i], "%8d ");
		else if (field_data_type[i] == FLOAT) strcpy(field_format[i], "%12.8f ");
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::allocate()
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

void AnalyzeAveTime::parse_fields(int narg, char **arg)
{
	int iarg = 6;
	while (iarg < narg) {
		if (strcmp(arg[iarg], "step") == 0) {
			addfield("Step", &AnalyzeAveTime::compute_step, SCALAR, INT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "x") == 0) {
			addfield("x", &AnalyzeAveTime::compute_x, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "y") == 0) {
			addfield("y", &AnalyzeAveTime::compute_y, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "z") == 0) {
			addfield("z", &AnalyzeAveTime::compute_z, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "vx") == 0) {
			addfield("vx", &AnalyzeAveTime::compute_vx, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "vy") == 0) {
			addfield("vy", &AnalyzeAveTime::compute_vy, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "vz") == 0) {
			addfield("vz", &AnalyzeAveTime::compute_vz, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "fx") == 0) {
			addfield("fx", &AnalyzeAveTime::compute_fx, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "fy") == 0) {
			addfield("fy", &AnalyzeAveTime::compute_fy, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "fz") == 0) {
			addfield("fz", &AnalyzeAveTime::compute_fz, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "wx") == 0) {
			addfield("wx", &AnalyzeAveTime::compute_wx, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "wy") == 0) {
			addfield("wy", &AnalyzeAveTime::compute_wy, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "wz") == 0) {
			addfield("wz", &AnalyzeAveTime::compute_wz, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "max_vx") == 0) {
			addfield("max_vx", &AnalyzeAveTime::compute_max_vx, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "max_vy") == 0) {
			addfield("max_vy", &AnalyzeAveTime::compute_max_vy, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "max_vz") == 0) {
			addfield("max_vz", &AnalyzeAveTime::compute_max_vz, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "min_vx") == 0) {
			addfield("min_vx", &AnalyzeAveTime::compute_min_vx, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "min_vy") == 0) {
			addfield("min_vy", &AnalyzeAveTime::compute_min_vy, SCALAR, FLOAT);
			iarg++;
		}
		else if (strcmp(arg[iarg], "min_vz") == 0) {
			addfield("min_vz", &AnalyzeAveTime::compute_min_vz, SCALAR, FLOAT);
			iarg++;
		}
		else if (strncmp(arg[iarg], "c_", 2) == 0) {
			int n = strlen(arg[iarg]) - 1;
			char *str = new char[n];
			int i;
			// copy compute id
			for (i = 2; i < n + 2; i++) {
				str[i - 2] = arg[iarg][i];
			}
			cid = modify->find_compute(str);
			field_index[nfields] = cid;
			if (compute[cid]->scalar_flag) {
				field_ncols[nfields] = compute[cid]->size_array_columns;
				addfield(arg[iarg], &AnalyzeAveTime::compute_compute, SCALAR, FLOAT);
			}
			else if (compute[cid]->vector_flag) {
				field_ncols[nfields] = compute[cid]->size_vector;
				addfield(arg[iarg], &AnalyzeAveTime::compute_compute, VECTOR, FLOAT);
			}
			else if (compute[cid]->array_flag) {
				field_ncols[nfields] = compute[cid]->size_array_columns;
				if (nfields == 0) nrows = compute[cid]->size_array_rows;
				if (nrows != compute[cid]->size_array_rows) error->all(FLERR, "For array type data, the number of rows have to be the same!");
				addfield(arg[iarg], &AnalyzeAveTime::compute_compute, ARRAY, FLOAT);
			}
			
			delete[] str;
			str = NULL;
			iarg++;
		}
		else break;
	} // while (iarg < nargs)
}

/* ----------------------------------------------------------------------
Add field to list of quantities to print
------------------------------------------------------------------------- */

void AnalyzeAveTime::addfield(const char *key, FnPtr func, int typeflag, int data_typeflag)
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

void AnalyzeAveTime::invoke_analyze() 
{
	int ntimestep = update->ntimestep;

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
			
			if (field_type[ifield] == SCALAR) {
				array_total[0][icol] = scalar;
				icol++;
			}
			else if (field_type[ifield] == VECTOR) {
				for (int i = 0; i < field_ncols[ifield]; i++) {
					array_total[0][icol+i] = vector[icol+i];
				}
				icol += field_ncols[ifield];
			}
			else if (field_type[ifield] == ARRAY) {
				for (int i = 0; i < nrows; i++) {
					for (int j = 0; j < field_ncols[ifield]; j++) {
						array_total[i][icol+j] = array[i][icol+j];
					}
				}
				icol += field_ncols[ifield];
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

void AnalyzeAveTime::write_array()
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
			if (flush_flag) fflush(file);
		} // if (procid == 0)
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::clean_array()
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

void AnalyzeAveTime::compute_compute()
{
	int cid; 
	cid = field_index[ifield];              // compute id

	if (field_type[ifield] == SCALAR) {
		if (compute[cid]->invoked_scalar != update->ntimestep) {
			scalar = compute[cid]->compute_scalar();
		}
		else scalar = compute[cid]->scalar;
	}
	else if (field_type[ifield] == VECTOR) {
		if (compute[cid]->invoked_vector != update->ntimestep) {
			compute[cid]->compute_vector();
		}
		for (int i = 0; i < compute[cid]->size_array_columns; i++) {
			vector[icol+i] = compute[cid]->vector[i];
		}
	}
	else if (field_type[ifield] == ARRAY) {
		if (compute[cid]->invoked_array != update->ntimestep) {
			compute[cid]->compute_array();
		}
		for (int i = 0; i < nrows; i++) {
			for (int j = 0; j < compute[cid]->size_array_columns; j++) {
				array[i][j+icol] = compute[cid]->array[i][j];
			}
		}
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_step()
{
	scalar = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_x()
{
	int *mask = particle->mask;
	double **x = particle->x;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += x[i][0];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_y()
{
	int *mask = particle->mask;
	double **x = particle->x;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += x[i][1];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_z()
{
	int *mask = particle->mask;
	double **x = particle->x;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += x[i][2];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_vx()
{
	int *mask = particle->mask;
	double **v = particle->v;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += v[i][0];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_vy()
{
	int *mask = particle->mask;
	double **v = particle->v;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += v[i][1];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_vz()
{
	int *mask = particle->mask;
	double **v = particle->v;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += v[i][2];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_fx()
{
	int *mask = particle->mask;
	double **f = particle->f;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += f[i][0];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_fy()
{
	int *mask = particle->mask;
	double **f = particle->f;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += f[i][1];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_fz()
{
	int *mask = particle->mask;
	double **f = particle->f;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += f[i][2];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_wx()
{
	int *mask = particle->mask;
	double **omega = particle->omega;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += omega[i][0];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_wy()
{
	int *mask = particle->mask;
	double **omega = particle->omega;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += omega[i][1];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_wz()
{
	int *mask = particle->mask;
	double **omega = particle->omega;
	int nlocal = particle->nlocal;

	scalar = 0.0;
	double temp = 0.0;
	int counter = 0;
	int total_counter;
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp += omega[i][2];
			counter++;
		}
	}

	MPI_Allreduce(&counter, &total_counter, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);
	if (total_counter > 0) scalar /= total_counter;
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_max_vx()
{
	int *mask = particle->mask;
	double **v = particle->v;
	int nlocal = particle->nlocal;

	double temp = -std::numeric_limits<double>::max();
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp = MAX(temp, v[i][0]);
		}
	}

	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_MAX, mworld);
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_max_vy()
{
	int *mask = particle->mask;
	double **v = particle->v;
	int nlocal = particle->nlocal;

	double temp = -std::numeric_limits<double>::max();
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp = MAX(temp, v[i][1]);
		}
	}

	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_MAX, mworld);
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_max_vz()
{
	int *mask = particle->mask;
	double **v = particle->v;
	int nlocal = particle->nlocal;

	double temp = -std::numeric_limits<double>::max();
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp = MAX(temp, v[i][2]);
		}
	}

	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_MAX, mworld);
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_min_vx()
{
	int *mask = particle->mask;
	double **v = particle->v;
	int nlocal = particle->nlocal;

	double temp = std::numeric_limits<double>::max();
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp = MIN(temp, v[i][0]);
		}
	}

	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_MIN, mworld);
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_min_vy()
{
	int *mask = particle->mask;
	double **v = particle->v;
	int nlocal = particle->nlocal;

	double temp = std::numeric_limits<double>::max();
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp = MIN(temp, v[i][1]);
		}
	}

	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_MIN, mworld);
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveTime::compute_min_vz()
{
	int *mask = particle->mask;
	double **v = particle->v;
	int nlocal = particle->nlocal;

	double temp = std::numeric_limits<double>::max();
	for (int i = 0; i < nlocal; i++) {
		if (groupbit & mask[i]) {
			temp = MIN(temp, v[i][2]);
		}
	}

	MPI_Allreduce(&temp, &scalar, 1, MPI_DOUBLE, MPI_MIN, mworld);
}

/* ---------------------------------------------------------------------- */
