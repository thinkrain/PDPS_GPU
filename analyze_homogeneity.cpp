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

#include "analyze_homogeneity.h"
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
#include "style_analyze.h"

using namespace PDPS_NS;

enum{INT, FLOAT};
enum{SCALAR, VECTOR, ARRAY};
enum{LOCAL_AVERAGE, CENTROID};

/* ---------------------------------------------------------------------- */

AnalyzeHomogeneity::AnalyzeHomogeneity(PDPS *ps, int narg, char **arg) : Analyze(ps, narg, arg)
{
	if (narg < 7) error->all(FLERR, "Illegal analyze ave/time command");

	nstart = 0;
	clean_ave_flag = 0;

	scalar_flag = vector_flag = 0;
	array_flag = 1;

	num_min = 0;

	nrows = ncolumns = 0;

	header_flag = 1;
	header_once_flag = 0;

	eta_s_flag = 0;
	eta_s = eta_r = eta = NULL;

	field_func = NULL;

	nfreq = atoi(arg[3]);

	if (!strcmp(arg[4], "local/average")) {
		method = LOCAL_AVERAGE;
	}
	else if (!strcmp(arg[4], "centroid/diff")) {
		method = CENTROID;
	}

	// parse field_name
	nfields_initial = narg - 5;
	allocate();

	nfields = 0;
	compute = modify->compute;
	analyze = modify->analyze;

	int iarg = 5;
	parse_fields(iarg, narg, arg);
	
	iarg += nfields;
	// parse options
	while (iarg < narg) {
		if (strcmp(arg[iarg], "file") == 0) {
			int n = strlen(arg[iarg+1]) + 1;
			fname = new char[n];
			strcpy(fname, arg[iarg+1]);
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
		else if (!strcmp(arg[iarg], "num_min")) {
			num_min = atoi(arg[iarg+1]);
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
		else error->all(FLERR, "Illegal analyz ave/homogeneity option");
	}

	if ((scalar_flag || vector_flag == 1) && (array_flag == 1)) {
		error->all(FLERR, "analyze ave/time field_name is not consisent with its style");
	}

	nrows = 1;              // it is always 1 for this class
	for (int i = 0; i < nfields; i++) {
		ncolumns += field_ncols[i];
	}

	memory->create(array, nrows, ncolumns, "AnalyzeHomogeneity: array");
	memory->create(array_total, nrows, ncolumns, "AnalyzeHomogeneity: array_total");
	memory->create(array_ave, nrows, ncolumns, "AnalyzeHomogeneity: array_ave");
	clean_array();
}

/* ---------------------------------------------------------------------- */

AnalyzeHomogeneity::~AnalyzeHomogeneity()
{
	if (array_flag == 1) {
		memory->destroy(array);
		memory->destroy(array_total);
		memory->destroy(array_ave);
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeHomogeneity::init()
{
	// update compute pointer
	if (modify->ncomputes > 0) compute = modify->compute;

	for (int i = 0; i < nfields; i++) {
		if (field_data_type[i] == INT) strcpy(field_format[i], "%8d ");
		else if (field_data_type[i] == FLOAT) strcpy(field_format[i], "%12.8f ");
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeHomogeneity::allocate()
{
	int n = nfields_initial;

	if (method == LOCAL_AVERAGE) {
		eta_s = new double[n];
		eta_r = new double[n];
		eta = new double[n];
		for (int i = 0; i < n; i++) {
			eta[i] = eta_s[i] = eta_r[i] = 0.0;
		}
	}

	field_format = new char*[n];
	// store field_names
	field_name = new char*[n];
	field_data_type = new int[n];
	field_type = new int[n];
	field_func = new FnPtr[n];
	field_nrows = new int[n];
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

void AnalyzeHomogeneity::parse_fields(int iarg, int narg, char **arg)
{
	while (iarg < narg) {
		if (strncmp(arg[iarg], "a_", 2) == 0) {
			int n = strlen(arg[iarg]) - 1;
			char *str = new char[n];
			int i;
			// copy analyze id
			for (i = 2; i < n + 2; i++) {
				str[i - 2] = arg[iarg][i];
			}
			aid = modify->find_analyze(str);
			if (aid == -1) {
				error->all(FLERR, "Cannot find the analyze id");
			}
			field_index[nfields] = aid;
			if (analyze[aid]->scalar_flag) {
				field_nrows[nfields] = analyze[aid]->nrows;
				field_ncols[nfields] = analyze[aid]->ncolumns;
				addfield(arg[iarg], &AnalyzeHomogeneity::compute_analyze, SCALAR, FLOAT);
			}
			else if (analyze[aid]->vector_flag) {
				field_nrows[nfields] = analyze[aid]->nrows;
				field_ncols[nfields] = analyze[aid]->ncolumns;
				addfield(arg[iarg], &AnalyzeHomogeneity::compute_analyze, VECTOR, FLOAT);
			}
			else if (analyze[aid]->array_flag) {
				field_nrows[nfields] = analyze[aid]->nrows;
				field_ncols[nfields] = analyze[aid]->ncolumns;
				addfield(arg[iarg], &AnalyzeHomogeneity::compute_analyze, ARRAY, FLOAT);
			}

			delete[] str;
			str = NULL;
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
			if (cid == -1) {
				error->all(FLERR, "Cannot find the compute id");
			}
			field_index[nfields] = cid;
			if (compute[cid]->scalar_flag) {
				field_ncols[nfields] = compute[cid]->size_array_columns;
				addfield(arg[iarg], &AnalyzeHomogeneity::compute_compute, SCALAR, FLOAT);
			}
			else if (compute[cid]->vector_flag) {
				field_ncols[nfields] = compute[cid]->size_vector;
				addfield(arg[iarg], &AnalyzeHomogeneity::compute_compute, VECTOR, FLOAT);
			}
			else if (compute[cid]->array_flag) {
				field_ncols[nfields] = compute[cid]->size_array_columns;
				if (nfields == 0) nrows = compute[cid]->size_array_rows;
				if (nrows != compute[cid]->size_array_rows) error->all(FLERR, "For array type data, the number of rows have to be the same!");
				addfield(arg[iarg], &AnalyzeHomogeneity::compute_compute, ARRAY, FLOAT);
			}

			delete[] str;
			str = NULL;
			iarg++;
		}
		else break;
		if (method == CENTROID) field_ncols[nfields - 1]++;
	} // while (iarg < nargs)
}

/* ----------------------------------------------------------------------
Add field to list of quantities to print
------------------------------------------------------------------------- */

void AnalyzeHomogeneity::addfield(const char *key, FnPtr func, int typeflag, int data_typeflag)
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

void AnalyzeHomogeneity::invoke_analyze()
{
	int ntimestep = update->ntimestep;

	// example: nevery = 10 nrepeat = 5 nfreq = 1000
	// time averages are calculated on 960 970 980 990 1000 steps
	if (ntimestep % nfreq == 0 && ntimestep >= nstart) {
		clean_array();
		
		// compute
		icol = 0;
		for (ifield = 0; ifield < nfields; ifield++) {
			(this->*field_func[ifield])();
			icol += field_ncols[ifield];
		}

		icol = 0;
		if (method == LOCAL_AVERAGE) {
			for (ifield = 0; ifield < nfields; ifield++) {
				for (int j = 0; j < field_ncols[ifield]; j++) {
					array_ave[0][icol+j] = (eta[icol+j] - eta_r[icol+j]) / (eta_s[icol+j] - eta_r[icol+j]);
				}
				icol += field_ncols[ifield];
			}
		}

		write_array();
		
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeHomogeneity::write_array()
{
	int i, j;

	// when output file is specified
	if (file != NULL) {
		// master processor
		if (procid == 0) {
			if (header_flag == 1) {
				fprintf(file, "TimeStep NumberOfRows\n");
				fprintf(file, "%d %d\n", update->ntimestep, nrows);
				fprintf(file, "Step ");
				for (ifield = 0; ifield < nfields; ifield++) {
					if (field_type[ifield] == SCALAR) fprintf(file, "%s ", field_name[ifield]);
					else if (field_type[ifield] == VECTOR || field_type[ifield] == ARRAY){
						for (j = 0; j < field_ncols[ifield]; j++) {
							fprintf(file, "%s[%d] ", field_name[ifield], j);
						}
					}
				}
				fprintf(file, "\n");
				if (header_once_flag == 1) header_flag = 0;
			}
			int temp = field_ncols[0];
			int curr_field = 0;
			int ivalue;
			double dvalue;
			for (i = 0; i < nrows; i++) {
				fprintf(file, "%d ", update->ntimestep);
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
				fflush(file);
			}
		} // if (procid == 0)
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeHomogeneity::clean_array()
{
	for (int i = 0; i < nrows; i++)
	for (int j = 0; j < ncolumns; j++) {
		array[i][j] = 0.0;
		array_total[i][j] = 0.0;
		array_ave[i][j] = 0.0;
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeHomogeneity::compute_analyze()
{
	int aid;
	aid = field_index[ifield];
	analyze = modify->analyze;
	if (analyze[aid]->scalar_flag) {
		 
	}
	else if (analyze[aid]->vector_flag) {
		 
	}
	else if (analyze[aid]->array_flag) {
		if (method == LOCAL_AVERAGE) {
			compute_local_average(analyze[aid]->array_ave, analyze[aid]->numAll_cell, analyze[aid]->nrows, analyze[aid]->ncolumns);
		}
		else if (method == CENTROID) {
			compute_centroid(analyze[aid]->array_ave, analyze[aid]->numAll_cell, analyze[aid]->nrows, analyze[aid]->ncolumns);
		}
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeHomogeneity::compute_compute()
{
	int cid;
	cid = field_index[ifield];              // compute id

}

/* ---------------------------------------------------------------------- */

void AnalyzeHomogeneity::compute_local_average(double **array_attr, int *numAll_attr, int rows, int columns)
{
	double mean, sigma, variance;
	int nvalids = 0;
	for (int j = 0; j < columns; j++) {
		mean = sigma = variance = 0.0;
		nvalids = 0;
		for (int i = 0; i < rows; i++) {
			if (numAll_attr[i] >= num_min) {
				mean += array_attr[i][j];
				nvalids += 1;
			}
		}
		mean /= nvalids;
		for (int i = 0; i < rows; i++) {
			if (numAll_attr[i] >= num_min) {
				variance += (array_attr[i][j] - mean)*(array_attr[i][j] - mean);
			}
		}
		variance /= nvalids;
		sigma = sqrt(variance);
		eta[icol+j] = sigma / mean;
		if (eta_s_flag == 0) {
			eta_s[icol + j] = eta[icol + j];
		}
	}
	eta_s_flag = 1;
}

/* ---------------------------------------------------------------------- */

void AnalyzeHomogeneity::compute_centroid(double **array_attr, int *numAll_attr, int rows, int columns)
{
	int nvalids, nvalids_total;
	int *mask = particle->mask;
	double **x = particle->x;
	int nlocal = particle->nlocal;
	// calculate global centroid first
	double centroid[3], centroid_total[3];

	if (columns != 3) error->all(FLERR, "Illegal field for centroid/diff method");

	for (int i = 0; i < 3; i++) {
		centroid[i] = 0.0;
		centroid_total[i] = 0.0;
	}

	nvalids = nvalids_total = 0;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			for (int j = 0; j < 3; j++) {
				centroid[j] += x[i][j];
			}
			nvalids++;
		}
	}
	MPI_Allreduce(&nvalids, &nvalids_total, 1, MPI_INT, MPI_SUM, mworld);
	MPI_Allreduce(&centroid[0], &centroid_total[0], 3, MPI_DOUBLE, MPI_SUM, mworld);

	if (nvalids_total > 0) {
		for (int i = 0; i < 3; i++) {
			centroid_total[i] /= nvalids_total;
		}
	}

	double dist;
	nvalids = 0;
	for (int i = 0; i < rows; i++) {
		dist = 0.0;
		if (numAll_attr[i] >= num_min) {
			for (int j = 0; j < columns; j++) {
				array_ave[0][icol + j] += array_attr[i][j] - centroid_total[j];
				dist += (array_attr[i][j] - centroid_total[j]) * (array_attr[i][j] - centroid_total[j]);
			}
			dist = sqrt(dist);
			array_ave[0][icol + columns] += dist;
			nvalids++;
		}
	}
	
	if (nvalids > 0) {
		for (int j = 0; j < ncolumns; j++) {
			array_ave[0][icol + j] /= nvalids;
		}
	}
}
