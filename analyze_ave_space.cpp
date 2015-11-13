/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "analyze_ave_space.h"
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

#define DELTA 4

enum{INT, FLOAT};
enum{SCALAR, VECTOR, ARRAY};
enum{LOWER, CENTER, UPPER, COORD};
enum{NUM, VOL, AREA};

/* ---------------------------------------------------------------------- */

AnalyzeAveSpace::AnalyzeAveSpace(PDPS *ps, int narg, char **arg) : Analyze(ps, narg, arg)
{
	if (narg < 6) error->all(FLERR, "Illegal analyze ave/time command");

	array_flag = 1;

	ave_flag = 0;
	clean_ave_flag = 0;
	num_flag = 0;
	area_flag = 0;
	vol_flag = 0;

	header_flag = 1;
	header_once_flag = 0;

	maxncells = 0;

	nrows = 0;
	ncolumns = 0;
	
	field_func = NULL;

	coord_cell = NULL;
	num_cell = NULL;
	numAll_cell = NULL;
	area_cell = NULL;
	vol_cell = NULL;

	// parse average operation parameters 
	nevery = atoi(arg[3]);
	nrepeat = atoi(arg[4]);
	nfreq = atoi(arg[5]);

	// initialization of dim[3] and originflag[3]
	ndims = 0;

	int iarg = 6;
	for (int i = 0; i < 3; i++) {
		dim[i] = -1;
		originflag[i] = -1;
	}

	// parse spatial cell information
	while (iarg < narg) {
		if (iarg + 3 > narg) break;	
		if (strcmp(arg[iarg],"x") == 0) dim[ndims] = 0;
		else if (strcmp(arg[iarg],"y") == 0) dim[ndims] = 1;
		else if (strcmp(arg[iarg],"z") == 0) dim[ndims] = 2;
		else break;

		if (dim[ndims] == 2 && domain->dim == 2) {
			error->all(FLERR, "Cannot use analyze ave/spatial z for 2D model");
		}

		if (strcmp(arg[iarg + 1], "lower") == 0) originflag[ndims] = LOWER;
		else if (strcmp(arg[iarg + 1], "center") == 0) originflag[ndims] = CENTER;
		else if (strcmp(arg[iarg + 1], "upper") == 0) originflag[ndims] = UPPER;
		else originflag[ndims] = COORD;
		if (originflag[ndims] == COORD) origin[ndims] = atof(arg[iarg + 1]);

		delta[ndims] = atof(arg[iarg+2]);

		ndims++;
		iarg += 3;
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
		if ((dim[0]+dim[1]) == 1) delta[2] = boxle[2];
		else if ((dim[0]+dim[1]) == 2) delta[2] = boxle[1];
		else if ((dim[0]+dim[1]) == 3) delta[2] = boxle[0];
	}
	for (int i = 0; i < 3; i++)	{
		inv_delta[i] = 1.0 / delta[i];
	}
	
	if (ndims == 0) error->all(FLERR,"Illegal analyze ave/space command");
	if (ndims == 2 && dim[0] == dim[1]) {
		error->all(FLERR, "Same dimension twice in analyze ave/space command");
	}
	if (ndims == 3 && (dim[0] == dim[1] || dim[1] == dim[2] || dim[0] == dim[2])) {
		error->all(FLERR, "Same dimension twice in analyze ave/space");
	}

	// initial field allocate
	nfields_initial = narg - iarg;
	allocate();

	// parse fields
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
		error->all(FLERR,"Illegal analyze ave/space average operation parameters");
	}
	if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq) {
		error->all(FLERR,"Illegal analyze ave/space average operation parameters");
	}

	for (int i = 0; i < nfields; i++) {
		ncolumns += field_ncols[i];
	}

	setup_cells();
}

/* ---------------------------------------------------------------------- */

AnalyzeAveSpace::~AnalyzeAveSpace()
{
	delete[] fname;
	fname = NULL;
	
	for (int i = 0; i < nfields_initial; i++) {
		delete[] field_name[i];
		field_name[i] = NULL;
	}
	delete[] field_name;
	field_name = NULL;
	
	memory->destroy(array);
	memory->destroy(array_total);
	memory->destroy(array_ave);
	memory->destroy(coord_cell);

	if (area_flag) {
		memory->destroy(area_cell);
	}
	if (vol_flag) {
		memory->destroy(vol_cell);
	}

	if (file && procid == 0) {
		fclose(file);
		file = NULL;
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveSpace::init()
{
	// update compute pointer
	if (modify->ncomputes > 0) compute = modify->compute;

	for (int i = 0; i < nfields; i++) {
		if (field_data_type[i] == INT) strcpy(field_format[i], "%8d ");
		else if (field_data_type[i] == FLOAT) strcpy(field_format[i], "%12.8f ");
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveSpace::setup_cells() 
{
	int i, j, k, idim;
	double area, vol;
	double *boxlo, *boxhi, *boxle;

	boxlo = domain->boxlo;
	boxhi = domain->boxhi;
	boxle = domain->boxle;

	// define cell origin
	ncells = 1;
	for (i = 0; i < 3; i++) cell[i] = 1;
	for (i = 0; i < ndims; i++) {
		if (originflag[i] == LOWER) origin[i] = boxlo[dim[i]];
		else if (originflag[i] == UPPER) origin[i] = boxhi[dim[i]];
		else if (originflag[i] == CENTER) {
			origin[i] = 0.5 * (boxlo[dim[i]] + boxhi[dim[i]]);
		}

		cell[i] = static_cast<int> (boxle[dim[i]] * inv_delta[i]);
		if (cell[i] == 0) cell[i] = 1;
		ncells *= cell[i];
	}

	nrows = ncells;
	if (ncells > maxncells) {
		maxncells = ncells;
		memory->grow(array, nrows, ncolumns, "AnalyzeAveSpace: array");
		memory->grow(array_total, nrows, ncolumns, "AnalyzeAveSpace: array_total");
		memory->grow(array_ave, nrows, ncolumns, "AnalyzeAveSpace: array_ave");
		memory->grow(coord_cell, nrows, ndims, "AnalyzeAveSpace: coord_cell");
		memory->grow(num_cell, nrows, "AnalyzeAveSpace: num_cell");
		memory->grow(numAll_cell, nrows, "AnalyzeAveSpace: numAll_cell");
		if (area_flag == 1) {
			memory->grow(area_cell, nrows, "AnalyzeAveSpace: area_cell");
		}
		if (vol_flag == 1) {
			memory->grow(vol_cell, nrows, "AnalyzeAveSpace: vol_cell");
		}
	}

	// standard volume for most of cells 
	if (area_flag == 1) {
		area = delta[0] * delta[1];
	}
	else if (vol_flag == 1) {
		vol = delta[0] * delta[1] * delta[2];
	}

	// cells at the upper edge needs to be recalculated
	int cell_id;
	for (k = 0; k < cell[2]; k++)
	for (j = 0; j < cell[1]; j++)
	for (i = 0; i < cell[0]; i++) {
		cell_id = k*cell[1]*cell[0] + j*cell[0] + i;

		// find coordinate of the cell
		if (i < cell[0] - 1) {
			coord_cell[cell_id][0] = boxlo[dim[0]] + (i + 0.5)*delta[0] - origin[0];
		}
		else  {
			coord_cell[cell_id][0] = 0.5*((boxlo[dim[0]] + (cell[0] - 1)*delta[0]) \
				                     + boxhi[dim[0]]) - origin[0];
		}
		if (ndims >= 2) {
			if (j < cell[1] - 1)
				coord_cell[cell_id][1] = boxlo[dim[1]] + (j + 0.5)*delta[1] - origin[1];
			else {
				coord_cell[cell_id][1] = 0.5*((boxlo[dim[1]] + (cell[1] - 1)*delta[1]) \
					                     + boxhi[dim[1]]) - origin[1];
			}
		}
		if (ndims == 3) {
			if (k < cell[2] - 1 && ndims == 3)
				coord_cell[cell_id][2] = boxlo[dim[2]] + (k + 0.5)*delta[2] - origin[2];
			else 
				coord_cell[cell_id][2] = 0.5*((boxlo[dim[2]] + (cell[2] - 1)*delta[2]) \
										 + boxhi[dim[2]]) - origin[2];
		}

		// find volume of the cell if required
		if (area_flag == 1) {
			area_cell[cell_id] = area;
			// if the cell is at the edge along dim[0]
			if (i == cell[0] - 1) {
				area_cell[cell_id] *= (inv_delta[0] * (boxhi[dim[0]] - (cell[0] - 1)*delta[0]));
			}
			// if the cell is at the edge along dim[1]
			if (j == cell[0] - 1 && ndims == 2) {
				area_cell[cell_id] *= (inv_delta[1] * (boxhi[dim[1]] - (cell[1] - 1)*delta[1]));
			}
		}
		else if (vol_flag == 1) {
			vol_cell[cell_id] = vol;
			// if the cell is at the edge along dim[0]
			if (i == cell[0] - 1) {
				vol_cell[cell_id] *= (inv_delta[0] * (boxhi[dim[0]] - (cell[0]-1)*delta[0]));
			}
			// if the cell is at the edge along dim[1]
			if (j == cell[1] - 1 && ndims >= 2) {
				vol_cell[cell_id] *= (inv_delta[1] * (boxhi[dim[1]] - (cell[1]-1)*delta[1]));
			}
			// if the cell is at the edge along dim[2]
			if (k == cell[2] - 1 && ndims == 3) {
				vol_cell[cell_id] *= (inv_delta[2] * (boxhi[dim[2]] - (cell[2]-1)*delta[2]));
			}
		}
	} // for (k = 0; k < cell[2]; k++)
	  // for (j = 0; j < cell[1]; j++)
	  // for (i = 0; i < cell[0]; i++)
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveSpace::count_num_cell() 
{
	int nlocal = particle->nlocal;
	int *mask = particle->mask;
	int cell_id, i;

	for (i = 0; i < ncells; i++) {
		num_cell[i] = 0.0;
		numAll_cell[i] = 0.0;
	}

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			cell_id = find_cell_id(i);
			num_cell[cell_id] += 1.0;
		}
	}
	
	MPI_Allreduce(&num_cell[0], &numAll_cell[0], ncells, MPI_INT, MPI_SUM, mworld);
}

/* ---------------------------------------------------------------------- */

int AnalyzeAveSpace::find_cell_id(int pid) 
{
	int cell_id;
	int c[3];
	double **x = particle->x;
	double *boxlo = domain->boxlo;
	double *boxhi = domain->boxhi;
	
	for (int i = 0; i < 3; i++) c[i] = 0;

	int idim;
	for (int i = 0; i < ndims; i++) {
		idim = dim[i];
		// particles may not be in the box
		if (x[pid][idim] < boxlo[idim]) c[i] = 0;
		else if (x[pid][idim] >= boxhi[idim]) c[i] = cell[i] - 1;
		else c[i] = static_cast<int> ((x[pid][idim] - boxlo[idim]) / delta[i]);
		// the last cell may have large size
		if (c[i] == cell[i]) c[i] = cell[i] - 1;
	}

	cell_id = c[2]*cell[1]*cell[0] + c[1]*cell[0] + c[0];

	return cell_id;	
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveSpace::allocate()
{
	int n = nfields_initial;

	field_ave_flag = new int[n];
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

void AnalyzeAveSpace::parse_fields(int iarg, int narg, char **arg)
{
	while (iarg < narg) {
		if (!strcmp(arg[iarg], "x")) {
			num_flag = 1;
			field_ave_flag[nfields] = NUM;
			addfield("x", &AnalyzeAveSpace::compute_x, SCALAR, FLOAT);
			iarg++;
		}
		else if (!strcmp(arg[iarg], "y")) {
			num_flag = 1;
			field_ave_flag[nfields] = NUM;
			addfield("y", &AnalyzeAveSpace::compute_y, SCALAR, FLOAT);
			iarg++;
		}
		else if (!strcmp(arg[iarg], "z")) {
			num_flag = 1;
			field_ave_flag[nfields] = NUM;
			addfield("z", &AnalyzeAveSpace::compute_z, SCALAR, FLOAT);
			iarg++;
		}
		else if (!strcmp(arg[iarg], "density/number")) {
			if (domain->dim == 2) {
				area_flag = 1;
				field_ave_flag[nfields] = AREA;
			}
			else {
				vol_flag = 1;
				field_ave_flag[nfields] = VOL;
			}
			addfield("Density/Number", &AnalyzeAveSpace::compute_density_number, SCALAR, FLOAT);
			iarg++;
		}
		else if (!strcmp(arg[iarg], "density/mass")) {
			if (domain->dim == 2) {
				area_flag = 1;
				num_flag = 1;
				field_ave_flag[nfields] = AREA;
			}
			else {
				vol_flag = 1;
				num_flag = 1;
				field_ave_flag[nfields] = VOL;
			}
			addfield("Density/Mass", &AnalyzeAveSpace::compute_density_mass, SCALAR, FLOAT);
			iarg++;
		}
		else if (!strcmp(arg[iarg], "density/particle")) {
			if (particle->density_flag == 0) {
				error->all(FLERR, "Illegal particle_style type for analyze density/particle");
			}
			num_flag = 1;
			field_ave_flag[nfields] = NUM;
			addfield("Density/Particle", &AnalyzeAveSpace::compute_density_particle, SCALAR, FLOAT);
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
				addfield(arg[iarg], &AnalyzeAveSpace::compute_compute, SCALAR, FLOAT);
			}
			else if (compute[cid]->vector_flag) {
				field_ncols[nfields] = compute[cid]->size_vector;
				addfield(arg[iarg], &AnalyzeAveSpace::compute_compute, VECTOR, FLOAT);
			}
			else if (compute[cid]->array_flag) {
				field_ncols[nfields] = compute[cid]->size_array_columns;
				if (nfields == 0) ncells = compute[cid]->size_array_rows;
				if (ncells != compute[cid]->size_array_rows) error->all(FLERR, "For array type data, the number of rows have to be the same!");
				addfield(arg[iarg], &AnalyzeAveSpace::compute_compute, ARRAY, FLOAT);
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

void AnalyzeAveSpace::addfield(const char *key, FnPtr func, int typeflag, int data_typeflag)
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

void AnalyzeAveSpace::invoke_analyze() 
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
	if ((end - ntimestep) % nevery == 0 && ntimestep >= start ) {
		if (ntimestep == start) clean_ave_flag = 1;
		else clean_ave_flag = 0;
		if (ntimestep == end) ave_flag = 1; // do average at this timestep
		else ave_flag = 0;

		clean_array();

		if (domain->box_change) setup_cells();
		if (num_flag) count_num_cell();

		// compute
		icol = 0;
		for (ifield = 0; ifield < nfields; ifield++) {
			(this->*field_func[ifield])();
			icol += field_ncols[ifield];
		}

		MPI_Allreduce(&array[0][0], &array_total[0][0], nrows*ncolumns, MPI_DOUBLE, MPI_SUM, mworld);

		icol = 0;
		for (ifield = 0; ifield < nfields; ifield++) {
			if (field_ave_flag[ifield] == NUM) {
				for (int icell = 0; icell < ncells; icell++) {
					for (int j = 0; j < field_ncols[ifield]; j++) {
						if (numAll_cell[icell] != 0) {
							array_total[icell][icol+j] /= numAll_cell[icell];
						}
					}
				}
			}
			else if (field_ave_flag[icol] == VOL) {
				for (int icell = 0; icell < ncells; icell++) {
					for (int j = 0; j < field_ncols[ifield]; j++) {
						array_total[icell][icol+j] /= vol_cell[icell];
					}
				}
			}
			else if (field_ave_flag[icol] == AREA) {
				for (int icell = 0; icell < ncells; icell++) {
					for (int j = 0; j < field_ncols[ifield]; j++) {
						array_total[icell][icol] /= area_cell[icell];
					}
				}
			}
			icol += field_ncols[ifield];
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

void AnalyzeAveSpace::write_array()
{
	int i, j;

	// when output file is specified
	if (file != NULL) {
		// master processor
		if (procid == 0) {
			if (header_flag == 1) {
				fprintf(file, "TimeStep NumberOfRows\n");
				fprintf(file, "%d %d\n", update->ntimestep, nrows);
				fprintf(file, "Cell_id ");
				for (int i = 0; i < ndims; i++) {
					if (dim[i] == 0) fprintf(file, "coord_x ");
					else if (dim[i] == 1) fprintf(file, "coord_y ");
					else if (dim[i] == 2) fprintf(file, "coord_z ");
				}
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
				fprintf(file, "%d ", i);
				for (j = 0; j < ndims; j++) {
					fprintf(file, "%f ", coord_cell[i][dim[j]]);
				}
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

void AnalyzeAveSpace::clean_array()
{
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncolumns; j++) {
			array[i][j] = 0.0;
			array_total[i][j] = 0.0;
			if (clean_ave_flag == 1) {
				array_ave[i][j] = 0.0;
			}
		}
	}
}

/* ---------------------------------------------------------------------- 
   Have not decided whether needs to consider to call compute result 
   or analyze result in ave/space. Because ave/space is simply for 
   particle-wise. The data to be averaged over spatial cell should be 
   related to particle, eg. x, v, f, density/number ...
   Otherwise, for data not related to particle, it does not make any 
   sense to use perform average over space. However, in the future, 
   some compute may compute particle-wise properties.
   Note: if we don't consider compute or analyze, nfields = ncolumns,
   which is a very important assumption!
---------------------------------------------------------------------- */

void AnalyzeAveSpace::compute_compute()
{

}

/* ---------------------------------------------------------------------- */

void AnalyzeAveSpace::compute_x()
{
	double **x = particle->x;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	int cell_id;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			cell_id = find_cell_id(i);
			array[cell_id][icol] += particle->x[i][0];
		}
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveSpace::compute_y()
{
	double **x = particle->x;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	int cell_id;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			cell_id = find_cell_id(i);
			array[cell_id][icol] += particle->x[i][1];
		}
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveSpace::compute_z()
{
	double **x = particle->x;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	int cell_id;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			cell_id = find_cell_id(i);
			array[cell_id][icol] += particle->x[i][2];
		}
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveSpace::compute_density_number()
{
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	int cell_id;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			cell_id = find_cell_id(i);
			array[cell_id][icol] += 1.0;
		}
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveSpace::compute_density_mass()
{
	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	int cell_id;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			cell_id = find_cell_id(i);
			if (particle->rmass_flag) array[cell_id][icol] += particle->rmass[i];
			else array[cell_id][icol] += particle->mass[type[i]];
		}
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeAveSpace::compute_density_particle()
{
	int *mask = particle->mask;
	int *type = particle->type;
	int nlocal = particle->nlocal;
	int cell_id;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			cell_id = find_cell_id(i);
			array[cell_id][icol] += particle->density[type[i]];
		}
	}
}
