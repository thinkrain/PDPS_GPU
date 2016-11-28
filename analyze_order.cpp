/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */
//#include <omp.h>
#include <limits>
#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "analyze_order.h"
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

//#define NUM_THREADS 12
#define PI 3.1416

using namespace PDPS_NS;

enum{X, V, MAX_V, F, STEP, COMPUTE};
enum{INT, FLOAT};
enum{SCALAR, VECTOR, ARRAY};
enum{ORIENTATIONAL, TRANSLATIONAL};
#define EPSILON 1.0e-10
/* ---------------------------------------------------------------------- */

AnalyzeOrder::AnalyzeOrder(PDPS *ps, int narg, char **arg) : Analyze(ps, narg, arg)
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

	num_cell = NULL;
	numAll_cell = NULL;

	r_cut = atof(arg[4]);
	nevery = atoi(arg[5]);
	nrepeat = atoi(arg[6]);
	nfreq = atoi(arg[7]);
	rid = domain->find_region(arg[8]);
	if (rid == -1) error->all(FLERR, "Cannot find the region id");

	nfields_initial = narg - 6;
	allocate();

	nfields = 0;
	compute = modify->compute;

	int iarg = 3;
	parse_fields(iarg, narg, arg);

	iarg = 9;
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

	memory->create(array, nrows, ncolumns, "AnalyzeOrder: array");
	memory->create(array_total, nrows, ncolumns, "AnalyzeOrder: array_total");
	memory->create(array_ave, nrows, ncolumns, "AnalyzeOrder: array_ave");
	clean_array();
}

/* ---------------------------------------------------------------------- */

AnalyzeOrder::~AnalyzeOrder()
{

	memory->destroy(numAll_cell);
	
	if (array_flag == 1) {
		memory->destroy(array);
		memory->destroy(array_total);
		memory->destroy(array_ave);
	}
}

/* ---------------------------------------------------------------------- */

void AnalyzeOrder::addfield(const char *key, FnPtr func, int typeflag, int data_typeflag)
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

void AnalyzeOrder::allocate()
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

void AnalyzeOrder::parse_fields(int iarg, int narg, char **arg)
{
		while (iarg < narg) {
		if (strcmp(arg[iarg], "orientational") == 0) {
			method = ORIENTATIONAL;
			field_ncols[nfields] = 1;
			addfield("orientational", &AnalyzeOrder::orientationalorder, SCALAR, FLOAT);
			iarg++;
		}
		if (strcmp(arg[iarg], "translational") == 0) {
			method = TRANSLATIONAL;
			field_ncols[nfields] = 1; 
			addfield("translational", &AnalyzeOrder::translationalorder, SCALAR, FLOAT);
			iarg++;
		}
		else break;
	} // while (iarg < nargs)

}

/* ---------------------------------------------------------------------- */

void AnalyzeOrder::init()
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
		memory->grow(tag2tracerID, max_ntags, "AnalyzeOrder: tag2tracerID");
	}*/
}

/* ---------------------------------------------------------------------- */


void AnalyzeOrder::invoke_analyze() 
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

		// sum all for average

 		for (int i = 0; i < nrows; i++){
		for (int j = 0; j < ncolumns; j++) 
			array_ave[i][j] = order;
		}

		// average
//		if (ave_flag == 1) {
//			for (int i = 0; i < nrows; i++) 
//			for (int j = 0; j < ncolumns; j++) {
//				array_ave[i][j] = array_ave[i][j] / nrepeat;
//			}
			write_array();

	}	
}


/* ---------------------------------------------------------------------- */

void AnalyzeOrder::write_array()
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

void AnalyzeOrder::clean_array()
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
/*				compute spherial harmonics                                */
/* ---------------------------------------------------------------------- */

void AnalyzeOrder::sphericalharmonics(double *Y, int m, int l, double thita, double fi)
{
	int i;
	double N;
	N=(2*l+1)/(4*PI);
	i=l-m;
	while(i>1)
	{	
		N=N*i;
		i--;
	}

	i=l+m;
	while(i>1)
	{
		N=N/i;
		i--;
	}
	N=sqrt(N);
	Y[0]=N*legendre(m,l,thita)*cos(m*fi);
	Y[1]=N*legendre(m,l,thita)*sin(m*fi);

}

/* ---------------------------------------------------------------------- */
/*				computing legendre polynomial
/* ---------------------------------------------------------------------- */
double AnalyzeOrder::legendre(int m,int l,double x)
{
	int i;
	double P;
	if (m<0)
	{
		P=legendre(-m,l,x);
		i=l+m;
		while(i>1)
		{
			P=P*i;
			i--;
		}
		i=l-m;
		while(i>1)
		{
			P=P/i;
			i--;
		}
		P=P*pow(-1.0,m);
		return P;
	}
	else if (l==0&&m==0)
		return 1;
	else if (l==1&&m==0)
		return x;
	else if (m>0)
	{
		if(1-x<0.0001)
			P=0;
		else
		{
			P=(l-m+1)*x*legendre(m-1,l,x)-(l+m-1)*legendre(m-1,l-1,x);
			P=P/sqrt(1-x*x);
		}
		return P;
	}
	else 
	{
		P=(2*l-1)*x*legendre(m,l-1,x)-(l+m-1)*legendre(m,l-2,x);
		P=P/(l-m);
		return P;
	}
}

/* ---------------------------------------------------------------------- */
/*						compute the orientationalorder					  */
/* ---------------------------------------------------------------------- */
void AnalyzeOrder::orientationalorder()
{
	int nparticles = particle->nparticles;

	int nlocal = particle->nlocal;
	double **x = particle->x;
	int *type = particle->type;
	int *mask = particle->mask;

	
	double rij, rijsq, r_cutsq;
//	double sum[NUM_THREADS], sum_none[NUM_THREADS];
	double sum;
	int sum_type, sum_none;					
	
	int *ilist;

	int **firstneigh;
	int *numneigh;


	ilist = neighbor->neighlist->ilist;
	firstneigh = neighbor->neighlist->firstneigh;
	numneigh = neighbor->neighlist->numneigh;
	r_cutsq = 2.5 * r_cut * r_cut;

//#pragma omp parallel
//	{

	int i, j, ii, jj, jnum, itype, jtype, id;
	int num_neigh;
	double xtmp, ytmp, ztmp, delx, dely, delz;
	double *thita, *fi;								// relative angle of neighbor particles
	double ***Y;									// Spherical Harnomics of each particle
//	double *order_local;							// orientational order metric of every point
//	order_local = (double *) memory->smalloc(nlocal * sizeof(double), "Local Orientational Order");
	thita = (double *)memory->smalloc(sizeof(double), "thita");
	fi = (double *)memory->smalloc(sizeof(double), "thita");
	int *jlist;
// Spherical Harnomics for each point 
	
	Y = (double ***) memory->smalloc(100 * sizeof(double **), "Local Y");
	for (jj = 0; jj < 100; jj++){
		Y[jj] = (double **) memory->smalloc(13 * sizeof(double *), "Local Y");
		for (int m = 0; m < 13; m++){
			Y[jj][m] = (double *) memory->smalloc(2 * sizeof(double), "Local Y");
			Y[jj][m][0] = 0;
			Y[jj][m][1] = 0;
		}
	}
	double Yr[13] = { 0.0 };
	double Yi[13] = { 0.0 };
	int num_bond = 0;
//	id = omp_get_thread_num();
	// loop all the particles and compute orientational order from the neighbor particles
//	for (ii = id, sum[id] = 0.0, sum_none[id] = 0.0 ; ii < inum; ii = ii + NUM_THREADS) {

	for (ii = 0; ii < nlocal; ii++ ){
		i = ilist[ii];
		if (!(mask[i] & groupbit))
			continue;
		itype = type[i];
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];
		num_neigh = 0;
		jlist = firstneigh[i];
		jnum = numneigh[i];
		for (jj = 0; jj < jnum; jj++) {
			j = jlist[jj];
			jtype = type[j];
			if (!(mask[j] & groupbit))
				continue;
			jtype = type[j];
			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			rijsq = (delx*delx + dely*dely + delz*delz);
			if (rijsq < r_cutsq){
				if (rijsq < EPSILON)
					continue;
				anglecalculate(i, j, thita, fi);
				for (int m = 0; m < 13; m++)
					sphericalharmonics(Y[num_neigh][m], m - 6, 6, *thita, *fi);
				num_neigh++;
			}
		}
//		order_local[i] = 0;

		// if no bond within cutoff distance, exclude the point from the order calculation
		if (num_neigh != 0){
		//	sum_none[id] = sum_none[id] + 1;
			// every local order is the average of the Y or each bond
			for (int m = 0; m < 13; m++){

				for (jj = 0; jj < num_neigh; jj++){
					Yr[m] += Y[jj][m][0];
					Yi[m] += Y[jj][m][1];
				}
				//Yr = Yr / num_neigh;
				//Yi = Yi / num_neigh;
				//order_local[i] = order_local[i] + Yr * Yr + Yi * Yi;
			}
			num_bond += num_neigh;
			//order_local[i] = sqrt(4 * PI / 13 * order_local[i]);
			//	sum[id] = sum[id] + order_local[i];
			//sum = sum + order_local[i];
	
	}
	}

	// order is the average of all the local order
	
//	memory->destroy(order_local);
	memory->destroy(Y);
	memory->destroy(thita);
	memory->destroy(fi);
//	}
	order = 0;

//	for (int i = 0; i < NUM_THREADS; i++){
//		order += sum[i]; 
//		num_none += sum_none[i];
//	}
	double Yrall[13], Yiall[13];
	int num_bondall;
	for (int m = 0; m < 13; m++){
		MPI_Allreduce(&Yr[m], &Yrall[m], 1, MPI_DOUBLE, MPI_SUM, mworld);
		MPI_Allreduce(&Yi[m], &Yiall[m], 1, MPI_DOUBLE, MPI_SUM, mworld);
	}
	MPI_Allreduce(&num_bond, &num_bondall, 1, MPI_INT, MPI_SUM, mworld);
	//MPI_Allreduce(&sum, &order, 1, MPI_DOUBLE, MPI_SUM, mworld);
	//MPI_Allreduce(&numlocal, &numall, 1, MPI_INT, MPI_SUM, mworld);
	for (int m = 0; m < 13; m++){
		Yrall[m] /= num_bondall;
		Yiall[m] /= num_bondall;
		order += Yrall[m] * Yrall[m] + Yiall[m] * Yiall[m];
	}
	order = sqrt(4 * PI / 13 * order);
	//order /= numall;

}


/* ---------------------------------------------------------------------- */
/*						compute the translationalorder					  */
/* ---------------------------------------------------------------------- */
void AnalyzeOrder::translationalorder()
{
	int nparticles = particle->nparticles;
	int i, j, ii, jj, inum, jnum, itype, jtype;
	int num_neigh;
	int nlocal = particle->nlocal;
	double **x = particle->x;
	int *type = particle->type;
	double xtmp, ytmp, ztmp, delx, dely, delz;
	double *order_local;
	double *thita, *fi;
	double ***Y;
	double rij, rijsq, rshell, r_cut, r_cutsq;
	int num_shell, i_shell;
	double shell[20], std_shell[20];
	double density, volume;
	order_local = (double *) memory->smalloc(nparticles * sizeof(double), "Local Orientational Order");
	volume = domain->regions[rid]->volume;


	density = nlocal / volume;
	rshell = 0.1;								// width of every shell used to calculate translational order
	num_shell = 20;								// number of shells
	r_cut = rshell * num_shell;
	r_cutsq = r_cut * r_cut;

	// std_shell is the shell distribution of random gas, i.e particles have equal probability to appear anywhere
	for (i = 0; i < 20; i++)
		std_shell[i] = 4 * PI / 3 * (pow((i + 1) * rshell, 3) - pow(i * rshell, 3));       
	order = 0;
	for (i = 0; i < nlocal; i++){
		itype = type[i];
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];
		num_neigh = 0;
		order_local[i] = 0;
		for (ii = 0; ii < 20; ii++)
			shell[ii] = 0;
		for (j = 0; j < nlocal; j++){
				jtype = type[j];
			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			rijsq = (delx*delx + dely*dely + delz*delz);
			rij = sqrt(rijsq);
			if (rijsq < r_cutsq){
				if (rijsq < EPSILON)
					continue;
				i_shell = static_cast<int> (rij / rshell);
				shell[i_shell] += 1;
			}
		}

		// this is the formula computing translational order 
		for (ii = 10; ii < 20; ii++)
			order_local[i] = order_local[i] + fabs(shell[ii] / std_shell[ii] - 1);
		order_local[i] = order_local[i] / 10;
		order += order_local[i];
	}
	order /= nlocal;
	memory->destroy(order_local);

}


/* ---------------------------------------------------------------------- */
/*					compute the angel between two points     			  */
/* ---------------------------------------------------------------------- */
void AnalyzeOrder::anglecalculate(int i, int j, double *thita, double *fi)
{
	double **x = particle->x;
	double rij, xi, xj, xij;
	xi = x[j][0] - x[i][0];
	xj = x[j][1] - x[i][1];
	rij = sqrt((x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) + (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]) + (x[i][2] - x[j][2]) * (x[i][2] - x[j][2]));
	xij = sqrt((x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) + (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]));
	thita[0] = fabs(x[i][2] - x[j][2]) / rij;
	if (xij < 0.000001)
		fi[0] = 0;
	else{
		if (xi >= 0){
			if (xj >= 0)
				fi[0] = asin(xj / xij);
			else
				fi[0] = asin(xj / xij) + 2 * PI;
		}		
		else    //  xi < 0 
			fi[0] = PI - asin(xj / xij);		
	}



}
