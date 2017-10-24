/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "domain.h"
#include "error.h"
#include "lattice.h"
#include "memory.h"
#include "parallel.h"
#include "particle.h"
#include "region.h"
#include "string.h"
#include "style_region.h"

#include "pdps_cuda.h"
#include "cuda_engine.h"
#include "device_launch_parameters.h"
#include "device_functions.h"

using namespace PDPS_NS;

#define DELTA 4         // For reallocation of memory, if the size of the array need to be expanded 

/* ----------------------------------------------------------------------
   Domain Class: Set dimension of simulation
			     Set initial simulation box
			     Set boundary conditions
			     Set periodic boundary algorithm
------------------------------------------------------------------------- */

Domain::Domain(PDPS *ps) : Pointers(ps)
{
	dim = 3;                     // default is 3-D
	maxarg = 0;
	for(int i = 0; i < 3; i++) {
		periodicity[i] = 1;      // default is periodic
		lsp[i] = 0;              // defualt is 0
		boundary[i][0] = -1;     
		boundary[i][1] = -1; 
		boxlo[i] = 0.0;  
		boxhi[i] = 0.0;
		boxle[i] = 0.0;
	}

    nregions = 0;
	regions = NULL;

	lattice = NULL;
	lsp_style = NULL;

	box_exist = 0;
	box_change = 0;

	nonperiodic = 0;
	xperiodic = 1;
	yperiodic = 1;
	zperiodic = 1;

	// check for fix deform
	deform_flag = deform_vremap = deform_groupbit = 0;

	h_rate[0] = h_rate[1] = h_rate[2] =
      h_rate[3] = h_rate[4] = h_rate[5] = 0.0;
	h_ratelo[0] = h_ratelo[1] = h_ratelo[2] = 0.0;
}

/* ---------------------------------------------------------------------- */

Domain::~Domain()
{
	for (int i = 0; i < nregions; i++) {
		delete regions[i];
		regions[i] = NULL;
	}
	memory->sfree(regions);
}

/* ----------------------------------------------------------------------
   Intialize
------------------------------------------------------------------------- */

void Domain::init()
{
	box_change = 0;
	for (int i = 0; i < nregions; i++) regions[i]->init();
}

/* ----------------------------------------------------------------------
   Define Regions
   ---------------------------------------------------------------------- */

void Domain::add_region(int narg, char **arg)
{
	if (narg < 3) { 
		error->all(FLERR,"Illegal region command");
	}

	// Reallocate memeory if index exceeds the maximum arguments
	if (nregions == maxarg) {
		maxarg += DELTA;
		regions = (Region **) \
			memory->srealloc(regions,maxarg*sizeof(Region *),"Domain: regions");
	}

	// create the region with a specific style
	if (0) return;

#define REGION_CLASS
#define RegionStyle(key,Class) \
	else if (strcmp(arg[1],#key) == 0) \
		regions[nregions] = new Class(ps,narg,arg);
#include "style_region.h"
#undef REGION_CLASS

	else error->all(FLERR, "Illegal region style");

	nregions++;
}

/* ----------------------------------------------------------------------
   Enforce PBC and modify particles's coordinate
------------------------------------------------------------------------- */

void Domain::pbc()
{
	int i, j;
	int nlocal = particle->nlocal;
	double **x = particle->x;
	double **v = particle->v;
	int *mask = particle->mask;

	for (i = 0; i < nlocal; i++) {
		for (j = 0; j < 3; j++) {
			if (periodicity[j] == 1) {
				while (x[i][j] < boxlo[j]) {
					x[i][j] += boxle[j];
				}
				while (x[i][j] > boxhi[j]) {
					x[i][j] -= boxle[j];
				}
			}
		}
	}
}

/* ----------------------------------------------------------------------
   Reset global & local boxes if box dimension is changed
------------------------------------------------------------------------- */

void Domain::reset_box()
{
	// perform shrink-wrapping
	// compute extent of particles on this proc
	// for triclinic, this is done in lamda space

	if (nonperiodic == 2) {
		int nlocal = particle->nlocal;
		double **x = particle->x;
		double lo[3], hi[3];
		double lo_all[3], hi_all[3];

		for (int i = 0; i < 3; i++) {
			lo[i] = boxlo[i];
			hi[i] = boxhi[i];
		}	
		// scan particle's outmost boundary
		for (int i = 0; i < nlocal; i++) {
			for (int j = 0; j < dim; j++) {
				lo[j] = MIN(x[i][j], lo[j]);
				hi[j] = MAX(x[i][j], hi[j]);
				
			}
		}

		MPI_Allreduce(lo,lo_all,3,MPI_DOUBLE,MPI_MIN,mworld);
		MPI_Allreduce(hi,hi_all,3,MPI_DOUBLE,MPI_MAX,mworld);

		// reset box domain
		for (int i = 0; i < 3; i++) {
			int flag = 0;
			if (periodicity[i] == 1) {
				if (lo_all[i] < boxlo[i]) {
					boxlo[i] = lo_all[i];
					flag = 1;
				}
				if (hi_all[i] > boxhi[i]) {
					boxhi[i] = hi_all[i];
					flag = 1;
				}
				if (flag == 1) boxle[i] = boxhi[i] - boxlo[i];
			}
		}
	} // if (nonperiodic == 2)

	set_global_box();
	set_local_box();
}

/* ----------------------------------------------------------------------
   Set boundary conditions
   ---------------------------------------------------------------------- */

void Domain::set_boundary(int narg, char** arg)
{
	if(narg != 3) {
		error->all(FLERR,"Illegal boundary command");
	}
	
	for (int i = 0; i < narg; i++) {
		int n = strlen(arg[i]);
		if (n < 1 || n > 2) 
			error->all(FLERR,"Illegal boundary command");
		if (n == 1) {
			if (arg[i][0] == 'p') {
				boundary[i][0] = 0;
				boundary[i][1] = 0;

			}
			else if (arg[i][0] == 'f') {
				boundary[i][0] = 1;
				boundary[i][1] = 1;
			}
			else if (arg[i][0] == 's') {
				boundary[i][0] = 2;
				boundary[i][1] = 2;
			}
		} // if (n == 1)
		else {
			for (int j = 0; j < n; j++) { 
				if (arg[i][j] == 'p') {
					boundary[i][j] = 0;
				}
				if (arg[i][j] == 'f') {
					boundary[i][j] = 1;
				}
				else if (arg[i][j] == 's') {
					boundary[i][j] = 2;
				}
			}
		} // if (n == 2)
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			if (boundary[i][j] == -1)
				error->all(FLERR,"Boundary conditions are not set completely");
		}
		if ((boundary[i][0] == 0 && boundary[i][1] !=0) ||
			(boundary[i][0] != 0 && boundary[i][1] ==0))
			error->all(FLERR,"Both side of boundary must be periodic");
	}

	// assign value to periodicity[], xperiodic, yperiodic and zperiodic 
	if (boundary[0][0] == 0 && boundary[0][1] == 0) {
		periodicity[0] = 1;
		xperiodic = 1;
	}
	else {
		xperiodic = 0;
		periodicity[0] = 0;
	}

	if (boundary[1][0] == 0 && boundary[1][1] == 0) {
		periodicity[1] = 1;
		yperiodic = 1;
	}
	else {
		yperiodic = 0;
		periodicity[1] = 0;
	}

	if (boundary[2][0] == 0 && boundary[2][1] == 0) {
		periodicity[2] = 1;
		zperiodic = 1;
	}
	else {
		zperiodic = 0;
		periodicity[2] = 0;
	}

	// assign value to nonperodic 
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 2; j++) {
		if (boundary[i][j] == 1)
			nonperiodic = 1;
		if (boundary[i][j] == 2) {
			nonperiodic = 2;
			break;
		}
	}
	
	if (domain->dim == 2 && (boundary[2][0] != 0 || boundary[2][1] != 0)) {
		error->all(FLERR,"z direction has to be periodic boundary condition for 2d simulation");
	}
}

/* ----------------------------------------------------------------------
   Set initial box (Called by "create_box" command)
   ---------------------------------------------------------------------- */

void Domain::set_initial_box() 
{
	int rid;

	
}

/* ----------------------------------------------------------------------
   Set global box (Called by "create_box" command)
   ---------------------------------------------------------------------- */

void Domain::set_global_box() 
{
	boxle[0] = boxhi[0] - boxlo[0];
	xle = boxle[0];
	boxle[1] = boxhi[1] - boxlo[1];
	yle = boxle[1];
	boxle[2] = boxhi[2] - boxlo[2];
	zle = boxle[2];
}

/* ----------------------------------------------------------------------
   Set local box (corresponding to each processor)
   ---------------------------------------------------------------------- */

void Domain::set_local_box() 
{
	int *procloc = parallel->procloc;
	int *procgrid = parallel->procgrid;
	double *xsplit = parallel->xsplit;
	double *ysplit = parallel->ysplit;
	double *zsplit = parallel->zsplit;

	sublo[0] = boxlo[0] + xle*xsplit[procloc[0]];
    if (procloc[0] < procgrid[0]-1)
      subhi[0] = boxlo[0] + xle*xsplit[procloc[0]+1];
    else subhi[0] = boxhi[0];

    sublo[1] = boxlo[1] + yle*ysplit[procloc[1]];
    if (procloc[1] < procgrid[1]-1)
      subhi[1] = boxlo[1] + yle*ysplit[procloc[1]+1];
    else subhi[1] = boxhi[1];

    sublo[2] = boxlo[2] + zle*zsplit[procloc[2]];
    if (procloc[2] < procgrid[2]-1)
      subhi[2] = boxlo[2] + zle*zsplit[procloc[2]+1];
    else subhi[2] = boxhi[2];

	for (int i = 0; i < 3; i++) 
		suble[i] = subhi[i] - sublo[i];
}

/* ----------------------------------------------------------------------
   Set lattice spacing information
   ---------------------------------------------------------------------- */

void Domain::set_lattice(int narg, char** arg)
{
	if (lattice) delete lattice;
	lattice = new Lattice(ps, narg, arg);
	
	/*
	int n = strlen(arg[0]) + 1;
	lsp_style = new char[n];
	strcpy(lsp_style,arg[0]);

	if(strcmp(arg[0],"spacing") == 0) {
		if (narg != 7) 
			error->all(FLERR,"Illegal lattice spacing arguments");
		for(int i = 0; i < 3; i++) {
			lsp[i] = atof(arg[2*(i+1)]);                // lattice spacing x 0.4 y 0.4 z 0.4
		}
	}
	else if (strcmp(arg[0],"number_density") == 0) {
		if (narg != 2)
			error->all(FLERR,"Illegal lattice number_density arguments");
		number_density = atof(arg[1]);
	}
	else {
		error->all(FLERR,"Illegal lattice command");
	}
	*/
}

/* ----------------------------------------------------------------------
   Find region index: no match, return -1
   ---------------------------------------------------------------------- */

int Domain::find_region(const char* rname)
{
	int rid;

	rid = -1;
	for (int i = 0; i < nregions; i++) {
		if (!strcmp(rname, regions[i]->name)) rid = i;       // scan for the region's name
	}
	return rid;
}

/* ----------------------------------------------------------------------
   format boundary string for output
   assume str is 9 chars or more in length
------------------------------------------------------------------------- */

void Domain::boundary_string(char *str)
{
	int m = 0;
	for (int idim = 0; idim < 3; idim++) {
		for (int iside = 0; iside < 2; iside++) {
			if (boundary[idim][iside] == 0) str[m++] = 'p';
			else if (boundary[idim][iside] == 1) str[m++] = 'f';
			else if (boundary[idim][iside] == 2) str[m++] = 's';
			else if (boundary[idim][iside] == 3) str[m++] = 'm';
		}
		str[m++] = ' ';
	}
	str[8] = '\0';
}

/* ----------------------------------------------------------------------
   Printf simulation box information
   ---------------------------------------------------------------------- */

void Domain::print_box(const char *str) 
{
	if (parallel->procid == 0) {
		if (screen) {
			fprintf(screen,"Simulation box is %s successfully:\n",str);
			fprintf(screen,"x: %g to %g\ny: %g to %g\nz: %g to %g\n",boxlo[0],boxhi[0],boxlo[1],boxhi[1],boxlo[2],boxhi[2]);
		}
		if (logfile) {
			fprintf(logfile,"Simulation box is %s successfully:\n",str);
			fprintf(logfile,"x: %g to %g\ny: %g to %g\nz: %g to %g\n",boxlo[0],boxhi[0],boxlo[1],boxhi[1],boxlo[2],boxhi[2]);
		}
	}

}

/* ----------------------------------------------------------------------
   GPU memcpy
   ---------------------------------------------------------------------- */

void Domain::gpuMemForRegions(void) 
{
	cudaMalloc(&devStyle, nregions * sizeof(int));
	cudaMalloc(&devRadiusCylinder, nregions * sizeof(double));
	cudaMalloc(&devHeight, nregions * sizeof(double));
	cudaMalloc(&devCoord1X, nregions * sizeof(double));
	cudaMalloc(&devCoord1Y, nregions * sizeof(double));
	cudaMalloc(&devCoord1Z, nregions * sizeof(double));
	cudaMalloc(&devCoord2X, nregions * sizeof(double));
	cudaMalloc(&devCoord2Y, nregions * sizeof(double));
	cudaMalloc(&devCoord2Z, nregions * sizeof(double));
	cudaMalloc(&devCoord3X, nregions * sizeof(double));
	cudaMalloc(&devCoord3Y, nregions * sizeof(double));
	cudaMalloc(&devCoord3Z, nregions * sizeof(double));
	cudaMalloc(&devCoord4X, nregions * sizeof(double));
	cudaMalloc(&devCoord4Y, nregions * sizeof(double));
	cudaMalloc(&devCoord4Z, nregions * sizeof(double));
	cudaMalloc(&devA, nregions * sizeof(double));
	cudaMalloc(&devB, nregions * sizeof(double));
	cudaMalloc(&devC, nregions * sizeof(double));
	cudaMalloc(&devD, nregions * sizeof(double));
	cudaMalloc(&devRotateFlag, nregions * sizeof(int));
	cudaMalloc(&devStartFlag, nregions * sizeof(int));
	cudaMalloc(&devEndFlag, nregions * sizeof(int));
	cudaMalloc(&devStart, nregions * sizeof(int));
	cudaMalloc(&devStableStart, nregions * sizeof(int));
	cudaMalloc(&devEnd, nregions * sizeof(int));
	cudaMalloc(&devStableEnd, nregions * sizeof(int));
	cudaMalloc(&devOmegaTarget, nregions * sizeof(double));
	cudaMalloc(&devCoord101X, nregions * sizeof(double));
	cudaMalloc(&devCoord101Y, nregions * sizeof(double));
	cudaMalloc(&devCoord101Z, nregions * sizeof(double));
	cudaMalloc(&devCoord102X, nregions * sizeof(double));
	cudaMalloc(&devCoord102Y, nregions * sizeof(double));
	cudaMalloc(&devCoord102Z, nregions * sizeof(double));
	cudaMalloc(&devAxisNormX, nregions * sizeof(double));
	cudaMalloc(&devAxisNormY, nregions * sizeof(double));
	cudaMalloc(&devAxisNormZ, nregions * sizeof(double));
	//
	hostStyle = (int *)malloc(nregions * sizeof(int));
	hostRadiusCylinder = (double *)malloc(nregions * sizeof(double));
	hostHeight = (double *)malloc(nregions * sizeof(double));
	hostCoord1X = (double *)malloc(nregions * sizeof(double));
	hostCoord1Y = (double *)malloc(nregions * sizeof(double));
	hostCoord1Z = (double *)malloc(nregions * sizeof(double));
	hostCoord2X = (double *)malloc(nregions * sizeof(double));
	hostCoord2Y = (double *)malloc(nregions * sizeof(double));
	hostCoord2Z = (double *)malloc(nregions * sizeof(double));
	hostCoord3X = (double *)malloc(nregions * sizeof(double));
	hostCoord3Y = (double *)malloc(nregions * sizeof(double));
	hostCoord3Z = (double *)malloc(nregions * sizeof(double));
	hostCoord4X = (double *)malloc(nregions * sizeof(double));
	hostCoord4Y = (double *)malloc(nregions * sizeof(double));
	hostCoord4Z = (double *)malloc(nregions * sizeof(double));
	hostA = (double *)malloc(nregions * sizeof(double));
	hostB = (double *)malloc(nregions * sizeof(double));
	hostC = (double *)malloc(nregions * sizeof(double));
	hostD = (double *)malloc(nregions * sizeof(double));
	hostRotateFlag = (int *)malloc(nregions * sizeof(int));
	hostStartFlag = (int *)malloc(nregions * sizeof(int));
	hostEndFlag = (int *)malloc(nregions * sizeof(int));
	hostStart = (int *)malloc(nregions * sizeof(int));
	hostStableStart = (int *)malloc(nregions * sizeof(int));
	hostEnd = (int *)malloc(nregions * sizeof(int));
	hostStableEnd = (int *)malloc(nregions * sizeof(int));
	hostOmegaTarget = (double *)malloc(nregions * sizeof(double));
	hostCoord101X = (double *)malloc(nregions * sizeof(double));
	hostCoord101Y = (double *)malloc(nregions * sizeof(double));
	hostCoord101Z = (double *)malloc(nregions * sizeof(double));
	hostCoord102X = (double *)malloc(nregions * sizeof(double));
	hostCoord102Y = (double *)malloc(nregions * sizeof(double));
	hostCoord102Z = (double *)malloc(nregions * sizeof(double));
	hostAxisNormX = (double *)malloc(nregions * sizeof(double));
	hostAxisNormY = (double *)malloc(nregions * sizeof(double));
	hostAxisNormZ = (double *)malloc(nregions * sizeof(double));
	//
	std::string keyCylinder = "cylinder";
	std::string keyPolygon = "polygon";
	std::cout<<"nregions="<<nregions<<std::endl;
	for (int i=0; i<nregions; i++) {
		// std::cout<<"********************************** regions[i]->style="<<regions[i]->style<<std::endl;
		if (!strncmp(regions[i]->style,keyCylinder.c_str(),keyCylinder.size())){
			hostStyle[i] = 1;
			hostRadiusCylinder[i] = regions[i]->radius;
			hostHeight[i] = regions[i]->height;
			hostCoord1X[i] = regions[i]->coords[0][0];
			hostCoord1Y[i] = regions[i]->coords[0][1];
			hostCoord1Z[i] = regions[i]->coords[0][2];
			hostCoord2X[i] = regions[i]->coords[1][0];
			hostCoord2Y[i] = regions[i]->coords[1][1];
			hostCoord2Z[i] = regions[i]->coords[1][2];
			hostA[i] = regions[i]->normal[0];
			hostB[i] = regions[i]->normal[1];
			hostC[i] = regions[i]->normal[2];
			// std::cout<<"hostStyle[i]="<<hostStyle[i]<<", hostRadiusCylinder[i]="<<hostRadiusCylinder[i]<<", hostHeight[i]="<<hostHeight[i]<<std::endl;
			// std::cout<<"hostCoord1[i]="<<hostCoord1X[i]<<", "<<hostCoord1Y[i]<<", "<<hostCoord1Z[i]<<std::endl;
			// std::cout<<"hostCoord2[i]="<<hostCoord2X[i]<<", "<<hostCoord2Y[i]<<", "<<hostCoord2Z[i]<<std::endl;
			// std::cout<<"hostA[i]="<<hostA[i]<<", "<<hostB[i]<<", "<<hostC[i]<<std::endl;
		} else if (!strncmp(regions[i]->style,keyPolygon.c_str(),keyPolygon.size())){
			hostStyle[i] = 2;
			hostCoord1X[i] = regions[i]->coords[0][0];
			hostCoord1Y[i] = regions[i]->coords[0][1];
			hostCoord1Z[i] = regions[i]->coords[0][2];
			hostCoord2X[i] = regions[i]->coords[1][0];
			hostCoord2Y[i] = regions[i]->coords[1][1];
			hostCoord2Z[i] = regions[i]->coords[1][2];
			hostCoord3X[i] = regions[i]->coords[2][0];
			hostCoord3Y[i] = regions[i]->coords[2][1];
			hostCoord3Z[i] = regions[i]->coords[2][2];
			hostCoord4X[i] = regions[i]->coords[3][0];
			hostCoord4Y[i] = regions[i]->coords[3][1];
			hostCoord4Z[i] = regions[i]->coords[3][2];
			hostA[i] = regions[i]->a;
			hostB[i] = regions[i]->b;
			hostC[i] = regions[i]->c;
			hostD[i] = regions[i]->d;
			hostRotateFlag[i] = regions[i]->rotate_flag;
			if (hostRotateFlag[i]==1) {
				hostStartFlag[i] = regions[i]->start_flag;
				hostEndFlag[i] = regions[i]->stop_flag;
				hostStart[i] = regions[i]->start;
				hostStableStart[i] = regions[i]->stable_start;
				hostEnd[i] = regions[i]->stop;
				hostStableEnd[i] = regions[i]->stable_end;
				hostOmegaTarget[i] = regions[i]->omega_target;
				hostCoord101X[i] = regions[i]->rot_axis->coords[0][0];
				hostCoord101Y[i] = regions[i]->rot_axis->coords[0][1];
				hostCoord101Z[i] = regions[i]->rot_axis->coords[0][2];
				hostCoord102X[i] = regions[i]->rot_axis->coords[1][0];
				hostCoord102Y[i] = regions[i]->rot_axis->coords[1][1];
				hostCoord102Z[i] = regions[i]->rot_axis->coords[1][2];
				hostAxisNormX[i] = regions[i]->rot_axis->normal[0];
				hostAxisNormY[i] = regions[i]->rot_axis->normal[1];
				hostAxisNormZ[i] = regions[i]->rot_axis->normal[2];
			}
			// std::cout<<"hostStyle[i]="<<hostStyle[i]<<std::endl;
			// std::cout<<"hostCoord1[i]="<<hostCoord1X[i]<<", "<<hostCoord1Y[i]<<", "<<hostCoord1Z[i]<<std::endl;
			// std::cout<<"hostCoord2[i]="<<hostCoord2X[i]<<", "<<hostCoord2Y[i]<<", "<<hostCoord2Z[i]<<std::endl;
			// std::cout<<"hostCoord3[i]="<<hostCoord3X[i]<<", "<<hostCoord3Y[i]<<", "<<hostCoord3Z[i]<<std::endl;
			// std::cout<<"hostCoord4[i]="<<hostCoord4X[i]<<", "<<hostCoord4Y[i]<<", "<<hostCoord4Z[i]<<std::endl;
			// std::cout<<"hostA[i]="<<hostA[i]<<", "<<hostB[i]<<", "<<hostC[i]<<", "<<hostD[i]<<std::endl;
			// std::cout<<"Flag[i]="<<hostRotateFlag[i]<<", "<<hostStartFlag[i]<<", "<<hostEndFlag[i]<<std::endl;
			// std::cout<<"Start[i]="<<hostStart[i]<<", "<<hostStableStart[i]<<std::endl;
			// std::cout<<"End[i]="<<hostEnd[i]<<", "<<hostStableEnd[i]<<std::endl;
			// std::cout<<"hostOmegaTarget[i]="<<hostOmegaTarget[i]<<", "<<hostOmegaTarget[i]<<std::endl;
			// std::cout<<"hostCoord101[i]="<<hostCoord101X[i]<<", "<<hostCoord101Y[i]<<", "<<hostCoord101Z[i]<<std::endl;
			// std::cout<<"hostCoord102[i]="<<hostCoord102X[i]<<", "<<hostCoord102Y[i]<<", "<<hostCoord102Z[i]<<std::endl;
			// std::cout<<"hostAxisNorm[i]="<<hostAxisNormX[i]<<", "<<hostAxisNormY[i]<<", "<<hostAxisNormZ[i]<<std::endl;
		} else {
			hostStyle[i] = 0;
		}
	}
	//
	cudaMemcpy(devStyle, hostStyle, nregions * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devRadiusCylinder, hostRadiusCylinder, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devHeight, hostHeight, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord1X, hostCoord1X, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord1Y, hostCoord1Y, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord1Z, hostCoord1Z, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord2X, hostCoord2X, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord2Y, hostCoord2Y, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord2Z, hostCoord2Z, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord3X, hostCoord3X, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord3Y, hostCoord3Y, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord3Z, hostCoord3Z, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord4X, hostCoord4X, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord4Y, hostCoord4Y, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord4Z, hostCoord4Z, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devA, hostA, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devB, hostB, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devC, hostC, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devD, hostD, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devRotateFlag, hostRotateFlag, nregions * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devStartFlag, hostStartFlag, nregions * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devEndFlag, hostEndFlag, nregions * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devStart, hostStart, nregions * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devStableStart, hostStableStart, nregions * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devEnd, hostEnd, nregions * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devStableEnd, hostStableEnd, nregions * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(devOmegaTarget, hostOmegaTarget, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord101X, hostCoord101X, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord101Y, hostCoord101Y, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord101Z, hostCoord101Z, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord102X, hostCoord102X, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord102Y, hostCoord102Y, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devCoord102Z, hostCoord102Z, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devAxisNormX, hostAxisNormX, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devAxisNormY, hostAxisNormY, nregions * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(devAxisNormZ, hostAxisNormZ, nregions * sizeof(double), cudaMemcpyHostToDevice);
	// std::cout<<"nregions="<<nregions<<std::endl;
	// exit(1);
}
