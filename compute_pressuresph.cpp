/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "compute_pressuresph.h"
#include "error.h"
#include "modify.h"
#include "particle.h"
#include "update.h"
#include "math.h"
#include "domain.h"
#include "parallel.h"

using namespace PDPS_NS;
#define PI 3.1416
/* ---------------------------------------------------------------------- */

ComputePressuresph::ComputePressuresph(PDPS *ps, int narg, char **arg) : Compute(ps, narg, arg)
{
	if (narg < 10) error->all(FLERR,"Illegal compute pressuresph command");

	scalar_flag = 1;
	pressuresph = 0.0;
	xtemp = atof(arg[3]);
	ytemp = atof(arg[4]);
	ztemp = atof(arg[5]);
	if (strcmp(arg[6], "Cubic") == 0)
		cubic_flag = 1;
	else if (strcmp(arg[6], "Quintic") == 0)
		quintic_flag = 1;
	rho_ref = atof(arg[7]);
	soundspeed = atof(arg[8]);
	h = atof(arg[9]);
	if (cubic_flag == 1){
		a2D = 10.0 / 7.0 / PI / h / h;
		a3D = 1.0 / PI / h / h / h;
	}
	else if (quintic_flag == 1){
		a2D = 7.0 / 4.0 / PI / h / h;
		a3D = 21.0 / 16.0 / PI / h / h / h;
	}

}

/* ---------------------------------------------------------------------- */

ComputePressuresph::~ComputePressuresph()
{

}

/* ---------------------------------------------------------------------- */

void ComputePressuresph::init()
{

}

/* ----------------------------------------------------------------------
   compute SPH pressure at the given point among given cutoff length
------------------------------------------------------------------------- */

double ComputePressuresph::compute_scalar()
{
	int i, j, k;
	double temp;
	double **x = particle->x;
	int *mask = particle->mask;
	int nlocal = particle->nlocal;
	double *rho = particle->rho;
	double q, rij, wf, p;
	double cutsq = h * h * 4;
	pressuresph = 0.0;
	wfsum = 0.0;
	double wfsumall = 0.0;
	for (int i = 0; i < nlocal; i++){
		if (mask[i] & groupbit) {
			rij = (x[i][0] - xtemp) * (x[i][0] - xtemp) + (x[i][1] - ytemp) * (x[i][1] - ytemp) + (x[i][2] - ztemp) * (x[i][2] - ztemp);
			if (rij < cutsq){
				wf = 0.0;
				q = sqrt(rij) / h;
				if (cubic_flag == 1){
					if (q < 1)
						wf = 1 - 1.5 * q * q + 0.75 * q * q * q;
					else
						wf = 0.25 * (2 - q) * (2 - q) * (2 - q);
				}
				else if (quintic_flag == 1)
					wf = (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (1 - q / 2.0) * (2 * q + 1);

				if (domain->dim == 3)
					wf = wf * a3D;
				else
					wf = wf * a2D;
				p = soundspeed * soundspeed * rho_ref / 7.0 * (pow((rho[i] / rho_ref), 7) - 1);
				wfsum += wf;
				pressuresph += p * wf;

			}
		}
	}
	
	//	communicate between processors to get the total average SPH pressure
	MPI_Allreduce(&wfsum, &wfsumall, 1, MPI_DOUBLE, MPI_SUM, mworld);
	MPI_Allreduce(&pressuresph, &scalar, 1, MPI_DOUBLE, MPI_SUM, mworld);

	if (wfsumall > 0.01)
		scalar /= wfsumall;
	else
		scalar = 0.0;

	return scalar;
}

