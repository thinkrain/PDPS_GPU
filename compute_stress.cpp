/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "compute_stress.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "modify.h"
#include "pair.h"
#include "particle.h"
#include "style_compute.h"
#include "update.h"

using namespace PDPS_NS;

enum{IK};

/* ---------------------------------------------------------------------- */

ComputeStress::ComputeStress(PDPS *ps, int narg, char **arg) : Compute(ps, narg, arg)
{
	if (narg < 3) error->all(FLERR,"Illegal compute stress command");
	
	scalar_flag = 0;
	vector_flag = 1;

	for (int i = 0; i < 6; i++) {
		stress[i] = 0.0;
	}

	if (!strcmp(arg[3],"IK")) {
		stress_flag = IK;
	}
	nevery = atoi(arg[4]);

	vector = new double[6];

	//fix_ave_id = modify->find_fix("ave/spatial");
	//if (fix_ave_id == -1) error->all(FLERR,"Cannot find fix ave/spatial id");

	fname = NULL;
	int n = strlen(arg[5]) + 1;
	fname = new char[n];
	strcpy(fname,arg[5]);
	file = fopen(fname,"w");
	if (file == NULL) {
		char str[128];
		sprintf(str,"Cannot open file \"%s\"",fname);

	}
	counter = 0;
}

/* ---------------------------------------------------------------------- */

ComputeStress::~ComputeStress()
{

}

/* ---------------------------------------------------------------------- */

void ComputeStress::init()
{
	//fixavespatial = modify->fix[fix_ave_id];
	for (int i = 0; i < 6; i++) {
		vector[i] = 0.0;
	}
}

/* ---------------------------------------------------------------------- */

double ComputeStress::compute_scalar()
{


	return 0;
}

/* ----------------------------------------------------------------------
   compute stress vector
------------------------------------------------------------------------- */

void ComputeStress::compute_vector()
{
	int ntimestep = update->ntimestep;

	counter++;
	if (stress_flag == IK) {
		Irving_Kirkwood();
	}

	if (counter % nevery == 0) {
		for (int i = 0; i < 6; i++) {
			vector[i] /= nevery;
		}
		fprintf(file,"%d %f %f %f %f %f %f\n",ntimestep,vector[0],vector[1],vector[2],vector[3],vector[4],vector[5]);
		for (int i = 0; i < 6; i++) 
			vector[i] = 0.0;
	}
}

/* ----------------------------------------------------------------------
   Irving-Kirkwood method: 
   stress = -n*<sum(mi*vi \otimes vi + 0.5*sum(sum(rij \otimes fij))  >
------------------------------------------------------------------------- */

void ComputeStress::Irving_Kirkwood()
{
	int nparticles = particle->nparticles;
	int *mask = particle->mask;
	double *mass = particle->mass;
	double **v = particle->v;
	double **x = particle->x;
	//double v[3];
	int ibin;
	int i, k;
	double volume;

	for (int i = 0; i < 6; i++) {
		stress[i] = 0.0;
		virial[i] = 0.0;
	}

	for (i = 0; i < force->npairs; i++)
	for (k = 0; k < 6; k++) {
		virial[k] += force->pair[i]->virial[k];
	}

	//volume = (domain->boxle[0] * domain->boxle[1] * domain->boxle[2]);
	volume = 1000;
	stress[0] = 0.5 / volume * virial[0];
	stress[1] = 0.5 / volume * virial[1];
	stress[2] = 0.5 / volume * virial[2];
	stress[3] = 0.5 / volume * virial[3];
	stress[4] = 0.5 / volume * virial[4];
	stress[5] = 0.5 / volume * virial[5];
	vector[0] = stress[0];
	vector[1] = stress[1];
	vector[2] = stress[2];
	vector[3] = stress[3];
	vector[4] = stress[4];
	vector[5] = stress[5];

			//v[0] = fixavespatial-
			/*stress[0] += mass[i] * v[i][0] * v[i][0];
			stress[1] += mass[i] * v[i][1] * v[i][1];
			stress[2] += mass[i] * v[i][2] * v[i][2];
			stress[3] += mass[i] * v[i][0] * v[i][1];
			stress[4] += mass[i] * v[i][0] * v[i][2];
			stress[5] += mass[i] * v[i][1] * v[i][2];*/
			//stress[1] 
}

/* ----------------------------------------------------------------------
   locate bin id for each particle
------------------------------------------------------------------------- */

int ComputeStress::check_bins(int id) 
{
	double **x = particle->x;
	int icell;

	/*if (dim[2] == 1) {
		icell = int((x[id][2] - zlo)/dh[2]);
		if (icell < 0) icell = 0;
		if (icell >= nbins[2]) icell = nbins[2] - 1;
	}*/

	return icell;
}
