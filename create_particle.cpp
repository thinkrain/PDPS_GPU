/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "create_particle.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "lattice.h"
#include "memory.h"
#include "neighbor.h"
#include "output.h"
#include "parallel.h"
#include "particle.h"
#include "particle_type.h"
#include "psmath.h"
#include "random_park.h"
#include "region.h"

using namespace PDPS_NS;
using namespace PsMath_NS;

#define GAP 1.0e-6           // define small GAP (GAP = GAP*spacing) set at two ends of each edge
#define MAX_TRY 10000        // maximum # of try to get valid random radius 

enum{ATOMIC, SINGLE, UNIFORM, GAUSSIAN};

CreateParticle::CreateParticle(PDPS *ps) : Pointers (ps) {}

/* ---------------------------------------------------------------------- */

void CreateParticle::command(int narg, char **arg) 
{
	if (domain->box_exist == 0) {
		error->all(FLERR,"Create_particle command before simulation box is defined");
	}
	random_no_overlap_flag = 0;

	dist_style = ATOMIC;
	lattice = NULL;
	random_radius = NULL;
	nparticles_previous = particle->nparticles;

	for (int i = 0; i < 3; i++) {
		sublo[i] = domain->sublo[i];
		subhi[i] = domain->subhi[i];
	}

	int iarg = 0;
	while (iarg < narg) {
		if (!strcmp(arg[iarg+1], "single")) {
			tid = atoi(arg[iarg]);
			double x = atof(arg[iarg+2]);
			double y = atof(arg[iarg+3]);
			double z = atof(arg[iarg+4]);
			if (particle->radius_flag) {
				if (narg != 6) error->all(FLERR, "Particle radius is needed for this type of \"particle_style\"");
				rsingle = atof(arg[iarg+5]);
			}
			create_single(x, y, z);
			iarg += 6;
		}
		else if (strcmp(arg[iarg+1], "number_density") == 0) {
			tid = atoi(arg[iarg]);
			rho = atof(arg[iarg+2]);
			rid = domain->find_region(arg[iarg+3]);
			if (rid == -1) error->all(FLERR, "Cannot find the region");
			if (particle->radius_flag) {
				if (narg <= 4) error->all(FLERR, "Particle radius is needed for this type of \"particle_style\"");
				if (!strcmp(arg[iarg+4], "single")) {
					dist_style = SINGLE;
					rsingle = atof(arg[iarg+5]);
					iarg += 6;
				}
				else if (!strcmp(arg[iarg+4], "uniform")) {
					dist_style = UNIFORM;
					rlo = atof(arg[iarg+5]);
					rhi = atof(arg[iarg+6]);
					seed = atoi(arg[iarg+7]);
					iarg += 8;
				}
				else if (!strcmp(arg[iarg+4], "gaussian")) {
					dist_style = GAUSSIAN;
					rlo = atof(arg[iarg+5]);
					rhi = atof(arg[iarg+6]);
					rmean = atof(arg[iarg+7]);
					rsigma = atof(arg[iarg+8]);
					seed = atoi(arg[iarg+9]);
					iarg += 10;
				}
			}
			create_number_density(tid, rid);
			iarg += 4;
		}
		else if (!strcmp(arg[iarg+1], "random_no_overlap") || !strcmp(arg[iarg+1], "random")) {
			tid = atoi(arg[iarg]);
			if (!strcmp(arg[iarg+1], "random_no_overlap")) random_no_overlap_flag = 1;
			else random_no_overlap_flag = 0;
			spacing = atof(arg[iarg+2]);
			n_target = atoi(arg[iarg+3]);
			rid = domain->find_region(arg[iarg+4]);
			if (rid == -1) error->all(FLERR, "Cannot find the region");
			if (particle->radius_flag) {
				if (!strcmp(arg[iarg+5], "single")) {
					dist_style = SINGLE;
					rsingle = atof(arg[iarg+6]);
					iarg += 7;
				}
				else if (!strcmp(arg[iarg+5], "uniform")) {
					iarg++;
					dist_style = UNIFORM;
					rlo = atof(arg[iarg+6]);
					rhi = atof(arg[iarg+7]);
					iarg += 8;
				}
				else if (!strcmp(arg[iarg+5], "gaussian")) {
					dist_style = GAUSSIAN;
					rlo = atof(arg[iarg+6]);
					rhi = atof(arg[iarg+7]);
					rmean = atof(arg[iarg+8]);
					rsigma = atof(arg[iarg+9]);
					iarg += 10;
				}
			}
			else {
				dist_style = ATOMIC;
			}
			seed = atoi(arg[iarg++]);
			iarg += 6;
			if (lattice == NULL) {
				char **lattice_arg = new char*[3];
				for (int i = 0; i < 3; i++) lattice_arg[i] = new char[128];
				strcpy(lattice_arg[0], domain->regions[rid]->name);
				sprintf(lattice_arg[1], "spacing");
				sprintf(lattice_arg[2], "%g", spacing);
				lattice = new Lattice(ps, 3, lattice_arg);
				for (int i = 0; i < 3; i++) delete[] lattice_arg[i];
				delete[] lattice_arg;
			}
			create_random();
		}
		else if (strcmp(arg[iarg + 1], "spacing") == 0){
			tid = atoi(arg[iarg]);
			rid = domain->find_region(arg[iarg + 2]);
			create_spacing(tid, rid);
			iarg += 3;
		}
		else error->all(FLERR, "Illegal create_particle command");
	}
		
	delete lattice;

	bigint nlocal = particle->nlocal;
	MPI_Allreduce(&nlocal, &particle->nparticles, 1, MPI_INT, MPI_SUM, mworld);
	
	nparticles_now = particle->nparticles - nparticles_previous;

	// Automatically update the number of particles in group "all"

	if (nparticles_now == 0) {
		error->all(FLERR,"No particle has been created, please check your input");
	}

	group->update_all(nparticles_now);
	
	if (particle->nparticles < 0 || particle->nparticles > (int) pow(2.0,30)) {
		error->all(FLERR,"Too many total particles");
	}

	if (parallel->procid == 0) {
		char str[128];
		sprintf(str,"%d of particles of type %d are created\n\n", nparticles_now, tid);
		output->print(str);
	}

	particle->add_tag();
}

/* ----------------------------------------------------------------------
   Creat single particle
---------------------------------------------------------------------- */

void CreateParticle::create_single(double x, double y, double z)
{ 
	double coord[3];

	coord[0] = x;
	coord[1] = y;
	coord[2] = z;

	if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
		coord[1] >= sublo[1] && coord[1] < subhi[1] &&
		coord[2] >= sublo[2] && coord[2] < subhi[2]) {
		if (particle->radius_flag) particle->ptype->create_particle(tid, coord, rsingle);
		else particle->ptype->create_particle(tid, coord);
	}
}

/* ----------------------------------------------------------------------
   Creat particles based on the lattice spacing (uniform)
---------------------------------------------------------------------- */

void CreateParticle::create_spacing(int tid, int rid)
{   
	double coord[3];

	xle = domain->regions[rid]->extent_xle;     // length along x coordinate
	xlo = domain->regions[rid]->extent_xlo;     // the lower bound along x coordinate
	xhi = domain->regions[rid]->extent_xhi;     // the upper bound along x coordinate
	yle = domain->regions[rid]->extent_yle;
	ylo = domain->regions[rid]->extent_ylo;
	yhi = domain->regions[rid]->extent_yhi;
	zle = domain->regions[rid]->extent_zle;          
	zlo = domain->regions[rid]->extent_zlo;    
	zhi = domain->regions[rid]->extent_zhi;

//	dx = domain->lsp[0];                        // store the lattice spacing
//	dy = domain->lsp[1];
//	dz = domain->lsp[2];
	dx = domain->lattice->cle[0];
	dy = domain->lattice->cle[1];
	dz = domain->lattice->cle[2];

	npx = int((xle - 2*GAP*xle)/dx) + 1;        // In order to avoid particle overlap on the box boundary, 
	npy =  int((yle - 2*GAP*yle)/dy) + 1;       // there are two small GAPs 1/4*dx set at two ends of each edge
	npz =  int((zle - 2*GAP*zle)/dz) + 1;

	// Creat particles
	for (int k = 0; k < npz; k++) 
	for (int j = 0; j < npy; j++) 
	for (int i = 0; i < npx; i++) 
	{
		coord[0] = xlo + GAP*xle + i*dx;
		coord[1] = ylo + GAP*yle + j*dy;
		coord[2] = zlo + GAP*zle + k*dz;

		// tell which processor it belongs to
		if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
			coord[1] >= sublo[1] && coord[1] < subhi[1] &&
			coord[2] >= sublo[2] && coord[2] < subhi[2]) {
			if (particle->radius_flag) particle->ptype->create_particle(tid, coord, rsingle);
			else particle->ptype->create_particle(tid, coord);
		}
	} 
}

/* ----------------------------------------------------------------------
   Creat particles based on the number of density (uniform)
---------------------------------------------------------------------- */

void CreateParticle::create_number_density(int tid, int rid)
{ 
	int i, j;
	double length, area, vol;
	int factors[3];
	double **coords, **rot_coords;
	double inv_rot[3][3];
	int nps;
	int n_possible, n_target;
	int dim = domain->dim;
	double radi;

	int *cell_num;
	int ncxy, ncxyz;
	int cox, coy, coz;
	int ncoxyz;
	double *cle;
	int *nc;

	Region *region = domain->regions[rid];

	extent_xle = region->extent_xle;     // length along x coordinate
	extent_xlo = region->extent_xlo;     // the lower bound along x coordinate
	extent_xhi = region->extent_xhi;     // the upper bound along x coordinate
	extent_yle = region->extent_yle;
	extent_ylo = region->extent_ylo;
	extent_yhi = region->extent_yhi;
	extent_zle = region->extent_zle;          
	extent_zlo = region->extent_zlo;    
	extent_zhi = region->extent_zhi;

	if (dist_style == UNIFORM || dist_style == GAUSSIAN) random_radius = new RanPark(ps, seed);

	// ----------- create particles on a line  --------------

	if (region->line_flag == 1) {
		double *normal = region->normal;
		coords = region->coords;
		length = region->length;

		n_possible = static_cast<int> (length*rho);
		n_target = n_possible;

		double dr = length*(1 - 2*GAP) / (n_target - 1);
		double pro_coord;
		double coord[3];

		for (i = 0; i < n_target; i++) {
			pro_coord = GAP*length + i*dr;
			coord[0] = coords[0][0] + pro_coord*normal[0];
			coord[1] = coords[0][1] + pro_coord*normal[1];
			coord[2] = coords[0][2] + pro_coord*normal[2];

			if (dist_style == SINGLE) radi = rsingle;
			else if (dist_style == UNIFORM) radi = get_uniform_radius();
			else if (dist_style == GAUSSIAN) radi = get_gaussian_radius();

			if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
			coord[1] >= sublo[1] && coord[1] < subhi[1] &&
			coord[2] >= sublo[2] && coord[2] < subhi[2]) {
				if (particle->radius_flag) particle->ptype->create_particle(tid, coord, radi);
				else particle->ptype->create_particle(tid, coord);
			}
		}
	}
	// ----------- create particles on a plane  --------------

	else if (region->plane_flag == 1) {
		double *normal = region->normal;
		nps = region->nps;
		coords = region->coords;
		rot_coords = region->rot_coords;
		for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			inv_rot[i][j] = region->inv_rot[i][j];
		}

		xlo = xhi = rot_coords[0][0];
		ylo = yhi = rot_coords[0][1];
		zlo = zhi = rot_coords[0][2];

		for (i = 1; i < nps; i++) {
			xlo = MIN(xlo, rot_coords[i][0]);
			xhi = MAX(xhi, rot_coords[i][0]);
			xle = xhi - xlo;
			ylo = MIN(ylo, rot_coords[i][1]);
			yhi = MAX(yhi, rot_coords[i][1]);
			yle = yhi - ylo;
			zlo = MIN(zlo, rot_coords[i][2]);
			zhi = MAX(zhi, rot_coords[i][2]);
			zle = zhi - zlo;
		}

		area = xle * yle;

		n_possible = static_cast<int> (area*rho);	
		n_target = static_cast<int> (region->area*rho);
		
		find_factors2(n_possible, factors);

		// In order to avoid particle overlap on the box boundary
		dx = xle * (1 - 2*GAP) / (factors[0] - 1);
		dy = yle * (1 - 2*GAP) / (factors[1] - 1);
		dz = 0.0;
		
		npx = factors[0];        
		npy = factors[1];
		npz = factors[2];

		double pro_coord[3];
		double coord[3];
		int flag;
		int counter1 = 0;
		int counter2 = 0;

		for (int k = 0; k < npz; k++) {
			for (int j = 0; j < npy; j++) {
				for (int i = 0; i < npx; i++) {
					counter1++;

					pro_coord[0] = xlo + GAP*xle + i*dx;
					pro_coord[1] = ylo + GAP*yle + j*dy;
					pro_coord[2] = zlo + GAP*zle + k*dz;

					// tell if the point locates inside the projected domain
					flag = region->inside_a_plane(pro_coord, nps, rot_coords);
					if (flag == 0) continue;

					Matrix_Prod_3D(coord, inv_rot, pro_coord);
					// tell which processor it belongs to

					if (dist_style == SINGLE) radi = rsingle;
					else if (dist_style == UNIFORM) radi = get_uniform_radius();
					else if (dist_style == GAUSSIAN) radi = get_gaussian_radius();

					if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
					coord[1] >= sublo[1] && coord[1] < subhi[1] &&
					coord[2] >= sublo[2] && coord[2] < subhi[2]) {
						if (particle->radius_flag) particle->ptype->create_particle(tid, coord, radi);
						else particle->ptype->create_particle(tid, coord);
						counter2++;
					}
					if (counter1 >= n_possible || counter2 >= n_target) {
						break;
					}
				}
				if (counter1 >= n_possible || counter2 >= n_target)  break;
			}
			if (counter1 >= n_possible || counter2 >= n_target)  break;
		} // for (int k = 0; k < npz; k++)
	} // if (region->plane_flag == 1)
	// ----------- create particles in volume --------------

	else if (region->volume_flag == 1) {
		xle = region->extent_xle;     // length along x coordinate
		xlo = region->extent_xlo;     // the lower bound along x coordinate
		xhi = region->extent_xhi;     // the upper bound along x coordinate
		yle = region->extent_yle;
		ylo = region->extent_ylo;
		yhi = region->extent_yhi;
		zle = region->extent_zle;          
		zlo = region->extent_zlo;    
		zhi = region->extent_zhi;

		vol = xle * yle * zle;
		n_possible = static_cast<int> (vol*rho);
		n_target = static_cast<int> (region->volume*rho);
		
		find_factors3(n_possible, factors);

		// In order to avoid particle overlap on the box boundary
		dx = dy = dz = 0.0;
		if (factors[0] > 1) dx = xle * (1 - 2*GAP) / (factors[0] - 1);
		if (factors[1] > 1) dy = yle * (1 - 2*GAP) / (factors[1] - 1);
		if (factors[2] > 1) dz = zle * (1 - 2*GAP) / (factors[2] - 1);

		npx = factors[0];        
		npy = factors[1];
		npz = factors[2];

		double coord[3];
		int flag;
		int counter1 = 0;
		int counter2 = 0;

		int icx, icy, icz;
		int cid;
		int exist_flag;
		int iter = 0;
		double dist, n[3];
		for (int k = 0; k < npz; k++) {
			for (int j = 0; j < npy; j++) {
				for (int i = 0; i < npx; i++) {
					counter1++;

					coord[0] = xlo + GAP*xle + i*dx;
					coord[1] = ylo + GAP*yle + j*dy;
					coord[2] = zlo + GAP*zle + k*dz;
					
					// tell if the point locates inside the projected domain
					flag = region->inside(coord);
					if (flag == 0) continue;

					if (dist_style == SINGLE) radi = rsingle;
					else if (dist_style == UNIFORM) radi = get_uniform_radius();
					else if (dist_style == GAUSSIAN) radi = get_gaussian_radius();

					// tell which processor it belongs to
					if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
					coord[1] >= sublo[1] && coord[1] < subhi[1] &&
					coord[2] >= sublo[2] && coord[2] < subhi[2]) {
						if (particle->radius_flag) particle->ptype->create_particle(tid, coord, radi);
						else particle->ptype->create_particle(tid, coord);
					    counter2++;
					}

					if (counter1 >= n_possible || counter2 >= n_target) {
						break;
					}
				}
				if (counter1 >= n_possible || counter2 >= n_target)  break;
			}
			if (counter1 >= n_possible || counter2 >= n_target)  break;
		} // for (int k = 0; k < npz; k++)
	} // else if (region->volume_flag == 1)

	if (random_radius) delete random_radius;
}

/* ----------------------------------------------------------------------
   Creat particles based on the number of density (uniform)
---------------------------------------------------------------------- */

void CreateParticle::create_random()
{
	int i, j, k, m;
	double coord[3];
	int *cell_num;
	int ncxy, ncxyz;
	int cox, coy, coz;
	int ncoxyz;
	double *cle;
	int *nc;
	int *periodicity;
	int inside_flag;

	Region *region = domain->regions[rid];
	periodicity = domain->periodicity;

	ncxyz = lattice->ncxyz;
	ncxy = lattice->ncxy;
	nc = lattice->nc;

	cell_num = lattice->cell_num;
	cle = lattice->cle;

	double temp;
	if (dist_style == ATOMIC) temp = spacing;
	else if (dist_style == SINGLE) temp = rsingle;
	else temp = rmean;
	cox = temp / cle[0];
	if (cox*cle[0] < temp) cox++;
	coy = temp / cle[1];
	if (coy*cle[1] < temp) coy++;
	coz = temp / cle[2];
	if (coz*cle[2] < temp) coz++;
	if (domain->dim == 2) coz = 0;
	ncoxyz = 2 * cox * 2 * coy * 2 * coz;
	if (domain->dim == 2) ncoxyz = 2 * cox * 2 * coy;
	
	RanPark *random_coord;
	random_coord = new RanPark(ps, seed);
	if (dist_style == UNIFORM || dist_style == GAUSSIAN) random_radius = new RanPark(ps, seed);

	double xlo, xhi, xle, ylo, yhi, yle, zlo, zhi, zle;

	xlo = region->extent_xlo;
	xhi = region->extent_xhi;
	xle = region->extent_xle;
	ylo = region->extent_ylo;
	yhi = region->extent_yhi;
	yle = region->extent_yle;
	zlo = region->extent_zlo;
	zhi = region->extent_zhi;
	zle = region->extent_zle;

	//cell_coord = lattice->cell_coord;

	//int n_possible = nc[0] / 4.0 * nc[1] /4.0 * nc[2] / 4.0;
	//if (domain->dim == 2) n_possible = nc[0] / 4.0 * nc[1] / 4.0;
	//if (n_target > n_possible) error->all(FLERR, "The targeted # of particles to be created may be too many");

	int cid;
	double x, y, z;
	int exist_flag;
	int iter;
	int max_iter = 0;
	int icx, icy, icz;
	double radi, dist; 
	double n[3];
	for (m = 0; m < n_target; m++) {
		exist_flag = 1;
		iter = 0;
		if (dist_style == SINGLE) radi = rsingle;
		else if (dist_style == UNIFORM) radi = get_uniform_radius();
		else if (dist_style == GAUSSIAN) radi = get_gaussian_radius();
		while (exist_flag > 0) {
			x = random_coord->uniform();
			y = random_coord->uniform();
			z = random_coord->uniform();
			coord[0] = x*xle + xlo;
			coord[1] = y*yle + ylo;
			coord[2] = z*zle + zlo;
			// the following two lines need further investigated  !!!!!!!
//			if (particle->radius_flag) {
//				dist = domain->regions[rid]->find_interaction_distance(n, coord);
//				if (fabs(dist) < radi) continue;
//			}

			if (domain->dim == 2) coord[2] = 0.0;
			cid = lattice->coord2cell(coord, icx, icy, icz);

			if (dist_style == UNIFORM || dist_style == GAUSSIAN) {
				cox = radi / cle[0];
				if (cox*cle[0] < radi) cox++;
				coy = radi / cle[1];
				if (coy*cle[1] < radi) coy++;
				coz = radi / cle[2];
				if (coz*cle[2] < radi) coz++;
				if (domain->dim == 2) coz = 0;

				ncoxyz = 2 * cox * 2 * coy * 2 * coz;
				if (domain->dim == 2) ncoxyz = 2 * cox * 2 * coy;
			}

			// the range cannot exceed the box boundary
			if (icx-cox < 0 || icx+cox >= nc[0]) continue;
			if (icy-coy < 0 || icy+coy >= nc[1]) continue;
			if (icz-coz < 0 || icz+coz >= nc[2]) continue;
			
			if (domain->dim == 3) {
				for (k = icz - coz; k <= icz + coz; k++) {
					for (j = icy - coy; j <= icy + coy; j++) {
						for (i = icx - cox; i <= icx + cox; i++)
						{
							cid = k*ncxy + j*nc[0] + i;
							exist_flag = cell_num[cid];
							if (exist_flag > 0) {
								break;
							}
						}
						if (exist_flag > 0) {
							break;
						}
					}
					if (exist_flag > 0) {
						break;
					}
				}
			}
			else if (domain->dim == 2) {
				for (j = icy - coy; j <= icy + coy; j++) {
					for (i = icx - cox; i <= icx + cox; i++)
					{
						cid = j*nc[0] + i;
						exist_flag = cell_num[cid];
						if (exist_flag > 0) {
							break;
						}
					}
					if (exist_flag > 0) {
						break;
					}
				}
			}

			iter++;
			if (iter > 100000) {
				bigint nlocal = particle->nlocal;
				int n_now;
				MPI_Allreduce(&nlocal, &n_now, 1, MPI_INT, MPI_SUM, mworld);
				int n_created = n_now - nparticles_previous;
				char str[1024];
				sprintf(str, "Too many attempts to find a valid position to \
					          create a particle in ramdom_no_overlap style. \
							 Only %d particles have been created", n_created);
				error->all(FLERR, str);
			}
			/*coord[0] = cell_coord[cid][0];
			coord[1] = cell_coord[cid][1];
			coord[2] = cell_coord[cid][2];*/

			inside_flag = region->inside(coord);
			if (inside_flag == 0) {
				exist_flag = 1;
				continue;
			}
		}
		max_iter = MAX(iter, max_iter);
		
		if (random_no_overlap_flag) {
			// scan the offsets and mark the cell
			if (domain->dim == 3) {
				for (k = icz - coz; k <= icz + coz; k++)
				for (j = icy - coy; j <= icy + coy; j++)
				for (i = icx - cox; i <= icx + cox; i++) {
					cid = k*ncxy + j*nc[0] + i;
					cell_num[cid]++;
					if (cell_num[cid] > 1) {
						error->all(FLERR, "Some cell has already been marked before");
					}
				}
			}
			else if (domain->dim == 2) {
				for (j = icy - coy; j <= icy + coy; j++)
				for (i = icx - cox; i <= icx + cox; i++) {
					cid = j*nc[0] + i;
					cell_num[cid]++;
					if (cell_num[cid] > 1) {
						error->all(FLERR, "Some cell has already been marked before");
					}
				}
			}
		}
		// tell which processor it belongs to
		if (coord[0] >= sublo[0] && coord[0] < subhi[0] && 
			coord[1] >= sublo[1] && coord[1] < subhi[1] && 
			coord[2] >= sublo[2] && coord[2] < subhi[2]) {
			if (particle->radius_flag) particle->ptype->create_particle(tid, coord, radi);
			else particle->ptype->create_particle(tid, coord);
		}
	}

	// The following block is used to visualize the radius distribution
	
	/*FILE *ftest;
	ftest = fopen("radius_sample.txt", "w");
	int ntest = 30;
	int *test_count = new int[ntest];
	double *test_vol = new double[ntest];
	for (int i = 0; i < ntest; i++) {
		test_count[i] = 0;
		test_vol[i] = 0.0;
	}
	double test_space = (rhi - rlo) / ntest;
	int index;
	double v_p = 0.0;
	double test_temp;
	for (int i = 0; i < particle->nlocal; i++) {
		if (particle->type[i] == 3) {
			test_temp = 3.1415926*particle->radius[i] * particle->radius[i];
			index = (particle->radius[i] - rlo) / test_space;
			v_p += test_temp;
			if (index == ntest) index--;
			test_count[index]++;
			test_vol[index] += test_temp;
		}
	}
	fprintf(ftest, "Total volume of particles is %f\n", v_p);
	for (int i = 0; i < ntest; i++) {
		fprintf(ftest, "%f %d %f\n", rlo + test_space*i, test_count[i], test_vol[i]);
	}*/
	
	delete random_coord;
	random_coord = NULL;
	if (random_radius) delete random_radius;
	//delete lattice;
	//lattice = NULL;
}

/* ----------------------------------------------------------------------
   Find factors for 3D cases so that npx * npy * npz > number of 
   particles to be created and make sure npx * npy * npz is the minimum 
   among possible choices
---------------------------------------------------------------------- */

void CreateParticle::find_factors3(int n, int *nbest)
{ 
	int nmin[3];
	int np[3];
	int nx,ny,nz;
	double temp[3];
	double c[3];
	int closest,exceed;

	c[0] = 1.0;
	c[1] = yle/xle;
	c[2] = zle/xle;

	temp[0] = pow(n/(c[1]*c[2]),1.0/3);
	temp[1] = temp[0] * c[1];
	temp[2] = temp[0] * c[2];
	
	for (int i = 0; i < 3; i++) {
		nmin[i] = static_cast<int> (temp[i]);
	}

	if (nmin[0] * nmin[1] * nmin[2] == n) closest = 0;
	else closest = n;

	np[0] = nbest[0] = nmin[0];
	np[1] = nbest[1] = nmin[1];
	np[2] = nbest[2] = nmin[2]; 

	// find the best factors
	// 1st try add "1" to any one of the components
	for (int i = 0; i < 3; i++) {
		np[i] = nmin[i] + 1;
		exceed = np[0] * np[1] * np[2] - n;
		if (exceed >= 0 && exceed < closest) {
			closest = exceed;
			nbest[0] = np[0];
			nbest[1] = np[1];
			nbest[2] = np[2];
		}
		// reset
		np[i] = nmin[i];
	}

	// 2nd try add "1" to any two of the components

	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			np[i] = nmin[i] + 1;
			np[j] = nmin[j] + 1;
			exceed = np[0] * np[1] * np[2] - n;
			if (exceed >= 0 && exceed < closest) {
				closest = exceed;
				nbest[0] = np[0];
				nbest[1] = np[1];
				nbest[2] = np[2];
			}
			np[i] = nmin[i];
			np[j] = nmin[j];
		}
	}

	// 3rd try add "1" to all components
	
	for (int i = 0; i < 3; i++)
		np[i] = nmin[i] + 1;

	exceed = np[0] * np[1] * np[2] - n;
	if (exceed >= 0 && exceed < closest) {
		closest = exceed;
		nbest[0] = np[0];
		nbest[1] = np[1];
		nbest[2] = np[2];
	}
}

/* ----------------------------------------------------------------------
   Find factors for 2D cases so that npx * npy * npz > number of 
   particles to be created and make sure npx * npy * npz is the minimum 
   among possible choices
---------------------------------------------------------------------- */

void CreateParticle::find_factors2(int n, int *nbest)
{ 
	int nmin[3];
	int np[3];
	int nx,ny,nz;
	double temp[3];
	double c[3];
	int closest,exceed;

	c[0] = 1.0;
	c[1] = yle/xle;
	c[2] = 0.0;

	temp[0] = sqrt(n/c[1]); 
	temp[1] = temp[0] * c[1];
	temp[2] = temp[0] * c[2];
	
	for (int i = 0; i < 3; i++) {
		nmin[i] = static_cast<int> (temp[i]);
	}
	nmin[2] = 1;

	if (nmin[0] * nmin[1] * nmin[2] == n) closest = 0;
	else closest = n;

	np[0] = nbest[0] = nmin[0];
	np[1] = nbest[1] = nmin[1];
	np[2] = nbest[2] = nmin[2]; 

	// find the best factors
	// 1st try add "1" to any one of the components
	for (int i = 0; i < 2; i++) {
		np[i] = nmin[i] + 1;
		exceed = np[0] * np[1] * np[2] - n;
		if (exceed >= 0 && exceed < closest) {
			closest = exceed;
			nbest[0] = np[0];
			nbest[1] = np[1];
			nbest[2] = np[2];
		}
		// reset
		np[i] = nmin[i];
	}

	// 2nd try add "1" to any two of the components

	for (int i = 0; i < 2; i++) {
		for (int j = i + 1; j < 2; j++) {
			np[i] = nmin[i] + 1;
			np[j] = nmin[j] + 1;
			exceed = np[0] * np[1] * np[2] - n;
			if (exceed >= 0 && exceed < closest) {
				closest = exceed;
				nbest[0] = np[0];
				nbest[1] = np[1];
				nbest[2] = np[2];
			}
			np[i] = nmin[i];
			np[j] = nmin[j];
		}
	}

	// 3rd try add "1" to all components
	
	for (int i = 0; i < 2; i++)
		np[i] = nmin[i] + 1;

	exceed = np[0] * np[1] * np[2] - n;
	if (exceed >= 0 && exceed < closest) {
		closest = exceed;
		nbest[0] = np[0];
		nbest[1] = np[1];
		nbest[2] = np[2];
	}
}

/* ----------------------------------------------------------------------
get radius in a gaussian distribution
---------------------------------------------------------------------- */

double CreateParticle::get_uniform_radius()
{
	double radius;

	radius = random_radius->uniform()*(rhi - rlo) + rlo;

	return radius;
}

/* ----------------------------------------------------------------------
    get radius in a gaussian distribution
---------------------------------------------------------------------- */

double CreateParticle::get_gaussian_radius()
{
	double radius, temp;

	temp = random_radius->gaussian();
	radius = temp*rsigma + rmean;
	int count = 0;
	while (radius < rlo || radius > rhi) {
		radius = random_radius->gaussian()*rsigma + rmean;
		count++;
		if (count > MAX_TRY) error->all(FLERR, "Could not find a random radius in the specified range");
	}

	return radius;
}
