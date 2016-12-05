/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "neighbor.h"
#include "create_particle.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "parallel.h"
#include "particle.h"
#include "particle_type.h"
#include "phy_const.h"
#include "random_park.h"
#include "style_particle.h"
#include "timer.h"

#include "pdps_cuda.h"
#include "cuda_engine.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
using namespace PDPS_NS;
using namespace PhyConst;

#define DELTA 10000
#define EPSILON 1.0e-6



__global__ void gputest(double *devCoordX){
	int pid = blockIdx.x * blockDim.x + threadIdx.x;
		devCoordX[pid] = -1;
}
//	convert from Aos to SoA
__global__ void gpuInterLeave(double *ArrayRaw, double *ArrayX, double *ArrayY, double * ArrayZ, const int nlocal){
	extern __shared__ double buf[];
	int pid = blockIdx.x * blockDim.x + threadIdx.x;
	if (pid < nlocal){
		buf[threadIdx.x * 3 + 0] = ArrayX[pid];
		buf[threadIdx.x * 3 + 1] = ArrayY[pid];
		buf[threadIdx.x * 3 + 2] = ArrayZ[pid];
	}
	__syncthreads();
	int g = blockIdx.x * blockDim.x * 3;
	int p = threadIdx.x;
	if (g + p < nlocal * 3)
		ArrayRaw[g + p] = buf[p];
	p += blockDim.x;
	if (g + p < nlocal * 3)
		ArrayRaw[g + p] = buf[p];
	p += blockDim.x;
	if (g + p < nlocal * 3)
		ArrayRaw[g + p] = buf[p];

}

//	convert from SoA to Aos
__global__ void gpuDeinterLeave(double *ArrayRaw, double *ArrayX, double *ArrayY, double * ArrayZ, const int nlocal){
	extern __shared__ double buf[];

	int g = blockIdx.x * blockDim.x * 3;
	int p = threadIdx.x;
	if (g + p < nlocal * 3)
		buf[p] = ArrayRaw[g + p];
	p += blockDim.x;
	if (g + p < nlocal * 3)
		buf[p] = ArrayRaw[g + p];
	p += blockDim.x;
	if (g + p < nlocal * 3)
		buf[p] = ArrayRaw[g + p];
	__syncthreads();

	int pid = blockIdx.x * blockDim.x + threadIdx.x;
	if (pid < nlocal){
		ArrayX[pid] = buf[threadIdx.x * 3 + 0];
		ArrayY[pid] = buf[threadIdx.x * 3 + 1];
		ArrayZ[pid] = buf[threadIdx.x * 3 + 2];
	}



}

// copy data between device
template<class TYPE> __global__ void gpuCopy(
	TYPE* __restrict Out,
	TYPE* __restrict In,
	const int  n)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n) Out[i] = In[i];
}


Particle::Particle(PDPS *ps) : Pointers(ps)
{


	x = NULL;
	v = NULL;
	f = NULL;
	tag = NULL;
	type = NULL;
	mass = NULL;
	mask = NULL;

	density = NULL;
	omega = NULL;
	radius = NULL;
	poro = NULL;
	volume = NULL;
	hlocal = NULL;
	rmass = NULL;
	torque = NULL;

	ptype = NULL;

	map_array = NULL;

	// Default value
	nparticles = 0;          // default number of particles
	nlocal = nghost = 0;
	nmax = 0;                // default number of materials
	nfirst = 0;
	ntypes = 0;              // default number of types
	//maxarg = 0;            // max argument to allocate memory

	vest = NULL;
	rho = NULL;
	drho = NULL;
	e = NULL;
	de = NULL;
	cv = NULL;

	// particle type flag
	atomic_flag = 1;
	sphere_flag = 0;
	rmass_flag = radius_flag = omega_flag = torque_flag = 0;
	ee_flag = rho_flag = cv_flag = vest_flag = 0;

	tag_enable = 1;
	map_style = 0;
	map_tag_max = 0;
	map_nhash = 0;

	// used by read_data class
	size_data_atom = 5;
	size_data_vel = 4;
	xcol_data = 3;

	nprimes = 38;
	primes = new int[nprimes];
	int plist[] = {5041,10007,20011,30011,40009,50021,60013,70001,80021,
				   90001,100003,110017,120011,130003,140009,150001,160001,
				   170003,180001,190027,200003,210011,220009,230003,240007,
				   250007,260003,270001,280001,290011,300007,310019,320009,
				   330017,340007,350003,362881,3628801};
	for (int i = 0; i < nprimes; i++) primes[i] = plist[i];

	particle_style = NULL;
	ptype = NULL;
	create_particle_type("atomic", 0, NULL);

	//	pointer to GPU device
	devCoordX = NULL;
	devCoordY = NULL;
	devCoordZ = NULL;
	devCoordXold = NULL;
	devCoordYold = NULL;
	devCoordZold = NULL;
	devVeloX = NULL;
	devVeloY = NULL;
	devVeloZ = NULL;
	devVestX = NULL;
	devVestY = NULL;
	devVestZ = NULL;
	devForceX = NULL;
	devForceY = NULL;
	devForceZ = NULL;
	devMask = NULL;
	devMass = NULL;
	devTag = NULL;
	devType = NULL;
	devRho = NULL;
	devRadius = NULL;
	devRmass = NULL;
	devDensity = NULL;
	devPoro = NULL;
	devVolume = NULL;

	devHostCoord = NULL;
	devHostVelo = NULL;
	devHostVest = NULL;
	devHostForce = NULL;
	devHostMask = NULL;
	devHostMassType = NULL;
	devHostTag = NULL;
	devHostType = NULL;
	devHostRho = NULL;
	devHostRadius = NULL;
	devHostRmass = NULL;
	devHostDensity = NULL;
	devHostPoro = NULL;
	devHostVolume = NULL;

	ptrHostCoord = NULL;
	ptrHostVelo = NULL;
	ptrHostVest = NULL;
	ptrHostForce = NULL;
	ptrHostMask = NULL;
	ptrHostMassType = NULL;
	ptrHostTag = NULL;
	ptrHostType = NULL;
	ptrHostRho = NULL;
	ptrHostRadius = NULL;
	ptrHostRmass = NULL;
	ptrHostDensity = NULL;
	ptrHostPoro = NULL;
	ptrHostVolume = NULL;
		

}

/* ---------------------------------------------------------------------- */

Particle::~Particle()
{
	memory->destroy(x);
	memory->destroy(v);
	memory->destroy(f);
	memory->destroy(type);
	memory->destroy(mask);

	delete[] mass;
	mass = NULL;
}

/* ---------------------------------------------------------------------- */

void Particle::init()
{
	if (nparticles == 0) {
		error->all(FLERR,"No particle has been created");
	}

  	ptype->init();
}

/* ---------------------------------------------------------------------- */

void Particle::create_particle_type(const char *style, int narg, char **arg)
{
	delete [] particle_style;
	if (ptype) delete ptype;

	if (0) return;

#define PARTICLE_CLASS
#define ParticleStyle(key,Class) \
	else if (strcmp(style,#key) == 0) ptype = new Class(ps,narg,arg);
#include "style_particle.h"
#undef ParticleStyle
#undef PARTICLE_CLASS

	else error->all(FLERR, "Invalid particle style");

	int n = strlen(style) + 1;
	particle_style = new char[n];
	strcpy(particle_style, style);
}



/* ----------------------------------------------------------------------
   Add tag for each created particle
------------------------------------------------------------------------- */

void Particle::add_tag()
{
	int maxtag = 0;
	for (int i = 0; i < nlocal; i++) 
		maxtag = MAX(maxtag,tag[i]);
	int maxtag_all;
	MPI_Allreduce(&maxtag,&maxtag_all,1,MPI_INT,MPI_MAX,mworld);

	// notag = # of particles with tag = 0 of one processor
	// notag_sum = total # of particles with tag = 0

	int notag = 0;
	for (int i = 0; i < nlocal; i++) {
		if (tag[i] == 0) {
			notag++;
		}
	}
	int notag_sum;
	MPI_Scan(&notag,&notag_sum,1,MPI_INT,MPI_SUM,mworld);

	// itag = 1st new tag that each processor should use

	int itag = maxtag_all + notag_sum - notag + 1;
	for (int i = 0; i < nlocal; i++) {
		if (tag[i] == 0) {
			tag[i] = itag;
			itag++;
		}
	}
}

/* ----------------------------------------------------------------------
   set a mass and flag it as set
   called from reading of data file
------------------------------------------------------------------------- */

void Particle::set_mass(const char *str)
{
	if (mass == NULL) error->all(FLERR,"Cannot set mass for this atom style");

	int itype;
	double mass_one;
	int n = sscanf(str,"%d %lg",&itype,&mass_one);
	if (n != 2) error->all(FLERR,"Invalid mass line in data file");

	if (itype < 1 || itype > ntypes)
		error->all(FLERR,"Invalid type for mass set");

	mass[itype] = mass_one;
	//mass_setflag[itype] = 1;

	if (mass[itype] <= 0.0) error->all(FLERR,"Invalid mass value");
}

/* ----------------------------------------------------------------------
   Set mass for each type of particle
------------------------------------------------------------------------- */

void Particle::set_mass(int narg, char** arg)
{
	int tid;

	// Need to tell if box exists
	if(mass == NULL) {
		allocate_type_arrays();
	}

	tid = atoi(arg[0]);                  // type id
	if (rmass_flag == 1) {
		double rm = atof(arg[1]);
		for (int i = 0; i < nlocal; i++) {
			if (type[i] == tid) rmass[i] = rm;
		}
	}
	else mass[tid] = atof(arg[1]);            // store mass

}

/* ----------------------------------------------------------------------
   Set density for each type of particle
------------------------------------------------------------------------- */

void Particle::set_density(int narg, char** arg)
{
	int tid;
	int gid = group->find_group(arg[0]);
	if (gid == -1) error->all(FLERR, "Cannot find the group id");

	int groupbit = group->bitmask[gid];

	if (sphere_flag == 0) error->all(FLERR, "Illegal particle style to call density command");

	if (!strcmp(arg[1], "set")) {
		if (narg != 3) error->all(FLERR, "Illegal radius command");
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				density[i] = atof(arg[2]);
			}
		}
	}


//	tid = atoi(arg[0]);
//	if (tid < 1) error->all(FLERR, "Illegal particle type");
	// Need to tell if box exists
//	if(density == NULL) {
//		allocate_type_arrays();
//	}

//	density[tid] = atof(arg[1]);
}

/* ----------------------------------------------------------------------
   Set radius for each type of particle
------------------------------------------------------------------------- */

void Particle::set_radius(int narg, char** arg)
{	
	if (sphere_flag == 0) error->all(FLERR, "Particle style is not correct to call radius command");
	if (narg < 1) error->all(FLERR, "Illegal radius command");

	int gid = group->find_group(arg[0]);
	if (gid == -1) error->all(FLERR, "Cannot find the group id");

	int groupbit = group->bitmask[gid];

	if (!strcmp(arg[1], "set")) {
		if (narg != 3) error->all(FLERR, "Illegal radius command");
		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & groupbit) {
				radius[i] = atof(arg[2]);
				volume[i] = 4.0 / 3 * PI * radius[i] * radius[i] * radius[i];
				rmass[i] = density[i] * volume[i];
				
			}
		}
	}
	else if (!strcmp(arg[1], "create")) {
		if (narg != 9) error->all(FLERR, "Illegal radius command");
		/* The following needs to be changed to global gaussian distribution
		RanPark *random;
		int seed;
		double rlo, rhi, rmean, rsigma;
		rlo = atof(arg[2]);
		rhi = atof(arg[3]);
		rmean = atof(arg[4]);
		rsigma = atof(arg[5]);
		seed = atoi(arg[4]);
		random = new RanPark(ps, seed);
		double num;
		int count = 0;
		for (int i = 0; i < nlocal; i++) {
			num = rlo - 1;
			count = 0;
			while (num < rlo || num > rhi) {
				num = (random->gaussian())*rsigma + rmean;
				count++;
				if (count > 100000) {
					error->all(FLERR, "Cannot generate raidus for the required distribution");
				}
			}
			radius[i] = num;
		}
		*/
	}
	else error->all(FLERR, "Illegal raidus command");
}

/* ----------------------------------------------------------------------
Set energy for each type of particle
------------------------------------------------------------------------- */
void Particle::set_energy(int narg, char** arg)
{
	int tid;

//	if (atomic_flag == 1) error->all(FLERR, "Illegal particle style to call density command");
	if (narg != 2) error->all(FLERR, "Illegal density command");

//	tid = atoi(arg[0]);
//	if (tid < 1) error->all(FLERR, "Illegal particle type");
	// Need to tell if box exists
	int gid = group->find_group(arg[0]);
	if (gid == -1) error->all(FLERR, "Cannot find the group id");
	int groupbit = group->bitmask[gid];
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			e[i] = atof(arg[1]);
		}
	}
}

/* ----------------------------------------------------------------------
Set rho for each type of particle
------------------------------------------------------------------------- */
void Particle::set_rho(int narg, char** arg)
{
	int tid;
	//	if (atomic_flag == 1) error->all(FLERR, "Illegal particle style to call density command");
	if (narg != 2) error->all(FLERR, "Illegal density command");

	//	tid = atoi(arg[0]);
	//	if (tid < 1) error->all(FLERR, "Illegal particle type");
	// Need to tell if box exists
	int gid = group->find_group(arg[0]);
	if (gid == -1) error->all(FLERR, "Cannot find the group id");
	int groupbit = group->bitmask[gid];
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			rho[i] = atof(arg[1]);
		}
	}
}
/* ----------------------------------------------------------------------
   allocate arrays of length ntypes
   only done after ntypes is set
------------------------------------------------------------------------- */

void Particle::allocate_type_arrays()
{
	//if (avec->mass_type) {
	if (rmass_flag) {
		density = new double[ntypes+1];
	}
	else {
		mass = new double[ntypes+1];
	}
}

/* ----------------------------------------------------------------------
   unpack n lines from Atom section of data file
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Particle::data_particles(int n, char *buf)
{
	int m,imagedata,xptr,iptr;
	double xdata[3],lamda[3];
	double *coord;
	char *next;

	next = strchr(buf,'\n');
	*next = '\0';
	int nwords = count_words(buf);
	*next = '\n';

	if (nwords != size_data_atom && nwords != size_data_atom + 3) {
		error->all(FLERR,"Incorrect atom format in data file");
	}

	char **values = new char*[nwords];

	// set bounds for my proc
	// if periodic and I am lo/hi proc, adjust bounds by EPSILON
	// insures all data atoms will be owned even with round-off

	double epsilon[3];

    epsilon[0] = domain->boxle[0] * EPSILON;
    epsilon[1] = domain->boxle[1] * EPSILON;
    epsilon[2] = domain->boxle[2] * EPSILON;
  
	double sublo[3],subhi[3];
 
    sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
    sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
    sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
  
	if (domain->xperiodic) {
		if (parallel->procloc[0] == 0) sublo[0] -= epsilon[0];
		if (parallel->procloc[0] == parallel->procgrid[0]-1) subhi[0] += epsilon[0];
	}
	if (domain->yperiodic) {
		if (parallel->procloc[1] == 0) sublo[1] -= epsilon[1];
		if (parallel->procloc[1] == parallel->procgrid[1]-1) subhi[1] += epsilon[1];
	}
	if (domain->zperiodic) {
		if (parallel->procloc[2] == 0) sublo[2] -= epsilon[2];
		if (parallel->procloc[2] == parallel->procgrid[2]-1) subhi[2] += epsilon[2];
	}

	// xptr = which word in line starts xyz coords
	// iptr = which word in line starts ix,iy,iz image flags

	xptr = xcol_data - 1;

	// loop over lines of atom data
	// tokenize the line into values
	// extract xyz coords and image flags
	// remap atom into simulation box
	// if atom is in my sub-domain, unpack its values

	for (int i = 0; i < n; i++) {
		next = strchr(buf,'\n');

		values[0] = strtok(buf," \t\n\r\f");
		if (values[0] == NULL)
		  error->all(FLERR,"Incorrect atom format in data file");
		for (m = 1; m < nwords; m++) {
		  values[m] = strtok(NULL," \t\n\r\f");
		  if (values[m] == NULL)
			error->all(FLERR,"Incorrect atom format in data file");
		}

		xdata[0] = atof(values[xptr]);
		xdata[1] = atof(values[xptr+1]);
		xdata[2] = atof(values[xptr+2]);
		//domain->remap(xdata,imagedata);
		coord = xdata;

		if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
			coord[1] >= sublo[1] && coord[1] < subhi[1] &&
			coord[2] >= sublo[2] && coord[2] < subhi[2])
			ptype->data_particle(xdata,values);

		buf = next + 1;
	}

	delete [] values;
}

/* ----------------------------------------------------------------------
   unpack n lines from Velocity section of data file
   check that atom IDs are > 0 and <= map_tag_max
   call style-specific routine to parse line
------------------------------------------------------------------------- */

void Particle::data_vels(int n, char *buf)
{
	int j,m,tagdata;
	char *next;

	next = strchr(buf,'\n');
	*next = '\0';
	int nwords = count_words(buf);
	*next = '\n';

	if (nwords != size_data_vel)
		error->all(FLERR,"Incorrect velocity format in data file");

	char **values = new char*[nwords];

	// loop over lines of atom velocities
	// tokenize the line into values
	// if I own atom tag, unpack its values

	for (int i = 0; i < n; i++) {
		next = strchr(buf,'\n');

		values[0] = strtok(buf," \t\n\r\f");
		for (j = 1; j < nwords; j++)
			values[j] = strtok(NULL," \t\n\r\f");

		tagdata = atoi(values[0]);
		if (tagdata <= 0 || tagdata > map_tag_max)
			error->one(FLERR,"Invalid atom ID in Velocities section of data file");
		if ((m = map(tagdata)) >= 0) {
			ptype->data_vel(m,&values[1]);
		}

		buf = next + 1;
	}

	delete [] values;
}

/* ----------------------------------------------------------------------
   count and return words in a single line
   make copy of line before using strtok so as not to change line
   trim anything from '#' onward
------------------------------------------------------------------------- */

int Particle::count_words(const char *line)
{
	int n = strlen(line) + 1;
	char *copy;
	memory->create(copy,n,"atom:copy");
	strcpy(copy,line);

	char *ptr;
	if (ptr = strchr(copy,'#')) *ptr = '\0';

	if (strtok(copy," \t\n\r\f") == NULL) {
		memory->destroy(copy);
		return 0;
	}
	n = 1;
	while (strtok(NULL," \t\n\r\f")) n++;

	memory->destroy(copy);
	return n;
}

/* ----------------------------------------------------------------------
   Allocate and initialize array table for global -> local map
   set map_tag_max = largest atom ID (may be larger than natoms)
   for array option:
     array length = 1 to largest tag of any atom
     set entire array to -1 as initial values
------------------------------------------------------------------------- */

void Particle::map_init()
{
	map_delete();

	if (tag_enable == 0) {
		error->all(FLERR,"Cannot create a particle map unless particles have IDs");
	}

	int max = 0;
	for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
	MPI_Allreduce(&max,&map_tag_max,1,MPI_INT,MPI_MAX,mworld);

	if (map_style == 1) {
		memory->create(map_array,map_tag_max+1,"atom:map_array");
		for (int i = 0; i <= map_tag_max; i++) map_array[i] = -1;
	}
	else {
		// map_nhash = max of atoms/proc or total atoms, times 2, at least 1000

		int nper = static_cast<int> (nparticles/parallel->nprocs);
		map_nhash = MAX(nper, nmax);
		if (map_nhash > nparticles) map_nhash = static_cast<int> (nparticles);
		if (parallel->nprocs > 1) map_nhash *= 2;
		map_nhash = MAX(map_nhash,1000);

		// map_nbucket = prime just larger than map_nhash

		int n = map_nhash/10000;
		n = MIN(n, nprimes-1);
		map_nbucket = primes[n];
		if (map_nbucket < map_nhash && n < nprimes-1) map_nbucket = primes[n+1];

		// set all buckets to empty
		// set hash to map_nhash in length
		// put all hash entries in free list and point them to each other

		map_bucket = new int[map_nbucket];
		for (int i = 0; i < map_nbucket; i++) map_bucket[i] = -1;

		map_hash = new HashElem[map_nhash];
		map_nused = 0;
		map_free = 0;
		for (int i = 0; i < map_nhash; i++) map_hash[i].next = i+1;
		map_hash[map_nhash-1].next = -1;
	}
}

/* ----------------------------------------------------------------------
   Clear global -> local map for all of my own and ghost atoms
------------------------------------------------------------------------- */

void Particle::map_clear()
{
	if (map_style == 1) {
		int nall = nlocal + nghost;
		for (int i = 0; i < nall; i++) map_array[tag[i]] = -1;
	} 
	else {
		int previous,global,ibucket,index;
		int nall = nlocal + nghost;
		for (int i = 0; i < nall; i++) {
			// search for key
			// if don't find it, done

			previous = -1;
			global = tag[i];
			ibucket = global % map_nbucket;
			index = map_bucket[ibucket];
			while (index > -1) {
				if (map_hash[index].global == global) break;
				previous = index;
				index = map_hash[index].next;
			}
			if (index == -1) continue;

			// delete the hash entry and add it to free list
			// special logic if entry is 1st in the bucket

			if (previous == -1) map_bucket[ibucket] = map_hash[index].next;
			else map_hash[previous].next = map_hash[index].next;

			map_hash[index].next = map_free;
			map_free = index;
			map_nused--;
		}
	}
	
}

/* ----------------------------------------------------------------------
   set global -> local map for all of my own and ghost atoms
   loop in reverse order so that nearby images take precedence over far ones
     and owned atoms take precedence over images
   this enables valid lookups of bond topology atoms
------------------------------------------------------------------------- */

void Particle::map_set()
{
	if (map_style == 1) {
		int nall = nlocal + nghost;
		for (int i = nall-1; i >= 0 ; i--) map_array[tag[i]] = i;
	} 
	else {
		int previous,global,ibucket,index;
		int nall = nlocal + nghost;
		if (nall > map_nhash) map_init();

		for (int i = nall-1; i >= 0 ; i--) {
		    // search for key
			// if found it, just overwrite local value with index

			previous = -1; 
			global = tag[i];
			ibucket = global % map_nbucket;
			index = map_bucket[ibucket];
			while (index > -1) {
				if (map_hash[index].global == global) break;
				previous = index;
				index = map_hash[index].next;
			}
			if (index > -1) {
				map_hash[index].local = i;
				continue;
			}

			// take one entry from free list
			// add the new global/local pair as entry at end of bucket list
			// special logic if this entry is 1st in bucket

			index = map_free;
			map_free = map_hash[map_free].next;
			if (previous == -1) map_bucket[ibucket] = index;
			else map_hash[previous].next = index;
			map_hash[index].global = global;
			map_hash[index].local = i;
			map_hash[index].next = -1;
			map_nused++;
		}
	}
}

/* ----------------------------------------------------------------------
   set global to local map for one atom
   for hash table option:
     global ID may already be in table if atom was already set
------------------------------------------------------------------------- */

void Particle::map_one(int global, int local)
{
	if (map_style == 1) map_array[global] = local;
	else {
		// search for key
		// if found it, just overwrite local value with index

		int previous = -1;
		int ibucket = global % map_nbucket;
		int index = map_bucket[ibucket];
		while (index > -1) {
			if (map_hash[index].global == global) break;
			previous = index;
			index = map_hash[index].next;
		}
		if (index > -1) {
			map_hash[index].local = local;
			return;
		}

		// take one entry from free list
		// add the new global/local pair as entry at end of bucket list
		// special logic if this entry is 1st in bucket

		index = map_free;
		map_free = map_hash[map_free].next;
		if (previous == -1) map_bucket[ibucket] = index;
		else map_hash[previous].next = index;
		map_hash[index].global = global;
		map_hash[index].local = local;
		map_hash[index].next = -1;
		map_nused++;
	}
}

/* ----------------------------------------------------------------------
   Free the array table for global to local mapping
------------------------------------------------------------------------- */

void Particle::map_delete()
{
	if (map_style == 1) {
		if (map_tag_max) memory->destroy(map_array);
	}
	else {
		if (map_nhash) {
			delete [] map_bucket;
			delete [] map_hash;
		}
		map_nhash = 0;
	}
	
	map_tag_max = 0;
}

/* ----------------------------------------------------------------------
   lookup global ID in hash table, return local index
------------------------------------------------------------------------- */

int Particle::map_find_hash(int global)
{
	int local = -1;
	int index = map_bucket[global % map_nbucket];
	while (index > -1) {
		if (map_hash[index].global == global) {
			local = map_hash[index].local;
			break;
		}
		index = map_hash[index].next;
	}
	return local;
}

/* ----------------------------------------------------------------------
   Set position for each particle
------------------------------------------------------------------------- */

void Particle::set_pos(int narg, char** arg)
{                


}

/* ---------------------------------------------------------------------- */

void Particle::lost_check()
{
	int ntotal;

	MPI_Allreduce(&nlocal, &ntotal, 1, MPI_INT, MPI_SUM, mworld);

	if (ntotal < nparticles) {
		char str[128];
		sprintf(str, "Particle lost from %d total particles to %d", particle->nparticles, ntotal);
		error->warning(FLERR, str);
	}
}

/* ----------------------------------------------------------------------
Delete particle
------------------------------------------------------------------------- */

void Particle::delete_particle(int n)
{
	int nlocal = particle->nlocal;
	int i;
	if (n < nlocal){
		tag[n] = tag[nlocal];
		type[n] = type[nlocal]; 
		mask[n] = mask[nlocal];
		for (i = 0; i < 3; i++) {
			x[n][i] = x[nlocal][i];
			v[n][i] = v[nlocal][i];
			f[n][i] = f[nlocal][i];
		
		}
	}
	if (atomic_flag == 1) 
		mass[n] = mass[nlocal];
	else {
		for (i = 0; i < 3; i++){
			omega[n][i] = omega[nlocal][i];
			torque[n][i] = torque[nlocal][i];
		}
		radius[n] = radius[nlocal];
		rmass[n] = rmass[nlocal];
		poro[n] = poro[nlocal];
		volume[n] = volume[nlocal];
		hlocal[n] = volume[nlocal];
	}

	particle->nlocal--;
}

/* ----------------------------------------------------------------------
Save particle
------------------------------------------------------------------------- */

void Particle::save_particle(int narg, char** arg)
{
	if (!strcmp(arg[0], "cylinder")) {
		if (narg != 6) error->all(FLERR, "Illegal save command");
		double coord[3], height, radius;
		coord[0] = atof(arg[1]);
		coord[1] = atof(arg[2]);
		coord[2] = atof(arg[3]);
		height = atof(arg[4]);
		radius = atof(arg[5]);
		int nlocal = particle->nlocal;
		int i;
		double rijsq;
		for (i = nlocal; i >= 0; i--) {
			rijsq = (x[i][0] - coord[0]) * (x[i][0] - coord[0]) + (x[i][1] - coord[1]) * (x[i][1] - coord[1]);
			if(rijsq < radius * radius && x[i][2] > coord[2] && x[i][2] < coord[2] + height)
				continue;
			else
				delete_particle(i);
		}
	}
}

//	transfer data to easy fetch memory
/* ---------------------------------------------------------------------- */
void Particle::PinHostArray(){
	cudaDeviceSynchronize();
	if (nmax){
			ptrHostCoord = cudaEngine->PinHost(devHostCoord, &(particle->x[0][0]), 3 * particle->nmax * sizeof(double));
			ptrHostForce = cudaEngine->PinHost(devHostForce, &(particle->f[0][0]), 3 * particle->nmax * sizeof(double));
			ptrHostVelo = cudaEngine->PinHost(devHostVelo, &(particle->v[0][0]), 3 * particle->nmax * sizeof(double));
			ptrHostVest = cudaEngine->PinHost(devHostVest, &(particle->vest[0][0]), 3 * particle->nmax * sizeof(double));
			ptrHostMassType = cudaEngine->PinHost(devHostMassType, particle->mass, (particle->ntypes + 1) * sizeof(double));
			ptrHostType = cudaEngine->PinHost(devHostType, particle->type, particle->nmax * sizeof(int));
			ptrHostTag = cudaEngine->PinHost(devHostTag, particle->tag, particle->nmax * sizeof(int));
			ptrHostMask = cudaEngine->PinHost(devHostMask, particle->mask, particle->nmax * sizeof(int));
	}
	if (rho) ptrHostRho = cudaEngine->PinHost(devHostRho, particle->rho, particle->nmax * sizeof(double));
	if (rmass)	ptrHostRmass = cudaEngine->PinHost(devHostRmass, particle->rmass, particle->nmax * sizeof(double));
	if (density) ptrHostDensity = cudaEngine->PinHost(devHostDensity, particle->density, particle->nmax * sizeof(double));
	if (radius) ptrHostRadius = cudaEngine->PinHost(devHostRadius, particle->radius, particle->nmax * sizeof(double));
	if (poro) ptrHostPoro = cudaEngine->PinHost(devHostPoro, particle->poro, particle->nmax * sizeof(double));
	if (volume) ptrHostVolume = cudaEngine->PinHost(devHostVolume, particle->volume, particle->nmax * sizeof(double));

}

/* ---------------------------------------------------------------------- */
void Particle::UnpinHostArray(){
	cudaDeviceSynchronize();
	if (ptrHostCoord)    cudaEngine->UnpinHost(&(particle->x[0][0]), ptrHostCoord, devHostCoord);
	if (ptrHostForce)    cudaEngine->UnpinHost(&(particle->f[0][0]), ptrHostForce, devHostForce);
	if (ptrHostVelo)     cudaEngine->UnpinHost(&(particle->v[0][0]), ptrHostVelo, devHostVelo);
	if (ptrHostVest)     cudaEngine->UnpinHost(&(particle->vest[0][0]), ptrHostVest, devHostVest);
	if (ptrHostRho)    cudaEngine->UnpinHost(particle->rho, ptrHostRho, devHostRho);
	if (ptrHostType)    cudaEngine->UnpinHost(particle->type, ptrHostType, devHostType);
	if (ptrHostTag)    cudaEngine->UnpinHost(particle->tag, ptrHostTag, devHostTag);
	if (ptrHostMask)    cudaEngine->UnpinHost(particle->mask, ptrHostMask, devHostMask);
	if (ptrHostMassType) cudaEngine->UnpinHost(particle->mass, ptrHostMassType, devHostMassType);
	if (ptrHostRmass) cudaEngine->UnpinHost(particle->rmass, ptrHostRmass, devHostRmass);
	if (ptrHostRadius) cudaEngine->UnpinHost(particle->radius, ptrHostRadius, devHostRadius);
	if (ptrHostDensity) cudaEngine->UnpinHost(particle->density, ptrHostDensity, devHostDensity);
	if (ptrHostPoro) cudaEngine->UnpinHost(particle->poro, ptrHostPoro, devHostPoro);
	if (ptrHostVolume) cudaEngine->UnpinHost(particle->volume, ptrHostVolume, devHostVolume);

}

//	transfer data From CPU to GPU
/* ---------------------------------------------------------------------- */
void Particle::TransferC2G(){

	enum cudaLimit 	limit;
	size_t *memory, *free, *total;
	memory = (size_t *)malloc(sizeof(size_t));
	free = (size_t *)malloc(sizeof(size_t));
	total = (size_t *)malloc(sizeof(size_t));
	int *device_num;
	device_num = (int *)malloc(sizeof(int));
	cudaError_t cudaStatus = cudaGetDevice(device_num);
	cudaStatus = cudaMemGetInfo(free, total);
	//cudaStatus = cudaDeviceGetLimit(memory, limit);



	int nlocal = particle->nlocal;
	if (nlocal == 0)
		return;

	const int BLK = 512;
	//	Therr are maximum 16 stream to run the code transfer concurrently
	int nStream = cudaEngine->StreamPool.size();
	vector<cudaStream_t> StreamPool = cudaEngine->StreamPool;
	vector<cudaEvent_t> Events;

	//	Transfer each terms on different streams and record the event
	cudaStream_t &Stream = StreamPool[0];
	cudaMemcpyAsync(devHostCoord, ptrHostCoord, 3 * nlocal * sizeof(double), cudaMemcpyHostToDevice, Stream);
	gpuDeinterLeave << < int(nlocal + BLK - 1) / BLK, BLK, BLK * 3 * sizeof(double), Stream >> >(
		devHostCoord, devCoordX, devCoordY, devCoordZ, nlocal);
	Events.push_back(cudaEngine->Event("X_C2G_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[1];
	cudaMemcpyAsync(devHostVelo, ptrHostVelo, 3 * nlocal * sizeof(double), cudaMemcpyHostToDevice, Stream);
	gpuDeinterLeave << < int(nlocal + BLK - 1) / BLK, BLK, BLK * 3 * sizeof(double), Stream >> >(
		devHostVelo, devVeloX, devVeloY, devVeloZ, nlocal);
	Events.push_back(cudaEngine->Event("V_C2G_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[2];
	cudaMemcpyAsync(devHostForce, ptrHostForce, 3 * nlocal * sizeof(double), cudaMemcpyHostToDevice, Stream);
	gpuDeinterLeave << < int(nlocal + BLK - 1) / BLK, BLK, BLK * 3 * sizeof(double), Stream >> >(
		devHostForce, devForceX, devForceY, devForceZ, nlocal);
	Events.push_back(cudaEngine->Event("F_C2G_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[3];
	cudaMemcpyAsync(devHostType, ptrHostType, nlocal * sizeof(int), cudaMemcpyHostToDevice, Stream);
	cudaMemcpyAsync(devHostMask, ptrHostMask, nlocal * sizeof(int), cudaMemcpyHostToDevice, Stream);
	cudaMemcpyAsync(devHostTag, ptrHostTag, nlocal * sizeof(int), cudaMemcpyHostToDevice, Stream);
	cudaMemcpyAsync(devHostMassType, ptrHostMassType, (ntypes + 1) * sizeof(double), cudaMemcpyHostToDevice, Stream);
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devType, devHostType, nlocal);
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devMask, devHostMask, nlocal);
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devTag, devHostTag, nlocal);
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devMass, devHostMassType, ntypes + 1);
	Events.push_back(cudaEngine->Event("T_C2G_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[4];
	cudaMemcpyAsync(devHostRadius, ptrHostRadius, nlocal * sizeof(double), cudaMemcpyHostToDevice, Stream);
	cudaMemcpyAsync(devHostRmass, ptrHostRmass, nlocal * sizeof(double), cudaMemcpyHostToDevice, Stream);
	cudaMemcpyAsync(devHostRho, ptrHostRho, nlocal * sizeof(double), cudaMemcpyHostToDevice, Stream);
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devRadius, devHostRadius, nlocal);
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devRmass, devHostRmass, nlocal);
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devRho, devHostRho, nlocal);
	Events.push_back(cudaEngine->Event("R_C2G_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[5];
	cudaMemcpyAsync(devHostDensity, ptrHostDensity, nlocal * sizeof(double), cudaMemcpyHostToDevice, Stream);
	cudaMemcpyAsync(devHostPoro, ptrHostPoro, nlocal * sizeof(double), cudaMemcpyHostToDevice, Stream);
	cudaMemcpyAsync(devHostVolume, ptrHostVolume, nlocal * sizeof(double), cudaMemcpyHostToDevice, Stream);
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devDensity, devHostDensity, nlocal);
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devPoro, devHostPoro, nlocal);
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devVolume, devHostVolume, nlocal);
	Events.push_back(cudaEngine->Event("D_C2G_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[6];
	cudaMemcpyAsync(devHostVest, ptrHostVest, 3 * nlocal * sizeof(double), cudaMemcpyHostToDevice, Stream);
	gpuDeinterLeave << < int(nlocal + BLK - 1) / BLK, BLK, BLK * 3 * sizeof(double), Stream >> >(
		devHostVest, devVestX, devVestY, devVestZ, nlocal);
	Events.push_back(cudaEngine->Event("V_C2G_PREV"));
	cudaEventRecord(Events.back(), Stream);

	//	Make sure they all stop before end the function
	for (int i = 0; i < Events.size(); i++)
		cudaStreamWaitEvent(StreamPool[0], Events[i], 0);

}

__global__ void add(int *a, int *b, int *c) {
	c[threadIdx.x] = a[threadIdx.x] + b[threadIdx.x];
}

//	transfer data From GPU to CPU
/* ---------------------------------------------------------------------- */
void Particle::TransferG2C(){

	//enum cudaLimit 	limit;
	//size_t *memory, *free, *total;
	//memory = (size_t *)malloc(sizeof(size_t));
	//free = (size_t *)malloc(sizeof(size_t));
	//total = (size_t *)malloc(sizeof(size_t));
	//int *device_num;
	//device_num = (int *)malloc(sizeof(int));
	//cudaError_t cudaStatus = cudaGetDevice(device_num);
	//cudaStatus = cudaMemGetInfo(free, total);
	////cudaStatus = cudaDeviceGetLimit(memory, limit);
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	int nlocal = particle->nlocal;
	if (nlocal == 0)
		return;

	const int BLK = 512;

	//	Therr are maximum 16 stream to run the code transfer concurrently
	int nStream = cudaEngine->StreamPool.size();
	vector<cudaStream_t> StreamPool = cudaEngine->StreamPool;
	vector<cudaEvent_t> Events;

	//	Transfer each terms on different streams and record the event
	cudaStream_t &Stream = StreamPool[0];
	gpuInterLeave << < int(nlocal + BLK - 1) / BLK, BLK, BLK * 3 * sizeof(double), Stream >> >(
		devHostCoord, devCoordX, devCoordY, devCoordZ, nlocal); 
	cudaMemcpyAsync(ptrHostCoord, devHostCoord, 3 * nlocal * sizeof(double), cudaMemcpyDeviceToHost, Stream);
	Events.push_back(cudaEngine->Event("X_G2C_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[1];
	gpuInterLeave << < int(nlocal + BLK - 1) / BLK, BLK, BLK * 3 * sizeof(double), Stream >> >(
		devHostVelo, devVeloX, devVeloY, devVeloZ, nlocal);
	cudaMemcpyAsync(ptrHostVelo, devHostVelo, 3 * nlocal * sizeof(double), cudaMemcpyDeviceToHost, Stream);
	Events.push_back(cudaEngine->Event("V_G2C_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[2];
	gpuInterLeave << < int(nlocal + BLK - 1) / BLK, BLK, BLK * 3 * sizeof(double), Stream >> >(
		devHostForce, devForceX, devForceY, devForceZ, nlocal);
	cudaMemcpyAsync(ptrHostForce, devHostForce, 3 * nlocal * sizeof(double), cudaMemcpyDeviceToHost, Stream);
	Events.push_back(cudaEngine->Event("F_G2C_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[3];
	cudaMemcpyAsync(ptrHostType, devType, nlocal * sizeof(int), cudaMemcpyDeviceToHost, Stream);
	cudaMemcpyAsync(ptrHostMask, devMask, nlocal * sizeof(int), cudaMemcpyDeviceToHost, Stream);
	cudaMemcpyAsync(ptrHostTag, devTag, nlocal * sizeof(int), cudaMemcpyDeviceToHost, Stream);
	cudaMemcpyAsync(ptrHostMassType, devMass, nlocal * sizeof(double), cudaMemcpyDeviceToHost, Stream);
	Events.push_back(cudaEngine->Event("T_G2C_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[4];
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devRadius, devHostRadius, nlocal);
	gpuCopy << < int(nlocal + BLK - 1) / BLK, BLK, 0, Stream >> > (devRmass, devHostRmass, nlocal);
	cudaMemcpyAsync(ptrHostRadius, devRadius, nlocal * sizeof(double), cudaMemcpyDeviceToHost, Stream);
	cudaMemcpyAsync(ptrHostRmass, devRmass, nlocal * sizeof(double), cudaMemcpyDeviceToHost, Stream);
	cudaMemcpyAsync(ptrHostRho, devRho, nlocal * sizeof(double), cudaMemcpyDeviceToHost, Stream);
	Events.push_back(cudaEngine->Event("R_G2C_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[5];
	cudaMemcpyAsync(ptrHostDensity, devDensity, nlocal * sizeof(double), cudaMemcpyDeviceToHost, Stream);
	cudaMemcpyAsync(ptrHostPoro, devPoro, nlocal * sizeof(double), cudaMemcpyDeviceToHost, Stream);
	cudaMemcpyAsync(ptrHostVolume, devVolume, nlocal * sizeof(double), cudaMemcpyDeviceToHost, Stream);
	Events.push_back(cudaEngine->Event("D_G2C_PREV"));
	cudaEventRecord(Events.back(), Stream);

	Stream = StreamPool[6];
	gpuInterLeave << < int(nlocal + BLK - 1) / BLK, BLK, BLK * 3 * sizeof(double), Stream >> >(
		devHostVest, devVestX, devVestY, devVestZ, nlocal);
	cudaMemcpyAsync(ptrHostVest, devHostVest, 3 * nlocal * sizeof(double), cudaMemcpyDeviceToHost, Stream);
	Events.push_back(cudaEngine->Event("V_G2C_PREV"));
	cudaEventRecord(Events.back(), Stream);

	//	Make sure they all stop before end the function
	for (int i = 0; i < Events.size(); i++)
		cudaStreamWaitEvent(StreamPool[0], Events[i], 0);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	time = time;
}