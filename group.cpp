/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "output.h"
#include "parallel.h"
#include "particle.h"
#include "region.h"
//#include "stdint.h"

using namespace PDPS_NS;
//using namespace std; 

#define MAXGROUP 31        // default maximum number of groups

Group::Group(PDPS *ps) : Pointers(ps)
{
	procid = parallel->procid;

	bitmask = NULL;
	name = NULL;
	glocal = NULL;
	gparticles = NULL;

	name = new char*[MAXGROUP];
	bitmask = new int[MAXGROUP];
	glocal = new int[MAXGROUP];
	gparticles = new int[MAXGROUP];

	for (int i = 0; i < MAXGROUP; i++) {
		name[i] = NULL;
		bitmask[i] = (1 << i);     // each group has bit "1" 
		glocal[i] = 0;
		gparticles[i] = 0;
	}
	
	// Create default group "all"
	ngroups = 0;    
	name[0] = new char[4];
	strcpy(name[ngroups],"all");	
	ngroups++;
}

/* ----------------------------------------------------------------------
					        Free All Memory
------------------------------------------------------------------------- */

Group::~Group()
{
	for (int i = 0; i < MAXGROUP; i++) {
		delete[] name[i];
		name[i] = NULL;
	}
	delete[] name;
	delete[] bitmask;
	delete[] glocal;
	delete[] gparticles;
	name = NULL;
	bitmask = NULL;
	glocal = NULL;
	gparticles = NULL;
}

/* ----------------------------------------------------------------------
   Create Group
------------------------------------------------------------------------- */

void Group::assign(int narg, char **arg)
{
	int igroup;
	int added_particles;
	int nparticles = particle->nparticles;
	
	added_particles = 0;

	//if (nparticles == 0) {
		//error->all(FLERR,"group command has to be after create_particle command"); 
	//}

	for (igroup = 0; igroup < ngroups; igroup++) {
		if(strcmp(arg[0],name[igroup]) == 0) break;
	}

	// create new group
	if (igroup == ngroups) {
		create(narg,arg,igroup);
	}

	// Default "all" group (no function here)
	if (narg == 1) {
		if (!strcmp(arg[0],"all") == 0) {
			error->all(FLERR,"The name of all group has to be 'all'");
		}
		if (procid == 0) {
			char str1[128];
			sprintf(str1,"%d particles added to group all\n",added_particles_all);
			char str2[128];
			sprintf(str2,"%d particles in group all now\n\n",gparticles[0]);
			output->print(str1);
			output->print(str2);
		}
	}

	// Create user-defined group
	if (narg > 1) {
		// Map region
		if (!strcmp(arg[1], "region")) {
			if (narg != 3) error->all(FLERR,"Illegal group region command");

			added_particles = map_region(narg,arg,igroup);
		} 

		// Combine multiple groups
		else if (!strcmp(arg[1], "union")) {
			if (narg == 2) error->all(FLERR,"Illega group type command");

			added_particles = map_union(narg,arg,igroup);
		} 

		// Map types
		else if (!strcmp(arg[1],"type")) {
	//		if (narg == 2) error->all(FLERR,"Illega group type command");

			added_particles = map_type(narg,arg,igroup);
		} 
		
		else if (!strcmp(arg[1], "dynamic")) {
			if (narg != 7) error->all(FLERR, "Illegal group dynamic command");
			
			int iarg = 3;
			while (iarg < narg) {
				if (!strcmp(arg[iarg], "region")) {
					map_region(narg, arg, igroup);
				}
				else if (!strcmp(arg[iarg], "every")) {
					
				}
				iarg += 2;
			}
			
		}
		else {
			error->all(FLERR,"Invalid group command");
		}

		glocal[igroup] += added_particles;

		bigint total_added_particles;
		total_added_particles = 0;
		MPI_Allreduce(&added_particles,&total_added_particles,1,MPI_INT,MPI_SUM,mworld);

		gparticles[igroup] += total_added_particles;
		
		//MPI_Allreduce(&glocal[igroup],&gparticles[igroup],1,MPI_INT,MPI_SUM,mworld);

		if (procid == 0) {
			char str1[128];
			sprintf(str1,"%d particles added to the group %s\n",total_added_particles,name[igroup]);
			char str2[128];
			sprintf(str2,"%d particles in group %s now\n\n",gparticles[igroup],name[igroup]);
			output->print(str1);
			output->print(str2);
		}
	} // if (narg > 1)
}

/* ---------------------------------------------------------------------- */

void Group::create(int narg, char **arg, int igroup)
{
	int n = strlen(arg[0]) + 1;
	name[igroup] = new char[n];
	strcpy(name[igroup],arg[0]);        // store group's name

	ngroups++;
}

/* ----------------------------------------------------------------------
   Update the group all
   ---------------------------------------------------------------------- */

void Group::update_all(int n) 
{
	added_particles_all = n;
	glocal[0] = particle->nlocal;
	gparticles[0] = particle->nparticles;
}

/* ---------------------------------------------------------------------- */

int Group::map_region(int narg, char **arg, int igroup)
{
	int nlocal = particle->nlocal;
	int *mask = particle->mask;
	double **x = particle->x;
	int rid;
	class Region **regions = domain->regions;
	int added = 0;

	rid = domain->find_region(arg[2]);
	if (rid == -1) {
		char str[128];
		sprintf(str, "Cannot find region %s", arg[2]);
		error->all(FLERR,str);
	}

	for (int i = 0; i < nlocal; i++) {
		if (regions[rid]->inside(x[i])) {
			if ((mask[i] & bitmask[igroup]) == 0) {
				mask[i] |= bitmask[igroup]; 
				added++;
			}
		}
	}

	return added;
}

/* ---------------------------------------------------------------------- */

int Group::map_type(int narg, char **arg, int igroup)
{
	int nlocal = particle->nlocal;
	int *mask = particle->mask;
	int *type = particle->type;
	int tid;
	int added = 0;

	for (int i = 0; i < nlocal; i++) {
		for (int iarg = 2; iarg < narg; iarg++) {
			tid = atoi(arg[iarg]);
			if (type[i] == tid) {
				if ((mask[i] & bitmask[igroup]) == 0) {
					mask[i] |= bitmask[igroup];
					added++;
				}
			}
		}
	} 

	return added;
}

/* ---------------------------------------------------------------------- */

int Group::map_union(int narg, char **arg, int igroup)
{
	int nlocal = particle->nlocal;
	int *mask = particle->mask;
	int gid;
	int added = 0;

	for (int i = 0; i < nlocal; i++) {
		for (int iarg = 2; iarg < narg; iarg++) {
			gid = find_group(arg[iarg]);
			if (gid == -1) {
				char str[128];
				sprintf(str, "Cannot find group %s", arg[iarg]);
				error->all(FLERR,str);
			}

			if ((mask[i] & bitmask[igroup]) == 0) { 
				if (bitmask[gid] == (mask[i] & bitmask[gid])) {
					mask[i] |= bitmask[igroup];
					added++;
				}
			}
		}
	}

	return added;
}

/* ----------------------------------------------------------------------
   Find group index: no match, return -1
------------------------------------------------------------------------- */

int Group::find_group(const char* gname)
{
	int gid;

	gid = -1;
	for (int i = 0; i < ngroups; i++) {
		if (!strcmp(gname, name[i])) {	
			gid = i;       // scan for the region's name
		}
	}

	return gid;
}
