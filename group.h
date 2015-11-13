/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_GROUP_H
#define PS_GROUP_H

#include "pointers.h"

namespace PDPS_NS {

class Group : protected Pointers {
public:
    // group
    int ngroups;                    // # of defined groups
    char **name;                    // name of each group                       
    int *glocal;                    // number of particles in a group of one processor
	int *gparticles;                // total number of particles in a group
    int *bitmask;                   // bit mask for group
	int added_particles_all;        // added particles for group all

    Group(class PDPS *);
    ~Group();
    //void assign(int, char **);    // assign atoms to a group
    void assign(int, char**);       // add flagged atoms to a group
	void create(int, char**, int);
    int find_group(const char*);    // find the group id
	void update_all(int);

private:
    int procid;
    int find_unused();

	void create_dynamic();

	int *dynamic;
	int *nevery;

	int map_union(int, char **, int);
	int map_region(int, char **, int);
	int map_type(int, char **, int);

};

}

#endif
