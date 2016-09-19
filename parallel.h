/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_PARALLEL_H
#define PS_PARALLEL_H

#include "pointers.h"

namespace PDPS_NS {

class Parallel : protected Pointers {
public:
	int procid, nprocs;                 // proc info
	int procgrid[3];                    // grid of processors (eg: 12 procs = 3 * 2 * 2)
	int user_procgrid[3];               // user reuqest for the grid of processors
	int procloc[3];                     // processor location (eg: 0, 0, 1)
	int **procfactors;                  // processor grid: pxyz = px * py * pz
	int procneigh[3][2];                // 6 neighbor processors, 0/1 = left/right
	int ***grid2proc;
	double *xsplit,*ysplit,*zsplit;     // fractional (0-1) sub-domain sizes

	double rcghost[3];                  // cutoffs used for acquiring ghost particles

	int ghost_velocity;					// 1 if ghost atoms have velocity, 0 if not

	Parallel(class PDPS *);
	virtual ~Parallel();

	void init();
	void setup();
	void exchange();
	void borders();
	virtual void set_proc_grid(int outflag = 1);
	void set_processors(int, char **);
	void cart_map(int, int *, int *, int [3][2], int ***);

	void forward_comm(int dummy = 0);   // forward comm of x and v
	void forward_force(int dummy = 0);   // forward comm of x and v
	void reverse_comm();                // reverse comm of forces
	void reverse_comm_pair(class Pair *);
	void forward_comm_pair(class Pair *);

protected:
	int style; 
	int map_style;
	int bordergroup;                    // only communicate this group in borders
	int maxneed[3];                     // max procs away any proc needs, per dim
	int maxswap;                        // max # of swaps memory is allocated for
	int nswap;                          // # of swaps to perform = sum of maxneed
	int recvneed[3][2];                 // # of procs away I recv atoms from (in one direction)
	int sendneed[3][2];                 // # of procs away I send atoms to (in one direction)
	
	int size_forward;					// # of per-atom datums in forward comm
	int size_reverse;					// # of datums in reverse comm
	int size_border;					// # of datums in forward border comm

	double *buf_send;				    // send buffer for all comm
	double *buf_recv;					// recv buffer for all comm
	int maxsend, maxrecv;				// current size of send/recv buffer
	int maxforward, maxreverse;			// max # of datums in forward/reverse comm

	int *firstrecv;                     // where to put 1st recv atom in each swap
	int **sendlist;                     // list of atoms to send in each swap
	int *maxsendlist;                   // max size of send list for each swap

	int *pbc_flag;						// general flag for sending atoms thru PBC
	int **pbc;							// dimension flags for PBC adjustments
	int *sendnum,*recvnum;				// # of atoms to send/recv in each swap
	int *sendproc,*recvproc;			// proc to send/recv to/from at each swap
	int *size_forward_recv;				// # of values to recv in each forward comm
	int *size_reverse_send;				// # to send in each reverse comm
	int *size_reverse_recv;				// # to recv in each reverse comm
	double *slablo,*slabhi;				// bounds of slab to send at each swap

	double **multilo, **multihi;        // bounds of slabs for multi-type swap
	double **rcghost_multi;             // cutghost on a per-type basis

	int comm_x_only, comm_f_only;       // 1 if only exchange x,f in for/rev comm
	
	void grow_send(int, int);           // reallocate send buffer
	void grow_recv(int);                // reallocate recv buffer
	void grow_list(int, int);           // reallocate one sendlist
	void grow_swap(int);                // grow swap and multi arrays
	void allocate_swap(int);            // allocate swap arrays
	void allocate_multi(int);           // allocate multi arrays
	void free_swap();                   // free swap arrays
	void free_multi();                  // free multi arrays

private:
	int find_factors3(int, int **);
	int find_factors2(int, int **);
	int best_factors(int, int **, int *, const int, const int, const int);
	int filter_user(int, int **, int, int *);
};

}

#endif
