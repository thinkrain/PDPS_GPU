/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_Pair_List_H
#define PS_Pair_List_H

#include "pointers.h"

#define BYTES 128
#define TABLE(A) (A / tbsize)
#define INDEX(A) (A % tbsize)

namespace PDPS_NS {

class PairList : protected Pointers {
public:
	int pair_id; 

	struct HashPair {
		char *key;
		int flag;
		int *ivalues;
		double *dvalues;
		int next;
	};

	// for HashPair
	int maxtable;                // max # of tables created
	int tbsize;                  // table size for each table
	int last_table;
	int last_index;              // last valid index
	int npairs;                  // # of valid pairs stored

	int nbuckets;
	int bucket_ntables;
	int *primes;
	int nprimes;
	HashPair **hash_pair;
	int **hash_bucket;

	int ifields;
	int dfields;

	int nbuilds;
	int nbuilds_total;

	int unpairs_nsend, unpairs_nrecv;           // # of undetermined pairs to send and recv
	
	int str_nsend, str_nrecv;                   // total str characters to send and recv

	int up_unpairs_nsend, up_unpairs_nrecv;     // # of updated undetermined pairs to send and recv
	int up_nsend, up_nrecv;                     // total elements to send and recv

	PairList(class PDPS *, int);
	~PairList();

	void init_hash(int);
	int insert_hash(char *);
	void set_hash();
	
	int find_hash(char *);              // find hash by the key string
	void set_zero(int);

	int set_pair(char *);
	void exchange_pair();

	int hash_repeated;
	
private:
	// processors information
	int procid, nprocs;
	int procneigh[3][2];
	int procgrid[3];

	// buf for sending and receving the strings of undetermined pairs
	char **buf_str_send;				        // send buffer for undetermined pairs
	char **buf_str_recv;			            // recv buffer for undetermined pairs from neighbor processors
	int max_unpairs_nsend, max_unpairs_nrecv;	// max. # of undetermined pairs to send and recv
	
	// buf for sending and receiving the updated pair info. 
	// for the undetermined pairs 
	double *buf_up_send;                        // send buffer for updating undetermined pairs of neighbor processors
	double *buf_up_recv;                        // recv buffer for updating undetermined pairs
	
	int max_up_nsend, max_up_nrecv;	            // max # of updated pairs to send and recv
	int up_size_one;

	int bits_limit;                             // used in "int key2int(char *)" function
	                                            // max # of bits for the integer converted from key

	// for hash_pair
	void add_tables(int);

	void clear_hash();

	

	unsigned int key2int(char *);

	// parallel communication for updating undetermined pairs
	
	void pack_pair_str(char *);
	void parse_str_recv();
	int pack_pair_values(int, double *);
	void unpack_pair_values();

	void grow_str_send(int);
	void grow_str_recv(int);
	void grow_up_send(int);
	void grow_up_recv(int);

	//--------- debug purpose ----------
	FILE *f_log;
	FILE *f_pairlist;
	FILE *f_unpairlist;
	FILE *f_up_unpairlist;

	void print_pairlist();
	void print_unpairlist_send();
	void print_unpairlist_recv(int);
	void print_up_unpairlist(int, int);
	// ----------------------------------
};

}

#endif
