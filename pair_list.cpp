/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h" 

#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair.h"
#include "pair_list.h"
#include "parallel.h"
#include "particle.h"

#include "update.h"

using namespace PDPS_NS;

#define BUFMIN 100
#define BUFEXTRA 100


/* ---------------------------------------------------------------------- */

PairList::PairList(PDPS *ps, int size) : Pointers(ps) 
{
	bucket_ntables = maxtable = 0;
	tbsize = size;
	last_table = last_index = 0;
	npairs = 0;
	nbuckets = 0;

	nbuilds = nbuilds_total = 0;

	hash_repeated = 0;
	bits_limit = 9;

	hash_pair = NULL;
	hash_bucket = NULL;

	// prime numbers
	nprimes = 38;
	primes = new int[nprimes];
	int plist[] = {5041,10007,20011,30011,40009,50021,60013,70001,80021,
                 90001,100003,110017,120011,130003,140009,150001,160001,
                 170003,180001,190027,200003,210011,220009,230003,240007,
                 250007,260003,270001,280001,290011,300007,310019,320009,
                 330017,340007,350003,362881,3628801};
	for (int i = 0; i < nprimes; i++) primes[i] = plist[i];

	// initialize prallel computing information
	for (int i = 0; i < 3; i++) {
		procgrid[i] = parallel->procgrid[i];
		for (int j = 0; j < 2; j++) {
			procneigh[i][j] = parallel->procneigh[i][j];	
		}
	}

	// undetermined pairs
	buf_str_send = buf_str_recv = NULL;
	max_unpairs_nsend = BUFMIN;
	memory->create(buf_str_send, max_unpairs_nsend, BYTES, "parallel: buf_str_send");
	max_unpairs_nrecv = BUFMIN;
	memory->create(buf_str_recv, max_unpairs_nrecv, BYTES, "parallel: buf_str_recv");

	// updated pairs
	up_unpairs_nsend = up_unpairs_nrecv = 0;
	up_nsend = up_nrecv = 0;
	buf_up_send = buf_up_recv = NULL;
	max_up_nsend = max_up_nrecv = 0;

	procid = parallel->procid;
	nprocs = parallel->nprocs;

	// ----------  debug purpose ------------
	f_log = f_pairlist = f_unpairlist = f_up_unpairlist = NULL;
	/*char str1[128], str2[128], str3[128];
	sprintf(str1, "p%d_pairlist.txt", procid);
	sprintf(str2, "p%d_unpairlist.txt", procid);
	sprintf(str3, "p%d_up_unpairlist.txt", procid);
	f_pairlist = fopen(str1,"w");
	f_unpairlist = fopen(str2,"w");
	f_up_unpairlist = fopen(str3,"w");*/
	// ---------------------------------------
}

/* ---------------------------------------------------------------------- */

PairList::~PairList()
{
	memory->destroy(buf_str_send);
	memory->destroy(buf_str_recv);
	
	if (buf_up_send) memory->destroy(buf_up_send);
	if (buf_up_recv) memory->destroy(buf_up_recv);

	delete[] primes;
	primes = NULL;

	if (bucket_ntables) {
		for (int i = 0; i < bucket_ntables; i++) {
			delete[] hash_bucket[i];
			hash_bucket[i] = NULL;
		}
		delete[] hash_bucket;
		hash_bucket = NULL;
	}

	if (maxtable) {
		for (int i = 0; i < maxtable; i++) {
			for (int j = 0; j < tbsize; j++) {
				delete[] hash_pair[i][j].key;
				hash_pair[i][j].key = NULL;
				if (ifields > 0) {
					delete[] hash_pair[i][j].ivalues;
					hash_pair[i][j].ivalues = NULL;
				}
				if (dfields > 0) {
					delete[] hash_pair[i][j].dvalues;
					hash_pair[i][j].dvalues = NULL;
				}
			}
			delete[] hash_pair[i];
			hash_pair[i] = NULL;
		}
		memory->destroy(hash_pair);
	}
}

/* ---------------------------------------------------------------------- 
         add tables
---------------------------------------------------------------------- */

void PairList::add_tables(int n)
{
	int i, j;
	int ifield;
	int ntables = maxtable;

	maxtable += n;

	hash_pair = (HashPair **) memory->srealloc(hash_pair, maxtable*sizeof(HashPair *), "PairList: hash_pair");

	for (i = ntables; i < maxtable; i++) {
		hash_pair[i] = NULL;
		hash_pair[i] = new HashPair[tbsize];
		// initialize hash_pair
		for (j = 0; j < tbsize; j++) {
			hash_pair[i][j].key = new char[BYTES];
			if (ifields > 0) {
				hash_pair[i][j].ivalues = new int[ifields];
			}
			if (dfields > 0) {
				hash_pair[i][j].dvalues = new double[dfields];
			}
		}
	}
}

/* ---------------------------------------------------------------------- 
         Initiate hash: allocation of hash_bucket and hash_pair
---------------------------------------------------------------------- */

void PairList::init_hash(int size)
{
	int i, j;

	int maxsize;

	nbuilds = 0;
	
	// # of elements stored for updating one undetermined pair:
	// index of undetermined pair + ifields = dfields
	up_size_one = 1 + ifields + dfields;

	if (buf_up_send == NULL) {
		max_up_nsend = BUFMIN*up_size_one;
		memory->create(buf_up_send, max_up_nsend, "PairList: buf_up_send");
	}
	if (buf_up_recv == NULL) {
		max_up_nrecv = BUFMIN*up_size_one;
		memory->create(buf_up_recv, max_up_nrecv, "PairList: buf_up_recv");
	}

	MPI_Allreduce(&size, &maxsize, 1, MPI_INT, MPI_MAX, mworld);

	int ntables = maxsize / tbsize;
	if (ntables*tbsize < maxsize) ntables++;

	ntables = ntables - maxtable;
	// at least create one page for each processor 
	// at the very beginning
	if (ntables > 0) {
		add_tables(ntables);
	}
	else if (ntables == 0 && maxtable == 0) {
		add_tables(1);
	}

	// first time to allocate hash_bucket
	if (nbuckets == 0) {
		nbuckets = 2 * size;

		if (nbuckets <= primes[0]) nbuckets = primes[0];
		else {
			for (i = 0; i < nprimes; i++) {
				if (nbuckets < primes[i]) break;
			}
			if (i < nprimes) nbuckets = primes[i];
			else nbuckets = primes[nprimes-1];
		}

		bucket_ntables = TABLE(nbuckets);
		if (bucket_ntables*tbsize < nbuckets) bucket_ntables++;
		hash_bucket = new int*[bucket_ntables];
		for (i = 0; i < bucket_ntables; i++) {
			hash_bucket[i] = NULL;
			hash_bucket[i] = new int[tbsize];
			// initialize hash_bucket
			for (j = 0; j < tbsize; j++) {
				hash_bucket[i][j] = -1;
			}
		}
	} // if (nbuckets == 0)
}

/* ---------------------------------------------------------------------- 
                    Insert or update hash_pair
---------------------------------------------------------------------- */

void PairList::set_hash()
{
	nbuilds++;
	nbuilds_total++;

	// initialize # of undetermined pairs to send and recv
	unpairs_nsend = unpairs_nrecv = 0;
	str_nsend = str_nrecv = 0;

	
	// only get pair info. from neighbor processors
	// when nprocs > 1 and not the first build
	// if (nbuilds_total > 1 && nprocs > 1) exchange_pair();
}

/* ---------------------------------------------------------------------- 
            Set one pair based on particles' global tags
   1. String is constructed by "itag-jtag"
   2. In the first build, always insert hash
   3. If ipair = -1, always insert hash, but one also needs to check
   neighbor processors to see if they possess this pair at last build
   4. If ipair != -1 & its flag is not equal to the last build index,
   it means this processor once has this pair but then lost it. One 
   also needs to check neighbor processors to see if any processor 
   owns it in the last build
   5. For all pairs need to be further checked, they are categorized 
   into "undetermined pairs". They are packed into buffer: buf_str_send.
   6. buf_str_send will be sent to neighbor processors to check, if any
   processor has a valid pair, an updated undetermined pairs will be sent
   back to update the hash_pair in this processor
---------------------------------------------------------------------- */

int PairList::set_pair(char *str)
{
	int ipair;
	int temp;
	int table, index;

	ipair = find_hash(str);
	
	if (ipair == -1) {
		if (nprocs > 1) {
			if (unpairs_nsend >= max_unpairs_nsend) grow_str_send(BUFMIN);
			// if two particles were in this processor, actually there is no need to pack it
			// however, currently we do not store particles' position on the last build
			pack_pair_str(str);           
		}
		ipair = insert_hash(str);
		if (ipair == -1) error->all(FLERR, "Failed to insert the hash");
	}
	else if (ipair > -1 && nprocs > 1) {
		table = TABLE(ipair);
		index = INDEX(ipair);
		if (hash_pair[table][index].flag < nbuilds_total-1) {
			if (unpairs_nsend >= max_unpairs_nsend) grow_str_send(BUFMIN);
			pack_pair_str(str);
		}
	}

	table = TABLE(ipair);
	index = INDEX(ipair);
	hash_pair[table][index].flag = nbuilds_total;   

	return ipair;
}

/* ---------------------------------------------------------------------- 
              (Insert hash and also assign initial values it)
   Mark this pair's validity in this processor of this build by letting
              hash_pair[table][index].flag = nbuilds_total 
---------------------------------------------------------------------- */

int PairList::insert_hash(char *str) 
{
	int ipair;
	int inserted_table, inserted_index;
	int ibucket, bucket_table, bucket_index;
	int previous;

	unsigned int num_key = key2int(str);

	ibucket = num_key % nbuckets;
	bucket_table = TABLE(ibucket);
	bucket_index = INDEX(ibucket);

	previous = -1;
	ipair = hash_bucket[bucket_table][bucket_index];
	int table, index;
	while (ipair > -1) {
		table = TABLE(ipair);
		index = INDEX(ipair);
		if (strcmp(hash_pair[table][index].key, str) == 0) break;
		previous = ipair;	
		ipair = hash_pair[table][index].next;
		hash_repeated++;   // debug purpose
	}
	if (ipair > -1) {
		error->all(FLERR, "It already exists, and there is no need to insert a new hash pair");
	}

	// take one entry from the next free item
    // if this entry is 1st in bucket, store it in hash_bucket[][] 
    // if not, do not update hash_bucket[][] but add this new pair 
	// by let previous pair point to it: hash_pair[][].next = previous

	// check to see if one needs add one more table to hash_pair
	if (npairs == maxtable*tbsize) add_tables(1);     
	ipair = npairs;
	inserted_table = TABLE(ipair);
	inserted_index = INDEX(ipair);

	// previous == -1: it is the first entry
	if (previous == -1) hash_bucket[bucket_table][bucket_index] = ipair;
	else {
		hash_pair[TABLE(previous)][INDEX(previous)].next = ipair;
	}

	// initialization
	int i;
	strcpy(hash_pair[inserted_table][inserted_index].key, str);
	hash_pair[inserted_table][inserted_index].flag = nbuilds_total;
	for (i = 0; i < ifields; i++) hash_pair[inserted_table][inserted_index].ivalues[i] = 0;
	for (i = 0; i < dfields; i++) hash_pair[inserted_table][inserted_index].dvalues[i] = 0.0;
	hash_pair[inserted_table][inserted_index].next = -1;

	npairs++;
	// update the last_table and last_index for hash_pair
	last_table = TABLE(npairs-1);
	last_index = INDEX(npairs-1);

	return ipair;
}

/* ---------------------------------------------------------------------- 
       (Lookup the pair by string and return its pair index)
         ipair = (table)*tbsize + index
         ipair = 1: either it exists in the current neighbor list or
		            was stored before;
		 ipair = -1: not exist at all
---------------------------------------------------------------------- */

int PairList::find_hash(char *str)
{
	int ipair;
	int previous = -1;
	int ibucket, bucket_table, bucket_index;

	size_t num_key = key2int(str);

	ibucket = num_key % nbuckets;
	bucket_table = TABLE(ibucket);
	bucket_index = INDEX(ibucket);

	ipair = hash_bucket[bucket_table][bucket_index];

	int table, index;
	while (ipair > -1) {
		table = TABLE(ipair);
		index = INDEX(ipair);
		if (strcmp(hash_pair[table][index].key, str) == 0) {
			break;
		}
		ipair = hash_pair[table][index].next;
	}

	return ipair;
}

void PairList::set_zero(int ipair) 
{
	int table, index;
	table = TABLE(ipair);
	index = INDEX(ipair);
	hash_pair[table][index].ivalues[0] = 0;
	hash_pair[table][index].dvalues[0] = 0.0;
	hash_pair[table][index].dvalues[1] = 0.0;
	hash_pair[table][index].dvalues[2] = 0.0;
	
}

/* ---------------------------------------------------------------------- 
                 (Convert key string to integer)
	             FNV-1a algorithm is adopted here
	      http://www.isthe.com/chongo/tech/comp/fnv/#FNV-1a
---------------------------------------------------------------------- */
static const size_t offset_basis = 2166136261U;     // for 32 bit
static const size_t fnv1a_prime = 16777619;         // 2^24 + 2^8 + 0x93    

unsigned int PairList::key2int(char *str) 
{
	int n = strlen(str);

	size_t hash = offset_basis;

	//unsigned int hash = 0;
	for (int i = 0; i < n; i++) {
		//hash = 101 * hash + str[i];
		hash ^= offset_basis*str[i];
		hash = hash * fnv1a_prime;
	}

	return hash;
}

/* ---------------------------------------------------------------------- 
        (Exchange and update pair info from neighbor processors)
   1. The MPI computing is followed a similar strategy in Parallel class
   2. If nprocs == 1, no exchange; if nprocs == 2, only send undetermined
   pairs (buf_str_send) to west; if nprocs > 2, also send buf_str_send to 
   east
   3. buf_str_recv could be received from west or east processor
   4. buf_up_send: based on buf_str_recv, send the corresponding pairs
   info to that processor via buf_up_send
   5. buf_up_recv: updated info. for the undetermined pairs (buf_str_send)
   in this processor. Copy their ivalues and dvalues to the hash_pair[][]
---------------------------------------------------------------------- */

void PairList::exchange_pair() 
{
	int i, j, m, nlocal;
	double lo, hi, value;

	MPI_Request request;
	MPI_Status status;

	str_nsend = unpairs_nsend * BYTES;

	// loop over dimensions
	for (int dim = 0; dim < 3; dim++) {
		// send/recv undertermined pairs in both directions
		// if 1 proc in dimension, no send/recv, set recv buf to send buf
		// if 2 procs in dimension, single send/recv
		// if more than 2 procs in dimension, send/recv to both neighbors
		
		// search for the west direction 
		if (procgrid[dim] > 1) {
			// send and receive string of undetermined pairs
			MPI_Sendrecv(&str_nsend,1,MPI_INT,procneigh[dim][0],0,
				         &str_nrecv,1,MPI_INT,procneigh[dim][1],0,mworld,&status);
			
			if (str_nrecv % BYTES != 0) {
				error->all(FLERR, "Total # of str characters received \
								  is not a multiple of BYTES, bug may exist");
			}
			unpairs_nrecv = str_nrecv/BYTES;
			if (unpairs_nrecv > max_unpairs_nrecv) {
				grow_str_recv(unpairs_nrecv-max_unpairs_nrecv);
			}
			
			MPI_Irecv(&buf_str_recv[0][0],str_nrecv,MPI_CHAR,procneigh[dim][1],0,
                mworld,&request);
			MPI_Send(&buf_str_send[0][0],str_nsend,MPI_CHAR,procneigh[dim][0],0,mworld);
			MPI_Wait(&request,&status);

			// parse the string of undetermined pairs and store its values if it exists
			parse_str_recv();

			// send and receive updated pair info
			MPI_Sendrecv(&up_nsend,1,MPI_INT,procneigh[dim][1],0,
				         &up_nrecv,1,MPI_INT,procneigh[dim][0],0,mworld,&status);
			
			if (up_nrecv % up_size_one != 0) {
				error->all(FLERR, "Total # of packed items is \
								  not multiple of up_size_one, bug may exist");
			}
			up_unpairs_nrecv = up_nrecv / up_size_one;	
			if (up_nrecv > max_up_nrecv) grow_up_recv(up_nrecv-max_up_nrecv);

			MPI_Irecv(&buf_up_recv[0], up_nrecv, MPI_DOUBLE, procneigh[dim][0], 0, 
                mworld, &request);
			MPI_Send(&buf_up_send[0], up_nsend, MPI_DOUBLE, procneigh[dim][1], 0, mworld);
			MPI_Wait(&request, &status);

			// copy updated pair info
			unpack_pair_values();

			// search for the east direction, if necessary
			if (procgrid[dim] > 2) {
				// send and receive string of undetermined pairs
				MPI_Sendrecv(&str_nsend,1,MPI_INT,procneigh[dim][1],0,
							 &str_nrecv,1,MPI_INT,procneigh[dim][0],0,mworld,&status);
				
				if (str_nrecv % BYTES != 0) {
					error->all(FLERR, "Total # of str characters received \
									  is not a multiple of BYTES, bug may exist");
				}
				unpairs_nrecv = str_nrecv/BYTES;
				if (unpairs_nrecv > max_unpairs_nrecv) {
					grow_str_recv((unpairs_nrecv-max_unpairs_nrecv)*BYTES);
				}
				
				MPI_Irecv(&buf_str_recv[0][0],str_nrecv,MPI_CHAR,procneigh[dim][0],0,
					mworld,&request);
				MPI_Send(&buf_str_send[0][0],str_nsend,MPI_CHAR,procneigh[dim][1],0,mworld);
				MPI_Wait(&request,&status);

				// parse the string of undetermined pairs and store its values if it exists
				parse_str_recv();
		
				// send and receive updated pair info
				MPI_Sendrecv(&up_nsend,1,MPI_INT,procneigh[dim][0],0,
							 &up_nrecv,1,MPI_INT,procneigh[dim][1],0,mworld,&status);
				
				if (up_nrecv % up_size_one != 0) {
					error->all(FLERR, "Total # of packed items is \
									  not multiple of up_size_one, bug may exist");
				}
				up_unpairs_nrecv = up_nrecv / up_size_one;
				if (up_nrecv > max_up_nrecv) grow_up_recv(up_nrecv-max_up_nrecv);
				MPI_Irecv(&buf_up_recv[0], up_nrecv, MPI_DOUBLE, procneigh[dim][1], 0, 
					mworld, &request);
				MPI_Send(&buf_up_send[0], up_nsend, MPI_DOUBLE, procneigh[dim][0], 0, mworld);
				MPI_Wait(&request, &status);

				// copy updated pair info
				unpack_pair_values();
			} // if (procgrid[dim] > 2) 
		} // if (procgrid[dim] > 1)
	} // for (int dim = 0; dim < 3; dim++)
}

/* ---------------------------------------------------------------------- 
                pack pair
---------------------------------------------------------------------- */

void PairList::pack_pair_str(char *str) 
{
	strcpy(buf_str_send[unpairs_nsend], str);

	unpairs_nsend++;
}

/* ---------------------------------------------------------------------- 
   Parse buf_str_recv, if it has the corresponding pair, store it in the 
   buf_up_send and send to the corresponding processor
---------------------------------------------------------------------- */

void PairList::parse_str_recv() 
{
	int i;
	int ipair;
	int table, index;

	// initialize # of updated undetermined pairs to send and recv
	up_unpairs_nsend = up_unpairs_nrecv = 0;
	up_nsend = up_nrecv = 0;
	for (i = 0; i < unpairs_nrecv; i++) {
		if (buf_str_recv[i] == NULL) error->all(FLERR, "buf_str_recv[i] == NULL");
		ipair = find_hash(buf_str_recv[i]);
		if (ipair == -1) continue;
		else {
			table = TABLE(ipair);
			index = INDEX(ipair);
			if (hash_pair[table][index].flag == nbuilds_total) {
				char str[128];
				sprintf(str, "Pair %s also exist in this processor", buf_str_recv[i]);
				error->all(FLERR, str);
			}
			else if (hash_pair[table][index].flag == nbuilds_total-1) {
				if (up_nsend >= max_up_nsend) grow_up_send(BUFMIN*up_size_one);
				buf_up_send[up_nsend++] = i;
				up_nsend += pack_pair_values(ipair,&buf_up_send[up_nsend]);
				up_unpairs_nsend++;
			}
		}
	}
}

/* ---------------------------------------------------------------------- 
   Pack pair's index in the undetermined pairs, ivalues and dvalues into
   buf_str_send
---------------------------------------------------------------------- */

int PairList::pack_pair_values(int ipair, double *buf) 
{
	int table, index;
	
	table = TABLE(ipair);
	index = INDEX(ipair);

	int m = 0;
	for (int i = 0; i < ifields; i++) {
		buf[m++] = hash_pair[table][index].ivalues[i];
	}
	for (int i = 0; i < dfields; i++) {
		buf[m++] = hash_pair[table][index].dvalues[i];
	}

	return m;
}

/* ---------------------------------------------------------------------- 
   Unpack pair from buf_up_recv and copy its values to hash_pair[][]
---------------------------------------------------------------------- */

void PairList::unpack_pair_values() 
{
	int un_index, ipair;
	int table, index;
	int m = 0;
	for (int i = 0; i < up_unpairs_nrecv; i++) {
		un_index = static_cast <int> (buf_up_recv[m++]);
		if (buf_str_send[un_index] == NULL) {
			error->all(FLERR, "Sth is wrong for buf_str_send");
		}
		ipair = find_hash(buf_str_send[un_index]);
		if (ipair == -1) error->all(FLERR, "Cannot find the pair in the buf_str_send");
		table = TABLE(ipair);
		index = INDEX(ipair);
		for (int j = 0; j < ifields; j++) {
			hash_pair[table][index].ivalues[j] = static_cast <int> (buf_up_recv[m++]);
		}
		for (int j = 0; j < dfields; j++) {
			hash_pair[table][index].dvalues[j] = buf_up_recv[m++];
		}
	}
}

/* ---------------------------------------------------------------------- */

void PairList::grow_str_send(int n) 
{
	max_unpairs_nsend += n;

	memory->grow(buf_str_send, max_unpairs_nsend, BYTES, "parallel: buf_str_send");
}

/* ---------------------------------------------------------------------- */

void PairList::grow_str_recv(int n) 
{
	max_unpairs_nrecv += n;

	memory->grow(buf_str_recv, max_unpairs_nrecv, BYTES, "parallel: buf_str_recv");
}

/* ---------------------------------------------------------------------- */

void PairList::grow_up_send(int n) 
{
	max_up_nsend += n;

	memory->grow(buf_up_send, max_up_nsend, "parallel: buf_str_send");
}

/* ---------------------------------------------------------------------- */

void PairList::grow_up_recv(int n) 
{
	max_up_nrecv += n;

	memory->grow(buf_up_recv, max_up_nrecv, "parallel: buf_str_recv");
}

/* ---------------------------------------------------------------------- */

void PairList::print_pairlist() 
{
	int table, index;

	fprintf(f_pairlist, "nbuilds_total = %d ntimestep = %d\n", nbuilds_total, update->ntimestep);
	for (int i = 0; i < npairs; i++) {
		table = TABLE(i);
		index = INDEX(i);
		fprintf(f_pairlist, "%s ", hash_pair[table][index].key);
		fprintf(f_pairlist, "%d ", hash_pair[table][index].flag);
		for (int j = 0; j < ifields; j++) {
			fprintf(f_pairlist, "%d ", hash_pair[table][index].ivalues[j]);
		}
		for (int j = 0; j < dfields; j++) {
			fprintf(f_pairlist, "%f ", hash_pair[table][index].dvalues[j]);
		}
		fprintf(f_pairlist, "%d ", hash_pair[table][index].next);
		fprintf(f_pairlist, "\n");
	}
	fprintf(f_pairlist, "\n");
}

/* ---------------------------------------------------------------------- */

void PairList::print_unpairlist_send() 
{
	fprintf(f_unpairlist, "nbuilds_total = %d \nbuf_str_send: unpairs_nsend \
						  = %d\n", nbuilds_total, unpairs_nsend);
	for (int i = 0; i < unpairs_nsend; i++) {
		fprintf(f_unpairlist, "%s\n", buf_str_send[i]);
	}
}

/* ---------------------------------------------------------------------- */

void PairList::print_unpairlist_recv(int neigh_proc) 
{
	fprintf(f_unpairlist, "nbuilds_total = %d \nbuf_str_recv: unpairs_nrecv = \
						  %d from processor %d\n", nbuilds_total, unpairs_nrecv, neigh_proc);
	for (int i = 0; i < unpairs_nrecv; i++) {
		fprintf(f_unpairlist, "%s\n", buf_str_recv[i]);
	}
	fprintf(f_unpairlist, "\n");
}

/* ---------------------------------------------------------------------- */

void PairList::print_up_unpairlist(int neigh_proc_send, int neigh_proc_recv) 
{
	int m;

	fprintf(f_up_unpairlist, "nbuilds_total = %d \nbuf_up_send: up_unpairs_nsend = %d to processor %d\n", 
		nbuilds_total, up_unpairs_nsend, neigh_proc_send);
	m = 0;
	for (int i = 0; i < up_unpairs_nsend; i++) {
		fprintf(f_up_unpairlist, "%s: ", buf_str_recv[static_cast<int> (buf_up_send[m])]);
		fprintf(f_up_unpairlist, "%d ", static_cast<int> (buf_up_send[m++]));
		for (int j = 0; j < ifields; j++) {
			fprintf(f_up_unpairlist, "%d ", static_cast<int> (buf_up_send[m++]));
		}
		for (int j = 0; j < dfields; j++) {
			fprintf(f_up_unpairlist, "%f ", buf_up_send[m++]);
		}
		fprintf(f_up_unpairlist, "\n");
	}

	fprintf(f_up_unpairlist, "nbuilds_total = %d \nbuf_up_recv: up_unpairs_nrecv = %d from processor %d\n", 
		nbuilds_total, up_unpairs_nrecv, neigh_proc_recv);
	m = 0;
	for (int i = 0; i < up_unpairs_nrecv; i++) {
		fprintf(f_up_unpairlist, "%s: ", buf_str_send[static_cast<int> (buf_up_recv[m])]);
		fprintf(f_up_unpairlist, "%d ", static_cast<int> (buf_up_recv[m++]));
		for (int j = 0; j < ifields; j++) {
			fprintf(f_up_unpairlist, "%d ", static_cast<int> (buf_up_recv[m++]));
		}
		for (int j = 0; j < dfields; j++) {
			fprintf(f_up_unpairlist, "%f ", buf_up_recv[m++]);
		}
		fprintf(f_up_unpairlist, "\n");
	}
	fprintf(f_up_unpairlist, "\n");
}
