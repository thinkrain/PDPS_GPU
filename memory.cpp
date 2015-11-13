/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
//#include "stdint.h"

#include "error.h"
#include "memory.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

Memory::Memory(PDPS *ps) : Pointers(ps) {}

/* ----------------------------------------------------------------------
						Safe Malloc
------------------------------------------------------------------------- */

void *Memory::smalloc(bigint nbytes, const char *name)
{
  if (nbytes == 0) return NULL;
  
  void *ptr = malloc(nbytes);

  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to allocate %d bytes for array %s",
            nbytes,name);
    error->all(FLERR,str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe realloc
------------------------------------------------------------------------- */

void *Memory::srealloc(void *ptr, bigint nbytes, const char *name)
{
	if (nbytes == 0) {
		//destroy(ptr);
		free(ptr);
		return NULL;
	}

	ptr = realloc(ptr,nbytes);
	if (ptr == NULL) {
		char str[128];
		sprintf(str,"Failed to reallocate %d bytes for array %s",
            nbytes,name);
		error->all(FLERR,str);
	}
	return ptr;
}

/* ----------------------------------------------------------------------
   safe free
------------------------------------------------------------------------- */

void Memory::sfree(void *ptr)
{
	if (ptr == NULL) return;
	free(ptr);
	ptr = NULL;      //  avoid random pointer
}

/* ----------------------------------------------------------------------
   erroneous usage of templated create/grow functions
------------------------------------------------------------------------- */

void Memory::fail(const char *name)
{
  char str[128];
  sprintf(str,"Cannot create/grow a vector/array of pointers for %s",name);
  error->all(FLERR,str);
}
