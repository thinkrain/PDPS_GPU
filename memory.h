/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_MEMORY_H
#define PS_MEMORY_H

#include "pointers.h"

namespace PDPS_NS {

class Memory : protected Pointers {
public:
	Memory(class PDPS *);
	void *smalloc(bigint n, const char *);
	void *srealloc(void *, bigint n, const char *);
	void sfree(void *);
	void fail(const char *);

/* ----------------------------------------------------------------------
   Create a 1d array
------------------------------------------------------------------------- */

template <typename TYPE>
	TYPE *create(TYPE *&array, int n, const char *name)
	{
		bigint nbytes = ((bigint) sizeof(TYPE)) * n;
		array = (TYPE *) smalloc(nbytes,name);
		return array;
	}

template <typename TYPE>
	TYPE **create(TYPE **&array, int n, const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   Grow or shrink 1d array
---------------------------------------------------------------------- */

template <typename TYPE>
	TYPE *grow(TYPE *&array, int n, const char *name)
	{
		if (array == NULL) return create(array,n,name);

		bigint nbytes = ((bigint) sizeof(TYPE)) * n;
		array = (TYPE *) srealloc(array,nbytes,name);
		return array;
	}

template <typename TYPE>
	TYPE **grow(TYPE **&array, int n, const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   Destroy a 1d array
------------------------------------------------------------------------- */

template <typename TYPE>
	void destroy(TYPE *array)
	{
		sfree(array);
	}

/* ----------------------------------------------------------------------
   Create a 2d array
------------------------------------------------------------------------- */

template <typename TYPE>
	TYPE **create(TYPE **&array, int n1, int n2, const char *name)
	{
		bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2;
		TYPE *data = (TYPE *) smalloc(nbytes,name);
		nbytes = ((bigint) sizeof(TYPE *)) * n1;
		array = (TYPE **) smalloc(nbytes,name);

		bigint n = 0;
		for (int i = 0; i < n1; i++) {
		array[i] = &data[n];
		n += n2;
		}
		return array;
	}

template <typename TYPE>
	TYPE ***create(TYPE ***&array, int n1, int n2, const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   Grow or shrink 1st dim of a 2d array
   last dim must stay the same
------------------------------------------------------------------------- */

template <typename TYPE>
	TYPE **grow(TYPE **&array, int n1, int n2, const char *name)
	{
		if (array == NULL) return create(array,n1,n2,name);

		bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2;
		TYPE *data = (TYPE *) srealloc(array[0],nbytes,name);
		nbytes = ((bigint) sizeof(TYPE *)) * n1;
		array = (TYPE **) srealloc(array,nbytes,name);

		bigint n = 0;
		for (int i = 0; i < n1; i++) {
			array[i] = &data[n];
			n += n2;
		}
		return array;
	}

template <typename TYPE>
	TYPE ***grow(TYPE ***&array, int n1, int n2, const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   Destroy a 2d array
------------------------------------------------------------------------- */

template <typename TYPE>
	void destroy(TYPE **array)
	{
		if (array == NULL) return;
		sfree(array[0]);
		sfree(array);
	}

/* ----------------------------------------------------------------------
   create a 3d array
------------------------------------------------------------------------- */

template <typename TYPE>
    TYPE ***create(TYPE ***&array, int n1, int n2, int n3, const char *name)
    {
		bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3;
		TYPE *data = (TYPE *) smalloc(nbytes,name);
		nbytes = ((bigint) sizeof(TYPE *)) * n1*n2;
		TYPE **plane = (TYPE **) smalloc(nbytes,name);
		nbytes = ((bigint) sizeof(TYPE **)) * n1;
		array = (TYPE ***) smalloc(nbytes,name);

		int i,j;
		bigint m;
		bigint n = 0;
		for (i = 0; i < n1; i++) {
			m = ((bigint) i) * n2;
			array[i] = &plane[m];
			for (j = 0; j < n2; j++) {
				plane[m+j] = &data[n];
				n += n3;
			}
		}
		return array;
    }

template <typename TYPE>
    TYPE ****create(TYPE ****&array, int n1, int n2, int n3, const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   Grow or shrink 1st dim of a 3d array
   last 2 dims must stay the same
------------------------------------------------------------------------- */

template <typename TYPE>
	TYPE ***grow(TYPE ***&array, int n1, int n2, int n3, const char *name)
	{
		if (array == NULL) return create(array,n1,n2,n3,name);

		bigint nbytes = ((bigint) sizeof(TYPE)) * n1*n2*n3;
		TYPE *data = (TYPE *) srealloc(array[0][0],nbytes,name);
		nbytes = ((bigint) sizeof(TYPE *)) * n1*n2;
		TYPE **plane = (TYPE **) srealloc(array[0],nbytes,name);
		nbytes = ((bigint) sizeof(TYPE **)) * n1;
		array = (TYPE ***) srealloc(array,nbytes,name);

		int i,j;
		bigint m;
		bigint n = 0;
		for (i = 0; i < n1; i++) {
			m = ((bigint) i) * n2;
			array[i] = &plane[m];
			for (j = 0; j < n2; j++) {
				plane[m+j] = &data[n];
				n += n3;
			}
		}
		return array;
	}

template <typename TYPE>
    TYPE ****grow(TYPE ****&array, int n1, int n2, int n3, const char *name) {fail(name);}

/* ----------------------------------------------------------------------
   destroy a 3d array
------------------------------------------------------------------------- */

template <typename TYPE>
	void destroy(TYPE ***array)
	{
		if (array == NULL) return;
		sfree(array[0][0]);
		sfree(array[0]);
		sfree(array);
	}

/* ----------------------------------------------------------------------
   Memory usage of arrays, including pointers
------------------------------------------------------------------------- */

template <typename TYPE>
	bigint usage(TYPE *array, int n)
	{
		bigint bytes = ((bigint) sizeof(TYPE)) * n;
		return bytes;
	}

template <typename TYPE>
	bigint usage(TYPE **array, int n1, int n2)
	{
		bigint bytes = ((bigint) sizeof(TYPE)) * n1*n2;
		bytes += ((bigint) sizeof(TYPE *)) * n1;
		return bytes;
	}

template <typename TYPE>
	bigint usage(TYPE ***array, int n1, int n2, int n3)
	{
		bigint bytes = ((bigint) sizeof(TYPE)) * n1*n2*n3;
		bytes += ((bigint) sizeof(TYPE *)) * n1*n2;
		bytes += ((bigint) sizeof(TYPE **)) * n1;
		return bytes;
	}

};
}
#endif
