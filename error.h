/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_ERROR_H
#define PS_ERROR_H

#include "pointers.h"

namespace PDPS_NS {

class Error : protected Pointers {
 public:
     Error(class PDPS *);

     void all(const char *, int, const char *);
     void one(const char *, int, const char *);
     void warning(const char *, int, const char *);
     void message(const char *, int, const char *, int = 1);
     void done();
};

}

#endif
