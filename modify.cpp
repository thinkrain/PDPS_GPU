/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"

#include "analyze.h"
#include "compute.h"
#include "fix.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "style_analyze.h"
#include "style_compute.h"
#include "style_fix.h"

using namespace PDPS_NS;
using namespace FixConst;

#define DELTA 4

/* ---------------------------------------------------------------------- */

Modify::Modify(PDPS *ps) : Pointers(ps)
{
	analyze = NULL;
	compute = NULL;
	fix = NULL;
	fmask = NULL;
	list_initial_integrate = list_final_integrate = NULL;
	list_pre_force = list_post_force = NULL;
	list_pre_integrate = list_post_integrate = NULL;
	list_end_of_step = NULL;

	nanalyzes = maxanalyze = 0;
	ncomputes = maxcompute = 0;
	nfixes = maxfix = 0;
	
	n_initial_integrate = n_final_integrate = 0;
	n_pre_force = n_post_force = 0;
	n_pre_integrate = n_post_integrate = 0;
	n_end_of_step = 0;
}

/* ---------------------------------------------------------------------- */

Modify::~Modify()
{
	int i;

	for (i = 0; i < nanalyzes; i++) {
		delete analyze[i];
		analyze[i] = NULL;
	}
	memory->sfree(analyze);
	analyze = NULL;

	for (i = 0; i < ncomputes; i++) {
		delete compute[i];
		compute[i] = NULL;
	}
	memory->sfree(compute);
	compute = NULL;

	for (i = 0; i < nfixes; i++) { 
		delete fix[i];
		fix[i] = NULL;
	}
	memory->sfree(fix);
	fix = NULL;
	
	memory->sfree(fmask);
	fmask = NULL;
}

/* ----------------------------------------------------------------------
   Add a New Analyze
------------------------------------------------------------------------- */

void Modify::add_analyze(int narg, char **arg) 
{
	if (nanalyzes == maxanalyze) {
		maxanalyze += DELTA;
		analyze = (Analyze **)
			memory->srealloc(analyze, maxanalyze*sizeof(Analyze *), "Modify: analyze");
	}

	if (0) return;

#define ANALYZE_CLASS
#define AnalyzeStyle(key,Class) \
    else if (strcmp(arg[2],#key) == 0) \
		analyze[nanalyzes] = new Class(ps,narg,arg);
#include "style_analyze.h"
#undef AnalyzeStyle
#undef ANALYZE_CLASS

	else error->all(FLERR, "Invalid analyze style");

	nanalyzes++;
}

/* ----------------------------------------------------------------------
   Delete an Analyze
------------------------------------------------------------------------- */

void Modify::delete_analyze(const char *name) 
{
	int ianalyze = find_analyze(name);
	if (ianalyze < 0) error->all(FLERR, "Could not find analyze ID to delete");
	delete analyze[ianalyze];
	
	// move other Analyzes down in list one slot

	for (int i = ianalyze+1; i < nanalyzes; i++) {
		analyze[i-1] = analyze[i];
	}
	nanalyzes--;
}

/* ----------------------------------------------------------------------
   Find an Analyze
------------------------------------------------------------------------- */

int Modify::find_analyze(const char *name) 
{
	int aid;
	aid = -1;
	for (int i = 0; i < nanalyzes; i++) {
		if (strcmp(name, analyze[i]->name) == 0) 
			aid = i;
	}
	return aid;
}

/* ----------------------------------------------------------------------
   Add a New Compute
------------------------------------------------------------------------- */

void Modify::add_compute(int narg, char **arg) 
{
	// if (narg < 3) error->...

	// Expand compute list
	if (ncomputes == maxcompute) {
		maxcompute += DELTA;
		compute = (Compute **) 
			memory->srealloc(compute,maxcompute*sizeof(Compute *),"Modify: compute");
	}

	if (0) return;

#define COMPUTE_CLASS
#define ComputeStyle(key,Class) \
	else if (strcmp(arg[2],#key) == 0) \
		compute[ncomputes] = new Class(ps,narg,arg);
#include "style_compute.h"
#undef ComputeStyle
#undef COMPUTE_CLASS

	else error->all(FLERR, "Invalid compute style");

	ncomputes++;
}

/* ----------------------------------------------------------------------
   Add a New Fix
------------------------------------------------------------------------- */

void Modify::add_fix(int narg, char **arg) 
{
	if (narg < 3) error->all(FLERR,"Illegal fix command");

	// Expand compute list
	if (nfixes == maxfix) {
		maxfix += DELTA;
		fix = (Fix **) 
			memory->srealloc(fix,maxfix*sizeof(Compute *),"Modify: compute");
		fmask = (int *) memory->srealloc(fmask,maxfix*sizeof(int),"Modify: fmask");
	}

	if (0) return;

#define FIX_CLASS
#define FixStyle(key,Class) \
    else if (strcmp(arg[2],#key) == 0) \
		fix[nfixes] = new Class(ps,narg,arg);
#include "style_fix.h"
#undef FixStyle
#undef FIX_CLASS

	else error->all(FLERR, "Illegal fix style");

	// set fix mask  (check FixConst in fix.h)
	fmask[nfixes] = fix[nfixes]->setmask();
	nfixes++;
}

/* ----------------------------------------------------------------------
   Delete a fix
------------------------------------------------------------------------- */

void Modify::delete_fix(const char *name)
{
	int ifix = find_fix(name);
	if (ifix < 0) error->all(FLERR, "Could not find fix ID to delete");
	delete fix[ifix];

	// move other fixes down in list one slot

	for (int i = ifix + 1; i < nfixes; i++) {
		fix[i - 1] = fix[i];
	}
	nfixes--;
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_integrate()
{
	for (int i = 0; i < n_pre_integrate; i++)
		fix[list_pre_integrate[i]]->pre_integrate();
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::initial_integrate()
{
  for (int i = 0; i < n_initial_integrate; i++)
    fix[list_initial_integrate[i]]->initial_integrate();
}

/* ----------------------------------------------------------------------
   2nd half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::final_integrate()
{
  for (int i = 0; i < n_final_integrate; i++)
    fix[list_final_integrate[i]]->final_integrate();
}

void Modify::post_integrate()
{
	for (int i = 0; i < n_post_integrate; i++) 
		fix[list_post_integrate[i]]->post_integrate();
}

/* ----------------------------------------------------------------------
   Find compute id
------------------------------------------------------------------------- */

int Modify::find_compute(const char *cname) 
{
	int cid;
	cid = -1;
	for (int i = 0; i < ncomputes; i++) {
		if (strcmp(cname,compute[i]->name) == 0) 
			cid = i;
	}
	return cid;
}

/* ----------------------------------------------------------------------
   Find compute id based on style
------------------------------------------------------------------------- */

int Modify::find_compute_style(const char *cstyle) 
{
	int cid;
	cid = -1;
	for (int i = 0; i < ncomputes; i++) {
		if (strcmp(cstyle,compute[i]->style) == 0) 
			cid = i;
	}
	return cid;
}

/* ----------------------------------------------------------------------
   Find fix id
------------------------------------------------------------------------- */

int Modify::find_fix(const char *fname) 
{
	int fid;
	fid = -1;
	for (int i = 0; i < nfixes; i++) {
		if (!strcmp(fname,fix[i]->name) == 0) 
			fid = i;
	}

	return fid;
}

/* ----------------------------------------------------------------------
   initialize all fixes and computes
------------------------------------------------------------------------- */

void Modify::init()
{
	int i;

	// create lists of fixes to call at each stage of run

	list_init(PRE_INTEGRATE,n_pre_integrate,list_pre_integrate);
	list_init(INITIAL_INTEGRATE,n_initial_integrate,list_initial_integrate);
	list_init(FINAL_INTEGRATE,n_final_integrate,list_final_integrate);
	list_init(PRE_FORCE, n_pre_force, list_pre_force);
	list_init(POST_FORCE,n_post_force,list_post_force);
	list_init(POST_INTEGRATE,n_post_integrate,list_post_integrate);
	list_init(END_OF_STEP,n_end_of_step,list_end_of_step);

	for (i = 0; i < nanalyzes; i++) analyze[i]->init();
	
	for (i = 0; i < ncomputes; i++) compute[i]->init();

	for (i = 0; i < nfixes; i++) fix[i]->init();
}

/* ----------------------------------------------------------------------
   initialize all fixes and computes
------------------------------------------------------------------------- */

void Modify::setup()
{
	int i;

	for (i = 0; i < nanalyzes; i++) analyze[i]->setup();

	for (int i = 0; i < nfixes; i++) fix[i]->setup();
}

/* ----------------------------------------------------------------------
   create list of fix indices for fixes which match mask
------------------------------------------------------------------------- */

void Modify::list_init(int mask, int &n, int *&list)
{
	delete [] list;

	n = 0;
	for (int i = 0; i < nfixes; i++) if (fmask[i] & mask) n++;
	list = new int[n];

	n = 0;
	for (int i = 0; i < nfixes; i++) if (fmask[i] & mask) list[n++] = i;
}

/* ----------------------------------------------------------------------
post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_force()
{
	for (int i = 0; i < n_pre_force; i++) {
		fix[list_pre_force[i]]->pre_force();
	}
}

/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_force()
{
	for (int i = 0; i < n_post_force; i++) {
		fix[list_post_force[i]]->post_force();
	}
}


/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::end_of_step()
{
	for (int i = 0; i < n_end_of_step; i++) {
		fix[list_end_of_step[i]]->end_of_step();
	}
}
/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::check_analyze()
{
	for (int i = 0; i < nanalyzes; i++) {
		analyze[i]->invoke_analyze();
	}
}
