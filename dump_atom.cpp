/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the PDPS directory.
------------------------------------------------------------------------- */

#include "string.h"

#include "domain.h"
#include "dump_atom.h"
#include "error.h"
#include "group.h"
#include "particle.h"
#include "update.h"

using namespace PDPS_NS;

/* ---------------------------------------------------------------------- */

DumpAtom::DumpAtom(PDPS *ps, int narg, char **arg) : Dump(ps, narg, arg)
{
	if (narg != 5) error->all(FLERR,"Illegal dump command");

	col_labels = NULL;
}

/* ---------------------------------------------------------------------- */

void DumpAtom::init_style()
{
	size_one = 5;

	pack_choice = &DumpAtom::pack_noscale_noimage;
	header_choice = &DumpAtom::header_item;
	write_choice = &DumpAtom::write_noimage;

	char *str;
    str = (char *) "%d %d %g %g %g";
    int n = strlen(str) + 2;
    format = new char[n];
    strcpy(format,str);
    strcat(format,"\n");

	// setup boundary string

	domain->boundary_string(boundstr);

	// setup column's labels string
	col_labels = (char *) "id type x y z";
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_header(bigint ndump)
{
	if (procid == 0) (this->*header_choice)(ndump);
}

/* ----------------------------------------------------------------------
   Choose a pack methdod
------------------------------------------------------------------------- */

void DumpAtom::pack(int *ids)
{
  (this->*pack_choice)(ids);
}

/* ----------------------------------------------------------------------
   Write to dump
------------------------------------------------------------------------- */

void DumpAtom::write_data(int n, double *mybuf)
{
	(this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::header_item(bigint ndump)
{
	fprintf(fp,"ITEM: TIMESTEP\n");
	fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
	fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
	fprintf(fp,"%d \n",ndump);
	fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
	fprintf(fp,"%g %g\n",boxlo[0],boxhi[0]);
	fprintf(fp,"%g %g\n",boxlo[1],boxhi[1]);
	fprintf(fp,"%g %g\n",boxlo[2],boxhi[2]);
	fprintf(fp,"ITEM: ATOMS %s\n",col_labels);
}

/* ---------------------------------------------------------------------- */

void DumpAtom::pack_noscale_noimage(int *ids)
{
	int m,n;

	int *tag = particle->tag;
	int *type = particle->type;
	int *mask = particle->mask;
	double **x = particle->x;
	int nlocal = particle->nlocal;

	m = n = 0;
	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			buf[m++] = tag[i];
			buf[m++] = type[i];
			buf[m++] = x[i][0];
			buf[m++] = x[i][1];
			buf[m++] = x[i][2];
			if (ids) ids[n++] = tag[i];
		}
	}
}

/* ---------------------------------------------------------------------- */

void DumpAtom::write_noimage(int n, double *mybuf)
{
	int m = 0;
	for (int i = 0; i < n; i++) {
		fprintf(fp,format,
            static_cast<int> (mybuf[m]), static_cast<int> (mybuf[m+1]),
            mybuf[m+2],mybuf[m+3],mybuf[m+4]);
		m += size_one;
	}
}
