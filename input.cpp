/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

// system library
#include "ctype.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"

// pdps library
// top-level class
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
//#include "multiscale.h"
#include "neighbor.h"
#include "output.h"
#include "parallel.h"
#include "particle.h"
#include "pair.h"
#include "read_data.h"
#include "run.h"
#include "update.h"

// accessory-level class
#include "create_box.h"
#include "create_particle.h"
#include "velocity.h"

using namespace PDPS_NS;

#define MAXLINE 2048       // Maximum string length
#define DELTA 4            // For reallocation of memory, if the size of the array need to be expand

/* ---------------------------------------------------------------------- */

Input::Input(PDPS *ps, int argc, char **argv) : Pointers(ps) 
{
	procid = parallel->procid;
	//procid = 0;

	// Initialization
	line = copy = work = NULL;
	line = new char[MAXLINE];
	copy = new char[MAXLINE];
	work = new char[MAXLINE];

	narg = maxarg = 0;
	arg = NULL;
	variable = NULL;
	command = NULL;
}

/* ---------------------------------------------------------------------- */

Input::~Input()
{
	delete line;
	line = NULL;
	delete copy;
	copy = NULL;
	delete work;
	work = NULL;

	memory->sfree(arg);
}

/* ---------------------------------------------------------------------- */

void Input::file()
{
	int m, n;

		// read a line from input script
		// if line ends in continuation char '&', concatenate next line(s)
		// n = length of line including str terminator, 0 if end of file         
		// m = position of last printable char in line or -1 if blank line
	while (1) {
		m = 0;
		if (procid == 0) {
			while (1) {
				if(fgets(&line[m],MAXLINE-m,infile) == NULL) n = 0;     // If the input script is read completely
				else {		
					n = strlen(line) + 1;                              // The end of string 
					if (screen) {
						fprintf(screen, "%s\n", line);
						//fflush(screen);
					}
					if (logfile) {
						fprintf(logfile, "%s\n", line);
					}
				}
				if(n == 0) break;
				m = n - 2;
				while (m >= 0 && isspace(line[m])) m--;                // Scan from the end of the string
				if(m < 0 || line[m] != '&') break;
			} // while (1) for one line
		} // if (procid == 0) 
		MPI_Bcast(&n,1,MPI_INT,0,mworld);    // Broad cast n 
		if (n == 0) break;
		MPI_Bcast(line,n,MPI_CHAR,0,mworld);
		parse();                                     // examine each key word in one line
		if (command == NULL) continue;
		if(execute_command()) {
			char str[MAXLINE];
			sprintf(str,"Unknown command: %s",line);
			error->all(FLERR,str);
		}
	} // while (1) for the whole script
}

/* ----------------------------------------------------------------------
   parse copy of command line
   strip comment = all chars from # on
   command = first word
   narg = # of args
   arg[] = individual args
   treat text between single/double quotes as one arg
------------------------------------------------------------------------- */

void Input::parse()
{
	//arg = NULL;
	strcpy(copy,line);

	// strip any # comment by replacing it with 0
	// do not strip # inside single/double quotes

	char *ptr = copy;
	while (*ptr) {
		if (*ptr == '#') {
			*ptr = '\0';
			break;
		}
		ptr++;
	}

	char *next;
	command = nextword(copy,&next);
	if (command == NULL) return;

	// point arg[] at each subsequent arg
	// nextword() inserts zeroes in copy to delimit args
	narg = 0;
	ptr = next;
	//free(arg);
	while (ptr) {
		if (narg == maxarg) {
			maxarg += DELTA;
			arg = (char **) memory->srealloc(arg,maxarg*sizeof(char *),"Input: arg");
		}
		arg[narg] = nextword(ptr,&next);
		if (!arg[narg]) break;
		narg++;
		ptr = next;
	}
}

/* ----------------------------------------------------------------------
   find next word in str
   insert 0 at end of word
   ignore leading whitespace
   return ptr to start of word
   return next = ptr after word or NULL if word ended with 0
   return NULL if no word in string
------------------------------------------------------------------------- */

char *Input::nextword(char *str, char **next)
{
  char *start,*stop;

  start = &str[strspn(str," \t\n\v\f\r")];
  if (*start == '\0') return NULL;

  stop = &start[strcspn(start," \t\n\v\f\r")];
  
  //n = stop - start;
  
  //strncpy(start,start,n);             // Extract the first word from the string
  if (*stop == '\0') *next = NULL;
  else *next = stop+1;
  *stop = '\0';
  return start;
}

/* ---------------------------------------------------------------------- */
int Input::execute_command()
{
    // flag = 1:  unkown command
	// flag = 0:  valid command or comment line

	int flag = 1;
	
	if(!strcmp(command, "analyze")) analyze();
	else if (!strcmp(command, "boundary")) boundary();
	else if (!strcmp(command, "compute")) compute();
	else if (!strcmp(command, "create_box")) create_box();
	else if (!strcmp(command, "create_particle")) create_particle();
	else if (!strcmp(command, "density")) density();
	else if (!strcmp(command, "dimension")) dimension();
	else if (!strcmp(command, "dump")) dump();
	else if (!strcmp(command, "dump_modify")) dump_modify();
	else if (!strcmp(command, "energy")) energy();
	else if (!strcmp(command, "fix")) fix();
	else if (!strcmp(command, "group")) group_command();
	else if (!strcmp(command, "integrate")) integrate();
	else if (!strcmp(command, "lattice")) lattice();
	else if (!strcmp(command, "mass")) mass();
	else if (!strcmp(command, "multiscale")) multi_command();
	else if (!strcmp(command, "multi_run")) multi_run();
	else if (!strcmp(command, "neigh_modify")) neigh_modify();
	else if (!strcmp(command, "neighbor")) neighbor_command();
	else if (!strcmp(command, "neighbor_slave")) neighbor_setslave();
	else if (!strcmp(command, "pair_coeff")) pair_coeff();
	else if (!strcmp(command, "pair_style")) pair_style();
	else if (!strcmp(command, "particle_style")) particle_style();
	else if (!strcmp(command, "processors")) processors();
	else if (!strcmp(command, "radius")) radius();
	else if (!strcmp(command, "rho")) rho();
	else if (!strcmp(command, "read_data")) read_data_command();
	else if (!strcmp(command, "region")) region();
	else if (!strcmp(command, "reset_timestep")) reset_timestep();
	else if (!strcmp(command, "run")) run();
	else if (!strcmp(command, "save")) save();
	else if (!strcmp(command, "timestep")) timestep();
	else if (!strcmp(command, "thermo")) thermo();
	else if (!strcmp(command, "thermo_style")) thermo_style();
	else if (!strcmp(command, "unanalyze")) unanalyze();
	else if (!strcmp(command, "undump")) undump();
	else if (!strcmp(command, "unfix")) unfix();
	else if (!strcmp(command, "units")) units();
	else if (!strcmp(command, "velocity")) velocity();

	else flag = 0;        // return 1, if there is no matched command

	// return if command was listed above
	if (flag) return 0; 

	// unrecognized command
	return -1;

}

/* ---------------------------------------------------------------------- */

void Input::analyze()
{
	modify->add_analyze(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::boundary()
{
	if (domain->box_exist)
		error->all(FLERR,"Boundary command after simulation box is defined");
	domain->set_boundary(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::compute()
{
	modify->add_compute(narg, arg);
}

void Input::create_box()
{
	CreateBox createbox(ps);
	createbox.command(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::create_particle()
{
	if (narg < 3) {
		error->all(FLERR,"Illegal create_particle command");
	}
  	CreateParticle createparticle(ps);              // initialize create_particle
	createparticle.command(narg,arg);         // 
}

/* ---------------------------------------------------------------------- */

void Input::density()
{
	particle->set_density(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::dimension()
{
	if (narg > 1) error->all(FLERR,"Illegal dimension command");
	domain->dim = atoi(arg[0]);
	if (domain->dim < 0 || domain->dim > 3) error->all(FLERR,"Illegal dimension number");
}

/* ---------------------------------------------------------------------- */

void Input::dump()
{
	output->add_dump(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::dump_modify()
{
	output->modify_dump(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::energy()
{
	particle->set_energy(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::fix()
{
	modify->add_fix(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::group_command()
{
	group->assign(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::integrate()
{
	update->create_integrate(narg, arg);
}


/* ---------------------------------------------------------------------- */

void Input::lattice()
{
	domain->set_lattice(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::mass()
{
  if (narg != 2) error->all(FLERR,"Illegal mass command");
  
    particle->set_mass(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::multi_command()
{
	//multiscale->set_style(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::multi_run()
{
	//multiscale->run(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::neigh_modify()
{
	if ((narg % 2) != 0) 
		error->all(FLERR,"Illegal neighbor_modify command");
	neighbor->modify_params(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::neighbor_command()
{
	neighbor->settings(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::pair_coeff()
{
	force->pair[force->npairs-1]->set_coeff(narg, arg);	
}

/* ---------------------------------------------------------------------- */

void Input::pair_style()
{
	force->create_pair(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::particle_style()
{
	if (narg < 1) error->all(FLERR, "Illegal particle_style command");
	if (domain->box_exist) {
		error->all(FLERR, "particle_style command should set before simulation box is defined");
	}
	particle->create_particle_type(arg[0], narg-1, &arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::processors()
{
	if (domain->box_exist) {
		error->all(FLERR, "Processors command after simulation box is defined");
	}
	parallel->set_processors(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::radius()
{
	particle->set_radius(narg, arg);
}

/* ---------------------------------------------------------------------- */
void Input::rho()
{
	particle->set_rho(narg, arg);
}
/* ---------------------------------------------------------------------- */

void Input::read_data_command()
{
	ReadData readdata(ps);
	readdata.command(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::region()
{
  //if (domain->box_exist)
    //error->all(FLERR,FLERR,"Boundary command after simulation box is defined");
	domain->add_region(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::reset_timestep()
{
	update->reset_timestep(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::run()
{
	Run run(ps);
	run.command(narg, arg);
}
/* ---------------------------------------------------------------------- */

void Input::save()
{
	particle->save_particle(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::timestep()
{
	if (narg != 1) 
		error->all(FLERR,"Illegal timestep command");
	update->dt = atof(arg[0]);
	if (update->dt < 0) 
		error->all(FLERR,"Time step cannot be negative");
}

/* ---------------------------------------------------------------------- */

void Input::thermo()
{
	//if (narg != 1) error->universe_one(screen, "Illegal timestep command");
	output->set_thermo(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::thermo_style()
{
	output->create_thermo(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::unanalyze()
{
	if (narg != 1) error->all(FLERR, "Illegal unanalyze command");
	modify->delete_analyze(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::undump()
{
	if (narg != 1) error->all(FLERR,"Illegal undump command");
	output->delete_dump(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::unfix()
{
	if (narg != 1) error->all(FLERR, "Illegal undump command");
	modify->delete_fix(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::units()
{
	if (narg > 1) {
		error->all(FLERR,"Too many arguments for units command");
	}
	update->set_units(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::velocity()
{
	// Error message
	if (narg < 5) error->all(FLERR, "Illegal velocity command");
	Velocity velocity(ps);
	velocity.command(narg, arg);
}

/* ---------------------------------------------------------------------- */

void Input::neighbor_setslave()
{
	neighbor->setslave(narg, arg);
}