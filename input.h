/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator
   
   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_INPUT_H
#define PS_INPUT_H

#include "stdio.h"
#include "pointers.h"

namespace PDPS_NS {

class Input : protected Pointers {
public:
	int narg;                    // # of command args
	char **arg;                  // parsed args for command
	class Variable *variable;    // defined variables

	Input(class PDPS *, int, char **);
	~Input();
	void file();                   // process input file

private:
	int procid;                     // proc ID
	char *command;               // ptr to current command
	int maxarg;                  // max # of args in arg
	char *line,*copy,*work;      // input line & copy of it


	void parse();                      // parse an input text line
	char *nextword(char *, char **);   // find next word in string, with quotes

	int execute_command();             // execute a single command
	 
	void analyze();
	void atom_style();               // PDPS commands
	void boundary();
	void compute();
	void create_box(); 
	void create_particle();
	void density();
	void dimension();
	void dump();
	void dump_modify();
	void energy();
	void fix();
	void group_command();
	void integrate();
	void lattice();
	void mass();
	void multi_command();
	void multi_run();
	void neigh_modify();
	void neighbor_command();
	void pair_coeff();
	void pair_style();
	void particle_style();
	void processors();
	void radius();
	void rho();
	void read_data_command();
	void region();
	void reset_timestep();
	void run();
	void save();
	void timestep();
	void thermo();
	void thermo_style();
	void units();
	void velocity();

	void communicate();

	void pair_modify();

	void pair_write();
	

	void restart();
	void run_style();

	void unanalyze();
	void uncompute();
	void undump();
	void unfix();
	void neighbor_setslave();
  
};

}

#endif

/* ERROR/WARNING messages:

E: Label wasn't found in input script

Self-explanatory.

E: Input line too long: %s

This is a hard (very large) limit defined in the input.cpp file.

E: Unknown command: %s

The command is not known to LAMMPS.  Check the input script.

E: Another input script is already being processed

Cannot attempt to open a 2nd input script, when the original file is
still being processed.

E: Cannot open input script %s

Self-explanatory.

E: Unbalanced quotes in input line

No matching end double quote was found following a leading double
quote.

E: Input line quote not followed by whitespace

An end quote must be followed by whitespace.

E: Invalid variable name

Variable name used in an input script line is invalid.

E: Substitution for illegal variable

Input script line contained a variable that could not be substituted
for.

E: Input line too long after variable substitution

This is a hard (very large) limit defined in the input.cpp file.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open logfile %s

The LAMMPS log file specified in the input script cannot be opened.
Check that the path and name are correct.

E: Angle_coeff command before simulation box is defined

The angle_coeff command cannot be used before a read_data,
read_restart, or create_box command.

E: Angle_coeff command before angle_style is defined

Coefficients cannot be set in the data file or via the angle_coeff
command until an angle_style has been assigned.

E: Angle_coeff command when no angles allowed

The chosen atom style does not allow for angles to be defined.

E: Angle_style command when no angles allowed

The chosen atom style does not allow for angles to be defined.

E: Atom_style command after simulation box is defined

The atom_style command cannot be used after a read_data,
read_restart, or create_box command.

E: Bond_coeff command before simulation box is defined

The bond_coeff command cannot be used before a read_data,
read_restart, or create_box command.

E: Bond_coeff command before bond_style is defined

Coefficients cannot be set in the data file or via the bond_coeff
command until an bond_style has been assigned.

E: Bond_coeff command when no bonds allowed

The chosen atom style does not allow for bonds to be defined.

E: Bond_style command when no bonds allowed

The chosen atom style does not allow for bonds to be defined.

E: Boundary command after simulation box is defined

The boundary command cannot be used after a read_data, read_restart,
or create_box command.

E: Dihedral_coeff command before simulation box is defined

The dihedral_coeff command cannot be used before a read_data,
read_restart, or create_box command.

E: Dihedral_coeff command before dihedral_style is defined

Coefficients cannot be set in the data file or via the dihedral_coeff
command until an dihedral_style has been assigned.

E: Dihedral_coeff command when no dihedrals allowed

The chosen atom style does not allow for dihedrals to be defined.

E: Dihedral_style command when no dihedrals allowed

The chosen atom style does not allow for dihedrals to be defined.

E: Dimension command after simulation box is defined

The dimension command cannot be used after a read_data,
read_restart, or create_box command.

E: Improper_coeff command before simulation box is defined

The improper_coeff command cannot be used before a read_data,
read_restart, or create_box command.

E: Improper_coeff command before improper_style is defined

Coefficients cannot be set in the data file or via the improper_coeff
command until an improper_style has been assigned.

E: Improper_coeff command when no impropers allowed

The chosen atom style does not allow for impropers to be defined.

E: Improper_style command when no impropers allowed

The chosen atom style does not allow for impropers to be defined.

E: KSpace style has not yet been set

Cannot use kspace_modify command until a kspace style is set.

E: Mass command before simulation box is defined

The mass command cannot be used before a read_data, read_restart, or
create_box command.

E: Min_style command before simulation box is defined

The min_style command cannot be used before a read_data, read_restart,
or create_box command.

E: Newton bond change after simulation box is defined

The newton command cannot be used to change the newton bond value
after a read_data, read_restart, or create_box command.

E: Package command after simulation box is defined

The package command cannot be used afer a read_data, read_restart, or
create_box command.

E: Package cuda command without USER-CUDA installed

The USER-CUDA package must be installed via "make yes-user-cuda"
before LAMMPS is built.

E: Pair_coeff command before simulation box is defined

The pair_coeff command cannot be used before a read_data,
read_restart, or create_box command.

E: Pair_coeff command before pair_style is defined

Self-explanatory.

E: Pair_modify command before pair_style is defined

Self-explanatory.

E: Pair_write command before pair_style is defined

Self-explanatory.

E: Processors command after simulation box is defined

The processors command cannot be used after a read_data, read_restart,
or create_box command.

E: Run_style command before simulation box is defined

The run_style command cannot be used before a read_data,
read_restart, or create_box command.

E: Units command after simulation box is defined

The units command cannot be used after a read_data, read_restart, or
create_box command.

*/
