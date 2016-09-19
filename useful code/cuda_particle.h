/* ----------------------------------------------------------------------
   PDPS - Particle Dynamics Parallel Simulator

   Copyright (2012) reserved by Lingqi Yang. 
   Email: ly2282@columbia.edu

   See the README file in the top-level PDPS directory.
------------------------------------------------------------------------- */

#ifndef PS_CUDA_PARTICLE_H
#define PS_CUDA_PARTICLE_H

#include "pointers.h"
#include "particle.h"

namespace PDPS_NS {

class CUDAParticle : protected Pointers {
public:
	

	CUDAParticle(class PDPS *);
	~CUDAParticle();

	static double LoadFactor;
	int   HashTableSize;
	int  *devMapArray;
	unsigned int *devHashKey;
	unsigned int *devHashVal;
	unsigned int *devTry;

	virtual void SetMap();
	void data_atoms(int n, char *buf);


	// permutation table for atom reordering
	unsigned long long *PermuKey;
	int *PermuTableA; // A: per-Atom
	int *PermuTableT; // T: Topology


	void Grow(int, int);
	void data_particles(int n, char *buf);
	virtual void PinHostArray();
	virtual void UnpinHostArray();
	// per-atom properties
	int *devTag;
	int *devType;
	int *devMask;
	int *devMolecule;
	double *devMass;
	double *devCoordX;
	double	*devCoordY;
	double	*devCoordZ;
	double	*devForceX;
	double	*devForceY;
	double	*devForceZ;
	double	*devVeloX;
	double	*devVeloY;
	double	*devVeloZ;
	double	*devVirialXX;
	double	*devVirialYY;
	double	*devVirialZZ;
	double	*devVirialXY;
	double	*devVirialYZ;
	double	*devVirialXZ;
	float	*devImageX;
	float	*devImageY;
	float	*devImageZ;
	double	*devEBond;
	double	*devEAngle;
	double	*devEDihed;
	double	*devEImprop;
	double	*devEPair;
	float *devRCoordX;
	float *devRCoordY;
	float *devRCoordZ;
	float *devRVeloX;
	float *devRVeloY;
	float *devRVeloZ;
	float4 *devMCoord;
	float4 *devMVelo;

	// pointer to pre-allocated device buffer
	double	*devHostCoord;
	double	*devHostVelo;
	double	*devHostForce;
	int	*devHostImage;
	int *devHostTag;
	int *devHostType;
	int *devHostMask;
	double	*devHostMassType;
	int	*devHostSpecial;
	int *devHostMolecule;

	int *devHostNBond, *devHostBondAtom, *devHostBondType;
	int *ptrHostNBond, *ptrHostBondAtom, *ptrHostBondType;

	// pointer to host memory, e.g. atom->x
	double	*ptrHostCoord;
	double	*ptrHostVelo;
	double	*ptrHostForce;
	int	*ptrHostImage;
	int *ptrHostTag;
	int *ptrHostType;
	int *ptrHostMask;
	double	*ptrHostMassType;
	int	*ptrHostSpecial;
	int *ptrHostMolecule;

	// topology
	static int	 MaxBondNum;
	size_t		 BondTablePitch;
	int			 BondTablePadding;
	int			*devNBond;
	int2		*devBonds;
	int2		*devBondsMapped;
	// *****those below are dummy entry***
	static int	 MaxAngleNum;
	size_t		 AngleTablePitch;
	int			 AngleTablePadding;
	int         *devNAngle;
	int			*devAngles,
				*devAngleType;
	static int	 MaxDihedNum;
	size_t		 DihedTablePitch;
	int			 DihedTablePadding;
	int         *devNDihed;
	int			*devDiheds,
				*devDihedType;
	static int	 MaxImpropNum;
	size_t		 ImpropTablePitch;
	int			 ImpropTablePadding;
	int         *devNImprop;
	int			*devImprops,
				*devImpropType;
	static int	 MaxExclNum;
	int			 ExclTablePadding;
	int4        *devNExcl;
	int         *devHostNExcl;
	int			*devExclTable; // table in column fasion

	protected:
		const static int    GrowthFactorI;
		const static double GrowthFactorM;
		int Allocated;
	//	CUDASortPlan<16, u64, int> Sorter;
};

}

#endif

