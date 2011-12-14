/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator 

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov 

   See the README file in the top-level LAMMPS directory. 

   ----------------------------------------------------------------------- 

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/ 

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany 

   See the README file in the USER-CUDA directory. 

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_cg_cmm_coul_debye_cuda.h"
#include "pair_cg_cmm_coul_debye_cuda_cu.h"
#include "cuda_data.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "cuda_neigh_list.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "cuda.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCGCMMCoulDebyeCuda::PairCGCMMCoulDebyeCuda(LAMMPS *lmp) : PairCGCMMCoulCut(lmp)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

	allocated2 = false;
	cg_type_double = NULL;
	cuda->shared_data.pair.cudable_force = 1;
	cuda->setSystemParams();
}

/* ----------------------------------------------------------------------
   remember pointer to arrays in cuda shared data
------------------------------------------------------------------------- */

void PairCGCMMCoulDebyeCuda::allocate()
{
	if(! allocated) PairCGCMMCoulCut::allocate();
	int n = atom->ntypes;
	if(! allocated2)
	{
		allocated2 = true;
		
  
  		memory->create(cg_type_double,n+1,n+1,"paircg:cgtypedouble");
  		
		cuda->shared_data.pair.cut     = cut_lj;
		cuda->shared_data.pair.cut_coul= cut_coul;
		cuda->shared_data.pair.coeff1  = lj1;
		cuda->shared_data.pair.coeff2  = lj2;
		cuda->shared_data.pair.coeff3  = lj3;
		cuda->shared_data.pair.coeff4  = lj4;
		cuda->shared_data.pair.coeff5  = cg_type_double;
		cuda->shared_data.pair.offset  = offset;
		cuda->shared_data.pair.special_lj  = force->special_lj;
		cuda->shared_data.pair.special_coul  = force->special_coul;
	}
  	for (int i = 1; i <= n; i++) {
      for (int j = i; j <= n; j++) {
        cg_type_double[i][j] = cg_type[i][j];
        cg_type_double[j][i] = cg_type[i][j];
      }
    }
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulDebyeCuda::compute(int eflag, int vflag)
{
	if (eflag || vflag) ev_setup(eflag,vflag);
	if(eflag) cuda->cu_eng_vdwl->upload();
	if(eflag) cuda->cu_eng_coul->upload();
	if(vflag) cuda->cu_virial->upload();

	Cuda_PairCGCMMCoulDebyeCuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom);

    if(not cuda->shared_data.pair.collect_forces_later)
    {
	  if(eflag) cuda->cu_eng_vdwl->download();
	  if(eflag) cuda->cu_eng_coul->download();
	  if(vflag) cuda->cu_virial->download();
    }
	
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulDebyeCuda::settings(int narg, char **arg)
{
	PairCGCMMCoulCut::settings(narg, arg);
	cuda->shared_data.pair.cut_global = (F_FLOAT) cut_lj_global;
	cuda->shared_data.pair.cut_coul_global = (F_FLOAT) cut_coul_global;
	cuda->shared_data.pair.kappa = (F_FLOAT) kappa;
}

/* ---------------------------------------------------------------------- */

void PairCGCMMCoulDebyeCuda::coeff(int narg, char **arg)
{
	PairCGCMMCoulCut::coeff(narg, arg);
	allocate();
}

void PairCGCMMCoulDebyeCuda::init_style()
{
	MYDBG(printf("# CUDA PairCGCMMCoulDebyeCuda::init_style start\n"); )
  // request regular or rRESPA neighbor lists

  int irequest;
 
  if (update->whichflag == 0 && strstr(update->integrate_style,"respa")) {

  } 
  else 
  {
  	irequest = neighbor->request(this);
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->cudable = 1;
    //neighbor->style=0; //0=NSQ neighboring
  }

  cuda->shared_data.pppm.qqrd2e=force->qqrd2e;
  cut_respa=NULL;
  if (force->newton) error->warning(FLERR,"Pair style uses does not use \"newton\" setting. You might test if \"newton off\" makes the simulation run faster.");

  MYDBG(printf("# CUDA PairCGCMMCoulDebyeCuda::init_style end\n"); )
}

void PairCGCMMCoulDebyeCuda::init_list(int id, NeighList *ptr)
{
	MYDBG(printf("# CUDA PairCGCMMCoulDebyeCuda::init_list\n");)
	PairCGCMMCoulCut::init_list(id, ptr);
	#ifndef CUDA_USE_BINNING
	// right now we can only handle verlet (id 0), not respa
	if(id == 0) cuda_neigh_list = cuda->registerNeighborList(ptr);
	// see Neighbor::init() for details on lammps lists' logic
	#endif
	MYDBG(printf("# CUDA PairCGCMMCoulDebyeCuda::init_list end\n");)
}

void PairCGCMMCoulDebyeCuda::ev_setup(int eflag, int vflag)
{
	int maxeatomold=maxeatom;
	PairCGCMMCoulCut::ev_setup(eflag,vflag);

  if (eflag_atom && atom->nmax > maxeatomold) 
	{delete cuda->cu_eatom; cuda->cu_eatom = new cCudaData<double, ENERGY_FLOAT, x > ((double*)eatom, & cuda->shared_data.atom.eatom , atom->nmax  );}

  if (vflag_atom && atom->nmax > maxeatomold) 
	{delete cuda->cu_vatom; cuda->cu_vatom = new cCudaData<double, ENERGY_FLOAT, yx > ((double*)vatom, & cuda->shared_data.atom.vatom , atom->nmax, 6  );}
	
}


