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

#ifdef COMPUTE_CLASS

ComputeStyle(displace/atom,ComputeDisplaceAtom)

#else

#ifndef LMP_COMPUTE_DISPLACE_ATOM_H
#define LMP_COMPUTE_DISPLACE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDisplaceAtom : public Compute {
 public:
  ComputeDisplaceAtom(class LAMMPS *, int, char **);
  ~ComputeDisplaceAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double **displace;
  char *id_fix;
  class Fix *fix;
};

}

#endif
#endif