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

#ifdef FIX_CLASS

FixStyle(GPU,FixGPU)

#else

#ifndef LMP_FIX_GPU_H
#define LMP_FIX_GPU_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGPU : public Fix {
 public:
  FixGPU(class LAMMPS *, int, char **);
  ~FixGPU();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);
  double memory_usage();

 private:
  int _gpu_mode;
  double _particle_split;
};

}

#endif
#endif
