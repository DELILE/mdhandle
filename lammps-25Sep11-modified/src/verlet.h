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

#ifdef INTEGRATE_CLASS

IntegrateStyle(verlet,Verlet)

#else

#ifndef LMP_VERLET_H
#define LMP_VERLET_H

#include "integrate.h"

namespace LAMMPS_NS {

class Verlet : public Integrate {
 public:
  Verlet(class LAMMPS *, int, char **);
  virtual ~Verlet() {}
  void init();
  void setup();
  void setup_minimal(int);
  void run(int);
  void cleanup();

 protected:
  int triclinic;                    // 0 if domain is orthog, 1 if triclinic
  int torqueflag,erforceflag;
  int e_flag,rho_flag;

  void force_clear();
};

}

#endif
#endif
