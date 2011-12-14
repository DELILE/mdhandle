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

FixStyle(rigid,FixRigid)

#else

#ifndef LMP_FIX_RIGID_H
#define LMP_FIX_RIGID_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRigid : public Fix {
 public:
  FixRigid(class LAMMPS *, int, char **);
  virtual ~FixRigid();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void initial_integrate(int);
  void post_force(int);
  virtual void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int, int);
  virtual double compute_scalar();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  void pre_neighbor();
  int dof(int);
  void deform(int);
  void reset_dt();
  double compute_array(int, int);

 protected:
  int me,nprocs;
  double dtv,dtf,dtq;
  double *step_respa;
  int triclinic;

  int nbody;                // # of rigid bodies
  int *nrigid;              // # of atoms in each rigid body
  double *masstotal;        // total mass of each rigid body
  double **xcm;             // coords of center-of-mass of each rigid body
  double **vcm;             // velocity of center-of-mass of each
  double **fcm;             // force on center-of-mass of each
  double **inertia;         // 3 principal components of inertia of each
  double **ex_space,**ey_space,**ez_space;
                            // principal axes of each in space coords
  double **angmom;          // angular momentum of each in space coords
  double **omega;           // angular velocity of each in space coords
  double **torque;          // torque on each rigid body in space coords
  double **quat;            // quaternion of each rigid body
  int *imagebody;           // image flags of xcm of each rigid body
  double **fflag;           // flag for on/off of center-of-mass force
  double **tflag;           // flag for on/off of center-of-mass torque
  double **langextra;       // Langevin thermostat forces and torques

  int *body;                // which body each atom is part of (-1 if none)
  double **displace;        // displacement of each atom in body coords

  double **sum,**all;       // work vectors for each rigid body
  int **remapflag;          // PBC remap flags for each rigid body

  int extended;             // 1 if any particles have extended attributes
  int dorientflag;          // 1 if particles store dipole orientation
  int qorientflag;          // 1 if particles store quat orientation

  int *eflags;              // flags for extended particles
  double **qorient;         // rotation state of ext particle wrt rigid body
  double **dorient;         // orientation of dipole mu wrt rigid body

  double tfactor;           // scale factor on temperature of rigid bodies
  int langflag;             // 0/1 = no/yes Langevin thermostat

  int tempflag;             // NVT settings
  double t_start,t_stop;
  double t_period,t_freq;
  int t_chain,t_iter,t_order;

  int pressflag;            // NPT settings
  double p_start,p_stop;
  double p_period,p_freq;
  int p_chain;

  class RanMars *random;
  class AtomVecEllipsoid *avec_ellipsoid;

                            // bitmasks for eflags
  int INERTIA_POINT,INERTIA_SPHERE,INERTIA_ELLIPSOID;
  int ORIENT_DIPOLE,ORIENT_QUAT;
  int OMEGA,ANGMOM,TORQUE;

  void no_squish_rotate(int, double *, double *, double *, double);
  void set_xv();
  void set_v();
};

}

#endif
#endif
