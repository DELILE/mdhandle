# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Wrapper to compiled serial LAMMPS library.  

Changes made to the LAMMPS SVN v25Sep11 to allow Python access via the LAMMPS library build to enable access to other types of atom-wise output in ordered way as is originally possible with coordinates only.

Methods other than *new* :meth:`~mdhandle.lammps.lammps.put_coords` and
:meth:`~mdhandle.lammps.lammps.lazy_extract_compute` remain unchanged.

See LAMMPS distribution for documentation for further information.

* http://lammps.sandia.gov/doc/Section_python.html


**Credits**:

* ``lammps.py`` file from LAMMPS v25Sep11 is largely unchanged. Modified to allow for ordered output of atomwise quantities.
* See changes and patch file within LAMMPS source for modifications.

"""

# -----------------------------------------------------------------------------

# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   http://lammps.sandia.gov, Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under 
#   the GNU General Public License.
#
#   See the README file in the top-level LAMMPS directory.
# -------------------------------------------------------------------------

import os
from ctypes import *

import mdhandle

LMPINT = 0
LMPDOUBLE = 1
LMPIPTR = 2
LMPDPTR = 3
LMPDPTRPTR = 4

class lammps:
  def __init__(self, args=None):

    # Need to account for location of compiled LAMMPS library
    # at the root of the mdhandle egg when installed.
    # _lammps_serial.so is located in site-packages folder due to HACK 
    # build below.
    location = os.path.dirname(
                    os.path.dirname(
                            os.path.dirname(os.path.abspath(__file__))))

    # attempt to load parallel library first, serial library next
    # could provide caller a flag to choose which library to load
    #
    # See: http://docs.python.org/library/ctypes.html#loading-shared-libraries
    #
    try:
      self.lib = CDLL(os.path.join(location, "_lammps.so"))
    except:
      try:
        self.lib = CDLL( os.path.join(location, "_lammps_serial.so"))
      except:
        raise StandardError,"Could not load LAMMPS dynamic library"

    # create an instance of LAMMPS
    # don't know how to pass an MPI communicator from PyPar
    # no_mpi call lets LAMMPS use MPI_COMM_WORLD
    # cargs = array of C strings from args
    
    if args:
      args.insert(0,"lammps.py")
      narg = len(args)
      cargs = (c_char_p*narg)(*args)
      self.lmp = c_void_p()
      self.lib.lammps_open_no_mpi(narg,cargs,byref(self.lmp))
    else:
      self.lmp = c_void_p()
      self.lib.lammps_open_no_mpi(0,None,byref(self.lmp))
      # could use just this if LAMMPS lib interface supported it
      # self.lmp = self.lib.lammps_open_no_mpi(0,None)

  def __del__(self):
    if self.lmp: self.lib.lammps_close(self.lmp)

  def close(self):
    self.lib.lammps_close(self.lmp)
    self.lmp = None

  def file(self,file):
    self.lib.lammps_file(self.lmp,file)

  def command(self,cmd):
    self.lib.lammps_command(self.lmp,cmd)

  def extract_global(self,name,type):
    if type == LMPDOUBLE:
      self.lib.lammps_extract_global.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_global(self.lmp,name)
      return ptr[0]
    if type == LMPINT:
      self.lib.lammps_extract_global.restype = POINTER(c_int)
      ptr = self.lib.lammps_extract_global(self.lmp,name)
      return ptr[0]
    return None

  def extract_atom(self,name,type):
    if type == LMPDPTRPTR:
      self.lib.lammps_extract_atom.restype = POINTER(POINTER(c_double))
      ptr = self.lib.lammps_extract_atom(self.lmp,name)
      return ptr
    if type == LMPDPTR:
      self.lib.lammps_extract_atom.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_atom(self.lmp,name)
      return ptr
    if type == LMPIPTR:
      self.lib.lammps_extract_atom.restype = POINTER(c_int)
      ptr = self.lib.lammps_extract_atom(self.lmp,name)
      return ptr
    return None

  def extract_compute(self,id,style,type):
    if type == 0:
      if style > 0: return None
      self.lib.lammps_extract_compute.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr[0]
    elif type == 1:
      self.lib.lammps_extract_compute.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr
    elif type == 2:
      self.lib.lammps_extract_compute.restype = POINTER(POINTER(c_double))
      ptr = self.lib.lammps_extract_compute(self.lmp,id,style,type)
      return ptr
    return None

  # in case of global datum, free memory for 1 double via lammps_free()
  # double was allocated by library interface function
  
  def extract_fix(self,id,style,type,i=0,j=0):
    if type == 0:
      if style > 0: return None
      self.lib.lammps_extract_fix.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_bix(self.lmp,id,style,type,i,j)
      result = ptr[0]
      self.lib.lammps_free(ptr)
      return result
    elif type == 1:
      self.lib.lammps_extract_fix.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,i,j)
      return ptr
    elif type == 2:
      self.lib.lammps_extract_fix.restype = POINTER(POINTER(c_double))
      ptr = self.lib.lammps_extract_fix(self.lmp,id,style,type,i,j)
      return ptr
    return None

  # free memory for 1 double or 1 vector of doubles via lammps_free()
  # for vector, must copy nlocal returned values to local c_double vector
  # memory was allocated by library interface function
  
  def extract_variable(self,name,group,type):
    if type == 0:
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      result = ptr[0]
      self.lib.lammps_free(ptr)
      return result
    if type == 1:
      self.lib.lammps_extract_global.restype = POINTER(c_int)
      nlocalptr = self.lib.lammps_extract_global(self.lmp,"nlocal")
      nlocal = nlocalptr[0]
      result = (c_double*nlocal)()
      self.lib.lammps_extract_variable.restype = POINTER(c_double)
      ptr = self.lib.lammps_extract_variable(self.lmp,name,group)
      for i in xrange(nlocal): result[i] = ptr[i]
      self.lib.lammps_free(ptr)
      return result
    return None

  def get_natoms(self):
    return self.lib.lammps_get_natoms(self.lmp)

  def get_coords(self):
    nlen = 3 * self.lib.lammps_get_natoms(self.lmp)
    coords = (c_double*nlen)()
    self.lib.lammps_get_coords(self.lmp,coords)
    return coords

  
  def put_coords(self, coords):
    """
    Modifies atom coordinates in memory of running LAMMPS instance.

    Assume coords is an array of c_double, as created by get_coords().

    Parameters
    ----------
    coords :  :class:`ctypes.c_couble`
        :class:`ctypes` array of atom coordinates.
        Array is length ``nlen``, the number of atoms in the LAMMPS simulation.

    """
    self.lib.lammps_put_coords(self.lmp, coords)

  def lazy_extract_compute(self, id, type):
    """
    Modification for mdhandle.
    
    Alternative to extract_compute(...) method which returns
    a copy of array rather than pointer.
    
    This was previously only possible for coordinates via get_coords(...) 
    method.

    Parameters
    ----------
    id : string
        LAMMPS name for variable to extract.
    type : int
        Defines the shape of the variable.
        { 0 -> 'scalar' |  1 -> 'vector' | 2 -> 'tensor' |  4 -> 'symm tensor'}

    """

    # See library.cpp for info.  Other styles are not implemented as they
    # were in extract_compute(...)
    style = 1
    
    # Scalars: IN LAMMPS, per-atom vector is equivalent to scalar in mdhandle.
    if type == 0:
      nlen = self.lib.lammps_get_natoms(self.lmp)
      vec = (c_double*nlen)()
      self.lib.lammps_lazy_extract_compute(self.lmp, vec, id,style,type)
      return vec

    # For non-scalars, LAMMPS uses `array` rather than `vec`

    # vectors
    elif type == 1:
      nlen = 9 * self.lib.lammps_get_natoms(self.lmp)
      array = (c_double*nlen)()
      self.lib.lammps_lazy_extract_compute(self.lmp, array, id,style,type)
      return array

    # tensors
    elif type == 2:
      nlen = 9 * self.lib.lammps_get_natoms(self.lmp)
      array = (c_double*nlen)()
      self.lib.lammps_lazy_extract_compute(self.lmp, array, id,style,type)
      return array

    # symm tensors
    elif type == 3:
      nlen = 6 * self.lib.lammps_get_natoms(self.lmp)
      array = (c_double*nlen)()
      self.lib.lammps_lazy_extract_compute(self.lmp, array, id,style,type)
      return array

    # None of conditions are satisfied - returning None
    else:
      return None

