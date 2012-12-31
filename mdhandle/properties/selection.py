# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Enable the ability to easily create selections of a sub-set of atoms within
LAMMPS output.

"""

import sys

import numpy as np

import mdhandle.vectors as vectors
from mdhandle.logger import Logger

# TODO: Create abstract base class that doesn't do anything aside from store
#       mask, and temporal_lock.
# TODO: Boolean operations on selections.
#       - Need to accommodate selections with different temporal_lock

# -----------------------------------------------------------------------------

logger = Logger()

# -----------------------------------------------------------------------------


class Selection(object):
    """
    Class used for selecting a subset of atoms from a snap.

    Parameters
    -----------
    masking_function : callable
        Callable taking :class:`mdhandle.snap.Snap` as input.
        Must return a boolean :class:`numpy.ndarray` (Nx1).

        Should take a single argument - the `Snap` object on which the
        selection is applied.
    temporal_lock : boolean, [default=True]
        If ``True``, selection remains fixed once it is defined.

    Attributes
    -----------
    temporal_lock : boolean
        If ``True``, :attr:`Selection.masked` only defined once by first call.
    masking_function : callable
        Function used to define the selection.
    masked : :class:`numpy.ndarray`
        Result from mask operation.  Atoms with ``True`` entries are included.
    first_time : boolean
        If ``True``, :meth:`Selection.calculate_mask` has not been called.

    """

    def __init__(self, masking_function, temporal_lock=True):
        self.temporal_lock = temporal_lock

        self.masking_function = masking_function
        self.masked = None
        self.first_time = True

    def calculate_mask(self, snap):
        """
        Creates mask for selecting atoms.

        Wrapper for masking_function with logic for controlling changes
        via :attr:`Selection.temporal_lock`.

        Parameters
        ----------
        snap : :class:`mdhandle.snap.Snap`
            Snap object for mask.

        Returns
        -------
        self.masked : Boolean :class:`numpy.ndarray` (Nx1)
            Atoms with true entries are included in the selection.

        """
        if self.first_time is True:
            self.masked = self.masking_function(snap)
            self.first_time = False
        elif self.temporal_lock is False:
            self.masked = self.masking_function(snap)

        return self.masked

# ------------------- Selection Functions -------------------------------------

def mask_exists(mask_array):
    """
    Selection function for the case where the mask array already exists
    from a previous calculation.
    
    Parameters
    ----------
    mask_array : `numpy.ndarray`
        Masking array in the correct format for selection.
    
    Returns
    -------
    masking_function : callable
        Function for selecting liquid atoms.
    
    """
    def masking_function(s):
        return mask_array
    return masking_function

def liquid_only(threshold, fluid_atom_type, coord_name='coord'):
    """
    Selects liquid atoms using coordination number.
    Liquid atoms have coordination numbers larger than the threshold.

    Parameters
    ----------
    threshold : int
        Lower threshold on coordination number for liquid.
        Typical value: 40 for simulation with cutoff of :math:`4.0 \sigma`.
    fluid_atom_type : iterable
        Types of atoms which are fluids in the simulation.
        Data type of entry should be `'int'`
        If only a single atom type is desired: ``fluid_atom_type=[1,]``
    coord_name : string, [default='coord']
        Name of coordination number in snap.

    Returns
    -------
    masking_function : callable
        Function for selecting liquid atoms.
    """
    fluid = fluid_only(fluid_atom_type)

    def masking_function(s):
        coord = s.get_scalar(coord_name, raw_data=True)
        return fluid(s)*(coord > threshold)

    return masking_function


def fluid_only(fluid_atom_types):
    """
    Function factory for selecting only fluid atoms.

    Parameters
    -----------
    fluid_atom_type : iterable
        Types of atoms which are fluids in the simulation.
        Data type of entry should be ``'int'``
        If only a single atom type is desired: ``fluid_atom_type=[1,]``

    Returns
    -------
    masking_function : callable
        Function for selecting fluid atoms.
    """
    def masking_function(s):
        mask = np.ones(s.meta['num_atoms'], dtype=bool)
        snap_type = s.get_scalar('type', raw_data=True)
        for i in fluid_atom_types:
            mask *= (snap_type == int(i))

        return mask

    return masking_function


def coordinate_less_than(coord_name, threshold):
    """
    Function factory for selectino of atoms with coordinate
    less than a given value.

    Parameters
    ----------
    coord_name : string, { 'x' | 'y' | 'z'}
        Name of coordinate to be tested.
    threshold : float
        Maximum value of coordinate

    Returns
    -------
    masking_function : callable
        Function for selecting atoms with coordinates less than ``threshold``.
    """
    def masking_function(s):
        coord_test = s.get_scalar(coord_name, raw_data=True)
        return coord_test < threshold

    return masking_function


def coordinate_greater_than(coord_name, threshold):
    """
    Function factory for selection of atoms with coordinates
    greater than a given value

    Parameters
    -----------
    coord_name : string, { 'x' | 'y' | 'z'}
        Name of coordinate to be tested.
    threshold : float
        Minimum value of coordinate

    Returns
    -------
    masking_function : callable
        Function for selecting atoms with coordinates greater than
        ``threshold``.
    """
    def masking_function(s):
        coord_test = s.get_scalar(coord_name, raw_data=True)
        return coord_test > threshold

    return masking_function


def box(corner_low, dimensions, inside=True):
    """
    Function factory for selecting atoms in a box defined by low
    corner anddimensions.

    Box is aligned to coordinate directions.

    Parameters
    -----------
    corner_low : :class:`numpy.ndarray`, (3x1)
        Coordinates of the low corner of the box.
    dimensions : :class:`numpy.ndarray`, (3x1)
        Length of box in three coordinate directions.
    inside : boolean, [default=True]
        If ``True``, selects atoms inside the box.

    Returns
    -------
    masking_function : callable
        Function for selecting atoms with coordinates inside 3D box.
    """
    # TODO: Generalize to box that is not aligned to coordinate directions.

    # Converting inputs to Numpy array if not already
    corner_low = np.array(corner_low)
    dimensions = np.array(dimensions)

    corner_high = corner_low + dimensions

    def masking_function(s):
        xyz = s.get_vector(('x', 'y', 'z'), raw_data=True)
        mask = np.prod(xyz < corner_high, axis=1)*\
               np.prod(xyz > corner_low,  axis=1)

        if inside is True:
            return mask
        else:
            return np.logical_not(mask)

    return masking_function


def sphere(centre, radius, inside=True):
    """
    Factory function for selection of atoms inside a
    sphere defined by centre and radius.

    Parameters
    -----------
    centre : :class:`numpy.ndarray`, (3x1)
        Coordinates of the centre of the sphere.
    radius : float
        Radius of the sphere.
    inside : boolean, [default=True]
        If ``True``, selects atoms inside the sphere.

    Returns
    --------
    masking_function : callable
        Function for selecting atoms with coordinates inside 3D sphere.
    """
    # Converting input vectors to Numpy array if not already
    centre = np.array(centre)

    
    def masking_function(s):
        xyz = s.get_vector(('x', 'y', 'z'), raw_data=True)
        mask = vectors.sq_magnitude(xyz-centre) < radius**2
        if inside is True:
            return mask
        else:
            return np.logical_not(mask)

    return masking_function


def opposite_selection(filter_function_input):
    """
    Factory function inverting existing selection mask.

    Parameters
    ----------
    filter_function_input : callable
        Existing selection function.

    Returns
    -------
    reversed_filter_function : callable
        Inverted ``filter_function_input``.
    """
    def reversed_filter_function(s):
        return np.logical_not(filter_function_input(s))

    return reversed_filter_function


# =============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
