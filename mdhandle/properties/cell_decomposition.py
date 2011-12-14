# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Binning is a usefull technique for efficiently calculating all
distances between a set of coordinates, when you are only interested in
the distances below a given cutoff. The algorithm consists of two major
steps:

1. Divide the given set of coordinates into bins on a regular grid in space.
2. Calculate the distances (or other usefull things) between coordinates
   in neighbouring bins.

**Credits**: Inspiration for binning algorithm from binning module in ``mol_mod`` module in ``MD-Tracks`` package

"""
#TODO: Handle case where atoms are outside the primary simulation cell.
#           - max and min bins are reliant on sim_cell lengths not
#             maximum spread of atoms.
#           - Include optional wrap flag to control application of PBC.

import sys

import numpy as np

from mdhandle.logger import Logger

# -----------------------------------------------------------------------------

logger = Logger()

# -----------------------------------------------------------------------------


class SparseBinning(object):
    """
    3D space in the simulation domain is divided into a sparse grid with
    atoms assigned into the appropriate bin using on a nearest-grid-point (NGP)
    algorithm.

    Each cell in the grid is called a bin.
    
    Each bin contains a set of atoms ids, or ``None`` if no atoms are present.
    
    All bins are uniquely defined by their indices ``(i,j,k)``.  

    This implementation works with sparse bins as a bin is only created
    as required.

    Parameters
    ----------
    snap : :class:`mdhandle.snap.Snap`
        Snap object to be sorted into cells.
    grid_spacing : :class:`numpy.ndhandle`
        Grid spacing along each coordinate direction.
        Grid spacing must be an even divisor of the primary simulation
        cell dimensions.
        ``[default=np.ones(3,dtype=np.float)]``
    use_selection : boolean    [default=False]
        If ``True``, selections are applied to snap object when getting atom
        coordinates.

    Attributes
    -----------
    snap : :class:`mdhandle.snap.Snap`
         Snap object.  Active dataset must be of type ``'atoms'``.
    grid_spacing : :class:`numpy.ndarray` (3x1)
        Grid spacing along coordinate directions:  ``[dx, dy, dz]``.
        Must be even divisor of the simulation domain in each direction.
    xyz : :class:`numpy.ndarray` (Nx3)
        Atomic positions.
    bins : dict
        Organized by tuples corresponding to cell location
        (i.e. ``(1, 1, 1 )``).  Each dictionary entries contains atom ids
        for atoms located in that grid cell.
        This is the sparse element of implementation as dictionary entries
        only created as needed.
    max_bins_idx : :class:`numpy.ndarray` (3x1)
        Largest grid cell index along each coordinate direction.
    min_bins_idx : :class:`numpy.ndarray` (3x1)
        Smallest grid cell index along each coordinate direction.
    nbins : :class:`numpy.ndarray` (3x1)
        Number of bins along each coordinate direction.

    Notes
    ------

    **Class Attributes**:

    * ``full_compare_indices`` : :class:`numpy.ndarray`
        Deltas required to move through all neighbouring cells.
    * ``N3_compare_indices`` : :class:`numpy.ndarray`
        Deltas required to move through neighbouring cells with the help
        of Newton's Third law such that only half interactions are required
        at each grid cell.

    """

    # -------------- Class variables  -----------------------------------------

    # Displacements to nearest neighbouring cells
    full_compare_indices = np.array([
            (-1, -1, -1), (-1, -1,  0), (-1, -1,  1),
            (-1,  0, -1), (-1,  0,  0), (-1,  0,  1),
            (-1,  1, -1), (-1,  1,  0), (-1,  1,  1),

            ( 0, -1, -1), ( 0, -1,  0), ( 0, -1,  1),
            ( 0,  0, -1), ( 0,  0,  0), ( 0,  0,  1),
            ( 0,  1, -1), ( 0,  1,  0), ( 0,  1,  1),

            ( 1, -1, -1), ( 1, -1,  0), ( 1, -1,  1),
            ( 1,  0, -1), ( 1,  0,  0), ( 1,  0,  1),
            ( 1,  1, -1), ( 1,  1,  0), ( 1,  1,  1),
            ], dtype=np.int)

    N3_compare_indices = np.array([
            (0, 0, 0), (1, 1, 1),
            (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (0, 1, 1), (1, 0, 1), (1, 1, 0),
            (0, 1, -1), (-1, 0, 1), (1, -1, 0),
            (1, 1, -1), (1, -1, -1), (1, -1, 1)
            ], dtype=np.int)

    def __init__(self, snap, grid_spacing=np.ones(3, np.float),
                             use_selection=True):
        self.snap = snap
        self.grid_spacing = grid_spacing

        if snap.grid_type != 'atoms':
            logger.error("Snap must be grid type 'atoms'.")

        # Testing that grid divides evenly into provide sim cell.
        test_divisors = np.mod(self.snap.sim_cell.lengths, grid_spacing)
        if np.all(test_divisors == 0) is False:
            logger.error('Grid spacing must be even divisor of sim cell')

        # Getting atom coordinates - apply selection (if use_selection == True)
        self.xyz = self.snap.get_vector( ('x', 'y', 'z'),
                                                raw_data=(not use_selection)  )
        self.xyz = self.snap.sim_cell.wrap_to_central(self.xyz)

        self.bins = {}
        # -1 from initial calculation as bins are zero-indexed
        self.max_bins_idx = np.floor(
                    (self.snap.sim_cell.origin + self.snap.sim_cell.lengths) /
                    self.grid_spacing).astype(np.int) - 1
        self.min_bins_idx = np.floor(self.snap.sim_cell.origin /
                                            self.grid_spacing).astype(np.int)
        self.nbins = np.floor(self.snap.sim_cell.lengths /
                                            self.grid_spacing).astype(np.int)

    def _ngp(self, r):
        """
        Calculate bin index for position ``r`` using  nearest-grid-point (NGP)
        algorithm.  The grid point is located at the centre of the bin.

        Atoms located exactly at the maximum edge of the simulation
        domain are moved back into the largest bin to avoid off-by-one error.

        Parameters
        ----------
        r : :class:`numpy.ndarray` (Nx3 or 1x3)
            Atom position.

        Returns
        --------
        indices : dict
            Return :class:`numpy.ndarray` with same shape as input position
            containing the bin index for each atom.

        References
        -----------
        * Hockney & Eastwood, Computer Simulation Using Particles, CRC Press,
          1988.

        """
        indices = np.floor(r / self.grid_spacing).astype(np.int)

        # For edge case where particle is EXACTLY at very MAXIMUM must
        # move into maximum .  Otherwise will create a bin that is out of range
        indices[ indices == self.max_bins_idx + 1] = self.max_bins_idx
        return indices

    def get_surrounding(self, r, 
                              deltas=full_compare_indices):
        """
        Iterate over all atoms in the surrounding bins as defined by
        deltas parameter.

        Parameters
        ----------
        r : :class:`numpy.ndarray` (3x1)
            Atom position in ``(x,y,z)`` format.
        deltas : :class:`numpy.ndarray`  [default=full_compare_indices]
            Array of tuples with format ``(dx,dy,dz)`` giving integer deltas
            to neighbouring cells. For example, ``(-1, -1, -1)`` looks in the 
            bin ``dx = -1``, ``dy = -1``, ``dz = -1`` bin from the bin
            containing the position defined by ``r``.

        Returns
        --------
        bin_total : list
            List of ``atom_id`` for atoms in bins surrouding position ``r``.

        """
        # Wrapping to cental cell as per PBC
        r = self.snap.sim_cell.wrap_to_central(r)

        center = self._ngp(r)
        bin_total = []
        for delta in deltas:
            question = center + delta

            # Apply PBC (controlled by self.snap.sim_cell.pbc)
            question = question + \
                                    (question > self.max_bins_idx)*\
                                    -1*self.nbins*self.snap.sim_cell.pbc +\
               (question < self.min_bins_idx)*self.nbins*self.snap.sim_cell.pbc

            # Skipping any bins that are outside min or max bins
            if (np.any(question > self.max_bins_idx) or
                                         np.any(question < self.min_bins_idx)):
                continue

            question = tuple(question)

            bin0 = self.bins.get(question)
            if bin0 is not None:
                bin_total.extend(bin0)

        return bin_total

    def run(self):
        """
        Run cell decomposition for :attr:`SparseBinning.snap` object. 
        
        Cell decomposition is stored in :attr:`SparseBinning.bins`.

        """

        # Loop over atoms and assgn to self.bins
        atom_id = (i for i in xrange(0, self.xyz.shape[0]))

        for idx in self._ngp(self.xyz):
            idx = tuple(idx)
            bin = self.bins.get(idx)
            if bin is None:
                bin = set()
                bin.add(atom_id.next())
                self.bins[idx] = bin
            else:
                bin.add(atom_id.next())

# =============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
