# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2

"""
Handling information about unit cell and periodic boundary conditions.

**Credits**: 

* `MD-Tracks` simulation cell implementation.

"""

# TODO: Incorporate other simulation cell metadata currently stored w/i Snap
#       such as xlo, ylo, etc.
# TODO: Extend beyond orthogonal simulation cell aligned to global XYZ coord.

# -----------------------------------------------------------------------------

import sys

import numpy as np

# -----------------------------------------------------------------------------


class SimCell(object):
    """
    Stores information regarding the simulation cell.

    **Assumptions**:

    * Simulation cell is orthogonal.
    * Simulation cell is 3D.
    * Simulation cell attributes are implictly read-only once created by
      ``__init__()``.

    Parameters
    ----------
    axes : iterable of float
        Elements of iterable define the size of the simulation cell
        in each coordinate direction.
        The simulation cell is assumed to be aligned to the fixed global XYZ
        coordinate directions.
    origin : :class:`numpy.ndarray`, (3x1)
        Origin of simulation domain relative to fixed global XYZ coordinates.
    pbc : boolean :class:`numpy.ndarray` (3x1)
        Defines the periodic boundary conditions along each 
        coordinate direction. If ``True``, PBC invoked at BOTH the
        positive and negative faces of simulation cell.
        [default=np.array[[True, True, True]]

    Attributes
    ----------
    pbc : :class:`numpy.ndarray`, (3x1)
        Defines periodic boundary conditions along each coordinate direction.
    origin : :class:`numpy.ndarray`, (3x1)
        Low corner of simulation cell relative to fixed XYZ frame.
    volume : float
        Volume of simulation cell.
    angles : :class:`numpy.ndarray`, (3x1)
        Corner angles of simulation cell.
    _matrix : :class:`numpy.ndarray`, (3x3)
        Diagonals are lenghts of the simulation cell.
    _boxL : :class:`numpy.ndarray`, (3x1)
        Length of full grid in each coordinate direction.
    _recip : :class:`numpy.ndarray`, (3x3)
        Inverse of :attr:`_matrix`.

    """

    def __init__(self, axes,
                       origin=np.zeros(3),
                       pbc=np.ones(3, dtype=np.bool)):

        # ASSERT: pbc is a vector with length 3.
        assert len(pbc) == 3

        if isinstance(pbc, list):
            pbc = np.array(pbc)

        assert isinstance(pbc, np.ndarray)
        assert pbc.dtype == np.dtype('bool')
        self.pbc = pbc

        # ASSERT: axes is a vector with length 3.
        assert len(axes) == 3
        self._matrix = np.diag(axes).astype(np.float)

        # ASSERT: origin is a vector with length 3.
        assert len(origin) == 3
        assert isinstance(origin, np.ndarray)
        self.origin = origin

        self._boxL = np.diag(self._matrix)

        # ASSERT: All lengths must be postive
        for i in self._boxL:
            assert i > 0

        self.volume = abs(np.linalg.det(self._matrix))
        # ASSERT: Sim cell must have positive volume
        assert self.volume > 0

        # Reciprocal of simulation cell values.
        self._recip = np.linalg.inv(self._matrix)

        # Angles defining simulation cell
        alpha = np.arccos(np.dot(self._matrix[:,1], self._matrix[:,2]) /\
                                            (self._boxL[1] * self._boxL[2]))
        beta = np.arccos(np.dot(self._matrix[:,2], self._matrix[:,0]) /\
                                            (self._boxL[2] * self._boxL[0]))
        gamma = np.arccos(np.dot(self._matrix[:,0], self._matrix[:,1]) /\
                                            (self._boxL[0] * self._boxL[1]))
        self.angles = np.array([alpha, beta, gamma], dtype=np.float)

    def to_fractional(self, full_xyz):
        """
        Returns vector collection of positions coordinates converted
        to fractional values of basic simulation cell dimensions.

        Parameters
        ----------
        full_xyz : :class:`numpy.ndarray` (Nx3)
            Vector collection of raw coordinates.

        Returns
        --------
        frac : :class:`numpy.ndarray` (Nx3)
            Returns vector stack of same size as ``full_xyz``
            normalized by the size of the simulation cell

        """
        frac = np.dot(full_xyz, self._recip)
        return frac

    def from_fractional(self, fractional):
        """
        Returns a vector collection with raw position coordinates
        from fractional positions normalized by simulation cell size.

        Parameters
        -----------
        fractional : :class:`numpy.ndarray` (shape=Nx3)
            Vector collection of fractional coordinates.

        Returns
        --------
        pos : :class:`numpy.ndarray`, (Nx3)
            Returns vector stack of same size as fractional
            with positions back into unmodified absolute coordinates within
            the simulation cell

        See Also
        ---------
        :meth:`to_fractional`

        """
        pos = np.dot(fractional, self._matrix.transpose())
        return pos

    def mic_distance(self, raw_distance):
        """
        Enforce periodic boundary conditions and the minimum
        image convention on interatomic distance (``raw_distance``).

        Parameters
        -----------
        raw_distance : :class:`numpy.ndarray`, (Nx3)
            Vector collection of raw interatomic distance measures.

        Returns
        --------
        min_img : :class:`numpy.ndarray`, (Nx3)
            Returns vector collection of the same shape as ``raw_distance`` 
            with fractional coordinates ranging between [-0.5,0.5] to satisfy
            the minimum image convention.

         """
        frac = self.to_fractional(raw_distance)
        frac -= self.pbc*frac.round()
        min_img = self.from_fractional(frac)
        return min_img

    def wrap_to_central(self, pos):
        """
        Wraps a vector collection of atom positions into the central periodic
        image or primary simulation cell.

        Parameters
        ----------
        pos : :class:`numpy.ndarray`, (Nx3)
            Vector collection of atom positions.

        Returns
        -------
        wrap : :class:`numpy.ndarray`, (Nx3)
            Returns atomic positions wrapped into the primary simulation
            cell, or periodic image.

        """
        wrap = pos + (pos < self.origin)*self._boxL*self.pbc \
                      - (pos > self._boxL)*self._boxL*self.pbc
        return wrap


#=============================================================================

def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
