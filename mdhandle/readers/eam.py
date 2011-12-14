# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Parses and plots EAM potential files (i.e. ``Pt_u3.eam`` distributed 
with LAMMPS).

Not suitable to alloy potential files.

**Credits** :

* MATLAB script by Arun K Subramaniyan
* Pt EAM functions (universal 3): SM Foiles et al, Phys Rev B, 33, 7983 (1986)

.. note::
   EAM files in LAMMPS are calibrated to function in units metal
   so distance is in angstroms, and energy is in eV.

"""

import sys

import numpy as np
import pylab as plt

# -----------------------------------------------------------------------------


class EAM(object):
    """
    Plots components of EAM potential.

    Designed for format of ``Pt_u3.eam`` distributed with LAMMPS.
    Not suitable for alloy EAM potential files.

    Parameters
    -----------
    eam_fn : string
        Name of EAM potential file.
    verbose : boolean
        If ``True``, verbose user I/O.

    Attributes
    -----------
    eam_fn : string
        EAM file name.
    verbose : boolean
        If ``True``, verbose user I/O.
    credit : string
        Credit for EAM potential data from EAM file.
    atomic_num : int
        Atomic number of material from EAM file.
    lattice_constant : float
        Lattice constant of material from EAM file.
    crystal_structure : string
        Crystal structure type from EAM file.
    num_embedding : int
        Number of entries for embedding function.
    num_elec_dens : int
        Number of entries for electron density function.
    dr_elec_dens : float
        Delta in radial distance between entries for electron density function.
    dr_pair : float
        Delta in radial distance between entries for pair energy function.
    embedding : :class:`numpy.ndarray`
        Embedding function values as a function of radial distance.
    pair : :class:`numpy.ndarray`
        Pair energy function as a function of radial distance.
    density : :class:`numpy.ndarray`
        Electron density as a function of radial distance.
    _header_size : int
        Number of lines in EAM file header.
    _number_cols : int
        Number of columns in EAM file.
    _f : file
        File handle for EAM file.

    """

    def __init__(self, eam_fn, verbose=False):
        # File info
        self._header_size = 3
        self._number_cols = 5

        self.eam_fn = eam_fn
        self.verbose = verbose

        # Initializing attributes.
        self._f = None                         
        self.cutoff = 0.
        self.atomic_num = 0.
        self.atomic_mass = 0.
        self.crystal_structure = ''
        self.lattice_constant = 0.
        self.credit = ''

        self.num_pair_points = 0.
        self.num_elec_dens = 0.    

        self.density = 0.
        self.embedding  = 0.
        self.pair = 0.

        self.dr_pair = 0.
        self.dr_elec_dens = 0.

    def read(self):
        """
        Read and store data from EAM file.

        """
        # Read header
        self._f = open(self.eam_fn, 'r')

        self.credit = self._f.readline()

        line2 = self._f.readline().split()
        self.atomic_num = line2[0]
        self.atomic_mass = line2[1]
        self.lattice_constant = line2[2]
        self.crystal_structure = line2[3]

        line3 = self._f.readline().split()

        # Number of embedding is equal to num_elec_dens
        self.num_embedding = line3[0]

        self.num_elec_dens = line3[0]
        self.dr_elec_dens = line3[1]

        self.num_pair_points = line3[2]
        self.dr_pair = line3[3]

        self.cutoff = line3[4]

        self._f.close()

        # Big block of data from file - broken up later
        data = np.loadtxt(self.eam_fn, dtype=np.float,
                                       skiprows=self._header_size)

        # Embedding energy is first block (200 rows, 5 cols read left to right)
        self.embedding = data[:self.num_embedding/self._number_cols, :].ravel()

        # Pair energy is second block (200 rows, 5 cols read left to right)
        # TODO: simplify into separate statements.
        self.pair = data[self.num_embedding/self._number_cols:
                            self.num_embedding/self._number_cols
                                + self.num_pair_points/self._number_cols,
                                                                     :].ravel()

        # Density is third block (200 rows, 5 cols read left to right)
        self.density = data[self.num_embedding/self._number_cols+
                            self.num_pair_points/self._number_cols:, :].ravel()

    def plot(self):
        """
        Plot EAM components ('pair', 'density', 'embedding')

        """
        # Number of samples along radius
        radius_scaffold = np.arange(0, self.num_embedding)

        plt.figure()

        plt.title('pair')
        plt.xlabel('$r$')
        plt.ylabel('$\sqrt{r \cdot \phi}$')
        plt.plot(radius_scaffold*self.dr_pair, self.pair)

        plt.figure()
        plt.title('density')
        plt.xlabel('$r$')
        plt.ylabel('$\rho$')
        plt.plot(radius_scaffold*self.dr_elec_dens, self.density)

        plt.figure()
        plt.title('embedding')
        plt.xlabel('$r$')
        plt.ylabel('$F$')
        plt.plot(radius_scaffold*self.dr_pair, self.embedding)


# =============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
