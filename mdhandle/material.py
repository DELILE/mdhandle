# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Hold information about materials in MD simulation such as constants used
to define interactomic potential functions.

"""
# TODO: Add EAM file handling for EAMMaterial

import sys

# -----------------------------------------------------------------------------

class LJMaterial(object):
    """
    Class for Lennard-Jones 12-6 material.

    All values should be given in non-dimensional LJ units.
    (See : :mod:`mdhandle.units`)

    Parameters
    -----------
    name : string
        Name of material (e.g. argon)
    sigma : float
        Molecular diameter.
        Corresponds to zero of LJ 12-6 function.
    epsilon : float
        Strength of LJ 12-6 interaction.
        Minimum potential energy of LJ 12-6 potential.
    mass : float
        Atomic mass.

    Examples
    ----------
    >>> argon = ljMaterial('Ar', 3.405, 119.8, 39.95)

    """
    def __init__(self, name, sigma, epsilon, mass):
        self.sigma = sigma
        self.epsilon = epsilon
        self.name = name
        self.mass = mass


class EAMMaterial(object):
    """
    Class for embedded-atom-model (EAM) material

    All values should be in either LAMMPS metal units, or SI.

    Parameters
    -----------
    name : string
        Name of material (e.g. 'Platinum')
    mass : float
        Atomic mass.
    lat_str : string, [default=None]
        Name of crystal structure.  If not crystalline, set to ``None``.
    lat_cons : float, [default=0.0]
        Lattice constant.
        EAM materials usually crystalline, so lattice constant applicable.
        If non-crystalline material, set to ``0.0``

    Examples
    ---------
    >>> platinum = eamMaterial('Pt', 3.9200, 'fcc', 195.078,)

    """
    def __init__(self, name, mass, lat_str=None, lat_cons=0.):
        self.name = name
        self.lat_cons = lat_cons
        self.lat_str = lat_str
        self.mass = mass

# =============================================================================


def main(argv=None):
    pass

if __name__ == "__main__":
    sys.exit(main())
