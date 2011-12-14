# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Plotting velocity distribution of snapshot agains the Maxwell-Boltzmann
distribution at that temperature

**Credits**:

* http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.maxwell.html#scipy.stats.maxwell
* http://mathworld.wolfram.com/MaxwellDistribution.html

"""

import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell

from mdhandle import units
from mdhandle import vectors
from mdhandle.logger import Logger

# ----------------------------------------------------------------------------

logger = Logger()

# ----------------------------------------------------------------------------


def plot_velocity_distribution(snap, temp):
    """
    Plots histogram of atomic speeds against analytical Maxwell-Boltzman
    distribution for input temperature, ``temp``.

    Parameters
    ----------
    snap : :class:`mdhandle.snap.Snap`
        :class:`~mdhandle.snap.Snap` object.
    temp : float
        Temperature for analytical Maxwell distribution

    Returns
    --------
    0
        Successfully completed velocity distribution plot.
    -1
        An error has occured and velocity distribution plot was not generated.

    """
    kb = units.KBOLTZ[snap.meta['units']]
    speed = vectors.magnitude(snap.get_vector(('vx', 'vy', 'vz')))*\
            units.VELOCITY[snap.meta['units']]**2

    if np.unique(snap.get_scalar('mass')).size > 1:
        logger.warning('Cannot handle polyatomic system')
        return -1

    # Taking mass as constant based on value at index 0
    mass = snap.get_scalar('mass')[0]*units.MASS[snap.meta['units']]

    # See: <http://mathworld.wolfram.com/MaxwellDistribution.html>
    # Energy factor is applied to bring to SI units for calculation.
    # if metal or other unit system.
    a = np.sqrt(kb*temp/mass * units.ENERGY[snap.meta['units']])
    scale = 1 / np.sqrt(a)

    # Maxwell-Boltzmann RV based on temperature
    maxwell_rv = maxwell(loc=0, scale=scale)

    x = np.linspace(0, np.minimum( maxwell_rv.dist.b, speed.max() ))

    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, maxwell_rv.pdf(x), linewidth=2)

    n, bins, patches = ax.hist(speed, 100, normed=True,
                                           facecolor='red', alpha=0.75)
    return 0

#=============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
