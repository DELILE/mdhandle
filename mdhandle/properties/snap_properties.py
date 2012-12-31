# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Calculates summary properties for a simulation snap.

"""

import sys

import numpy as np

from mdhandle import units
from mdhandle import vectors
from mdhandle.logger import Logger
import mdhandle.utilities as util

# TODO: Consider refactoring into Snap module.
# TODO: Add optional mass weighting throughout.
# TODO: Add function to add new calculations  
# TODO: register if snap selections change between calcultions.

# -----------------------------------------------------------------------------

logger = Logger()

# -----------------------------------------------------------------------------


class SnapProperties(object):
    """
    Performs summary calculations on snap file
    (e.g. scalar pressure, mass density, etc.).

    Unlike calculations in :mod:`~mdhandle.properties.binned_properties`
    results are for all atoms in the snap (with any selections applied) and are
    not calculated on a grid.

    Parameters
    ----------
    snap : :class:`mdhandle.snap.Snap`
        Snap object used for calculations.

    Attributes
    -----------
    snap : :class:`mdhandle.snap.Snap`
    volume : float
        Volume occupied by snap atoms.
    epair : float
        Total epair in snap.
    ke : float
        Total kinetic energy in snap.
    count : float
        Total number of atoms in snap.
    mass : float
        Total mass in snap.
    force : float
        Total force on snap atoms.
    com : :class:`numpy.ndarray`, (3x1)
        Snap centre of mass position.
    com_velocity : :class:`numpy.ndarray`, (3x1)
        Snap centre of mass velocity.
    temp : float
        Snap temperature.
    stress : :class:`numpy.ndarray`, (6,1)
        Total stress tensor on snap.
    pressure : float
        Scalar pressure acting on snap.
    _KBOLTZ : float
        Boltzman constant in units consistent with :attr:`SnapProperties.snap`.

    """

    def __init__(self, snap):
        self.snap = snap
        self.snap.gather_data()

        self.volume = None
        self._KBOLTZ = units.KBOLTZ[snap.meta['units']]

    def set_volume(self, volume):
        """
        Sets the volume over which density calculations are defined.

        May be defined by user based on another calculation, but cannot be
        easily set by maximum coordinates of atoms within snap because may be
        large amounts of empty space within the system.

        Parameters
        ----------
        volume : float
            Volume of system used for density calculations.
            Units should be consistent with snap units.

        """
        self.volume = float(volume)

    def get_volume_density(self, prop_name):
        """
        Calculates the density of given ``prop_name`` per unit
        volume (:attr:`SnapProperties.volume`).

        Parameters
        -----------
        prop_name : string
            Name of snap property.

        Returns
        -------
        vol_dens : float
            Volume density of total snap property provided.
            Returns ``0.`` if ``prop_name`` doesn't exist in 
            :class:`SnapProperties``.
            
        """
        if self.volume is None or self.volume == 0:
            return 0.

        if hasattr(self, prop_name):
            return getattr(self, prop_name) / self.volume
        else:
            logger.warning('Snap property %s does not exist.' % prop_name)
            return 0.

    def get_mass_norm(self, prop_name):
        """
        Calculates the value of given property per-unit-mass
        in the snap.

        Parameters
        -----------
        prop_name : string
            Name of snap property.

        Returns
        -------
        mass_norm : float
            Per-unit-mass value of property provided.
            Returns ``0.`` if ``prop_name`` doesn't exist in 
            :class:`SnapProperties``.

        """
        if hasattr(self, 'mass') is not True:
            self.mass_total()

        if self.mass != 0:
            if hasattr(self, prop_name):
                return getattr(self, prop_name) / self.mass
            else:
                logger.warning('Snap property %s does not exist.' % prop_name)
                return 0.
        else:
            return 0

    def epair_total(self, epair_name='epair'):
        """
        Total potential energy (epair) in snap.

        Parameters
        -----------
        epair_name : string, [default='epair']
            Name of epair within snap.

        Returns
        --------
        self.epair : float
            Total potential energy.

        """
        self.epair = self.snap.get_scalar(epair_name).sum()
        return self.epair

    def count_total(self):
        """
        Counts number of atoms in snap.

        Returns
        -------
        self.count : int
            Total number of atoms in snap.

        """
        self.count = self.snap.get_scalar('id').size
        return self.count

    def ke_total(self, from_velocity=False, ke_name='ke',
                                            velocity_name='velocity',
                                            mass_name='mass'):
        """
        Total kinetic energy in snap.

        Parameters
        ----------
        from_velocity : boolean [default='False']
            Calculate kinetic energy from atomwise velocity.
        ke_name : string, [default='ke']
            Name of per-atom kinetic energy in snap.
        velocity_name : string [default='velocity']
            Name of atomwise velocity in snap.
        mass_name : string [default='mass']
            Name of atomwise mass in snap.

        Returns
        -------
        self.ke : float
            Total kinetic energy.

        """
        
        vel_comps = self.snap.vectors[velocity_name]['comps']
        
        if from_velocity is True or ke_name not in self.snap.col_mapping:

            if ke_name not in self.snap.col_mapping:
                logger.user_message(
                             'Kinetic energy not in snapshot. Using velocity.')

            u = self.snap.meta['units']
            units_factor = units.MASS[u]*units.VELOCITY[u]**2 / units.ENERGY[u]

            self.ke = units_factor * (0.5*self.snap.get_scalar((mass_name))*
                       vectors.sq_magnitude(self.snap.get_vector(vel_comps))
                                     ).sum(axis=0)
        else:
            self.ke = self.snap.get_scalar(ke_name).sum()

        return self.ke

    def mass_total(self):
        """
        Total mass of snapshot.

        Returns
        -------
        self.mass : float
            Total mass in snap.

        """
        self.mass = self.snap.get_scalar('mass').sum()
        return self.mass

    def force_total(self):
        """
        Net force on snap.

        Returns
        -------
        self.force : :class:`numpy.ndarray`, (3x1)
            Net force on snap.

        """
        f = self.snap.get_vector(('fx', 'fy', 'fz'))
        f = util.clean(f, ~np.isfinite(f), 0.)
        self.force = f.sum(axis=0)
        return self.force

    def other_total(self, prop_name):
        """
        Total of a general scalar property.

        Do not use for vector properties.

        Parameters
        ----------
        prop_name : string
            Name of property to be totaled.

        Returns
        -------
        other_total : float
            Total value in snap.  If ``prop_name`` does not exist
            in :attr:`SnapProperties.snap` returns ``0.``.

        """
        if prop_name in self.snap.col_mapping:
            return self.snap.get_scalar(prop_name).sum()
        else:
            logger.warning('Snap property %s does not exist.' % prop_name)
            return 0.

    def temp_total(self, correct=True):
        """
        Calculates temperature of snap.

        Parameters
        ----------
        correct : boolean, [default=True]
            If ``True``, streaming velocity is removed from the temperature
            calculation

        Returns
        -------
        self.temp : float
            Temperature of snap.

        """
        if not hasattr(self, 'count'):
            self.count_total()
        if not hasattr(self, 'com_velocity'):
            self.com_pos_vel()

        vel = self.snap.get_vector(('vx', 'vy', 'vz'))
        mass = self.snap.get_scalar('mass')

        u = self.snap.meta['units']
        units_factor = units.MASS[u]*units.VELOCITY[u]**2 / units.ENERGY[u]

        if correct is False:
            self.temp = units_factor * 2.0*( 0.5*mass*
                                             (vel*vel).sum(axis=1)).sum() / \
                                             (3*self.count*self._KBOLTZ)

        # Correcting for average velocity of selection
        else:
            vel_corr = vel - self.com_velocity
            self.temp = units_factor * 2.0*(0.5*mass*
                                                (vel_corr*vel_corr).sum(axis=1)
                                            ).sum() / \
                                            (3*self.count*self._KBOLTZ)

        return self.temp

    def stress_total(self):
        """
        Total stress tensor on snap.

        Returns
        -------
        self.stress : :class:`numpy.ndarray`, (6x1)
            Total stress tensor on snap, or ``0.`` if run into problems.

        """
        if self.volume is None or self.volume == 0:
            return 0.

        s = self.snap.get_symm_tensor(
                                 self.snap.symm_tensors['stress']['comps'])
        s = util.clean(s, ~np.isfinite(s), 0.)

        self.stress = s.sum(axis=0) / self.volume
        return self.stress

    def pressure_total(self):
        """
        Scalar pressure in snap based on stress tensor.

        Returns
        -------
        self.pressure : float
            Scalar pressure in snap or 0. if ran into problems.

        """
        if not hasattr(self, 'stress'):
            self.stress_total()

        self.pressure = -(1./3.)*(self.stress[0]+self.stress[3]+self.stress[5])
        return self.pressure

    def com_pos_vel(self, weight='mass'):
        """
        Calculates the position and centre of mass (CoM) of snap.
        
        Atomic position can be weighted by mass or by any other atomwise
        scalar.

        Parameters
        -----------
        weight : string
            Name of scalar used to weigh atom positions in calculation of 
            centre of mass. [default='mass']

        Returns
        -------
        com : :class:`numpy.ndarray`
            Snap centre of mass position.
        com_velocity : :class:`numpy.ndarray`
            Snap centre of mass velocity.

        """
        if not hasattr(self, weight):
            if weight == 'mass':
                weighting_total = self.mass_total()
            else:
                if weight in self.snap.meta['col_mapping']:
                    weighting_total = self.other_total(weight)
                else:
                    logger.error('Weighting scalar %s is not present' % weight)
                    raise Exception

        xyz = self.snap.get_vector(('x', 'y', 'z'))
        weighting_scalar = self.snap.get_scalar(weight)
        atom_vel = self.snap.get_vector(('vx', 'vy', 'vz'))

        self.com = (xyz.T*weighting_scalar.T).T.sum(axis=0) / weighting_total
        self.com_velocity = ((atom_vel.T*weighting_scalar.T).T.sum(axis=0)) \
                                                               / weighting_total
        return self.com, self.com_velocity



# =============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
