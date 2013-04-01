# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2

"""
Handle units and unit conversions for output from LAMMPS in lj, metal
and si units.

For description of unit systems in LAMMPS see
http://lammps.sandia.gov/doc/units.html


"""

# TODO: Integrate conversion cosntants and helper classes
# TODO: Add other LAMMPS unit systems (e.g. 'real', 'cms', etc.)
# TODO: Source of constants --> linkup with values used in LAMMPS.

import sys

import numpy as np

# -----------------------------------------------------------------------------

UNITS_ALLOWED = ['lj', 'metal', 'si']

# -----------------------------------------------------------------------------
#    Constants
# -----------------------------------------------------------------------------

AMU = 1.66053886e-27            # [kg/amu]
eV = 1.60217646e-19             # [J/eV]

# Boltzmann's constant for each unit type - can be used as is within 
# calculations in each unit type when other unit correction factors are used.
KBOLTZ = {}
KBOLTZ['si'] = 1.3806504e-23             # [J/K]
KBOLTZ['lj'] = 1.0                       # non-dimensional LAMMPS LJ units
KBOLTZ['metal'] = KBOLTZ['si'] / eV      # [eV / K]


AVAGADRO = 6.0221415e23                  # non-dimensional
LIGHT_SPEED = {}
LIGHT_SPEED['si'] = 2.99792458e8         # [m/s]

# -----------------------------------------------------------------------------
#    Unit Conversion Factors
# -----------------------------------------------------------------------------

# Used to produce multiplicative factors to produce the right units.
# The purpose of each factor is to bring the units into form where 
# calculations can be done.  For 'LJ' and 'si', this isn't necessary, so 1.0 is # used throughout.  However, for 'metal' care must be used when converting to
# SI.

# See mdhandle.properties.snap_properties for examples of use.

MASS = {}
MASS['lj'] = 1.0        # non-dimensional LAMMPS LJ units
MASS['si'] = 1.0        # [kg]
MASS['metal'] = AMU     # [kg/AMU] - LAMMPS metal units

LENGTH = {}
LENGTH['lj'] = 1.0                # non-dimensional LAMMPS LJ units
LENGTH['si'] = 1.0                # [m]
LENGTH['metal'] = 1e-10           # [m / Ang]

TIME = {}
TIME['lj'] = 1.0            # non-dimensional LAMMPS LJ units
TIME['si'] = 1.0            # [s]
TIME['metal'] = 1e-12       # [s / ps] - LAMMPS metal units

FORCE = {}
FORCE['lj'] = 1.0                                # non-dimensional LJ units
FORCE['si'] = 1.0                                # [N]
FORCE['metal'] = eV / LENGTH['metal']           # [J/eV * Ang / m]

ENERGY = {}
ENERGY['lj'] = 1.0                    # non-dimensional LAMMPS LJ units
ENERGY['si'] = 1.0                    # [J]
ENERGY['metal'] = eV                  # [J / eV] - LAMMPS metal units

VELOCITY = {}
VELOCITY['lj'] = 1.0                                 # non-dimensional LJ units
VELOCITY['si'] = 1.0                                 # [m/s]
VELOCITY['metal'] = LENGTH['metal'] / TIME['metal']  # [m / ang * ps/s]

PRESSURE = {}
PRESSURE['lj'] = 1.0                      # non-dimensional LAMMPS LJ units
PRESSURE['si'] = 1.0                      # [Pa]
PRESSURE['metal'] = 1e5                   # [Pa / bar] (convert to si)

# -----------------------------------------------------------------------------
#    Helper Unit Conversion Classes - Useful from the cmd line
# -----------------------------------------------------------------------------

class LJUnits(object):
    """
    Helper object used to convert LJ units to SI units

    Each method returns a parameter converted to SI units from
    non-dimensional LJ units as defined by the sigma, epsilon and mass.

    Parameters
    ----------
    sigma : float, [default=3.405]
        :math:`\sigma` parameter for 12-6 LJ potential in angstroms.
    epsilon_kb : float, [default=119.8]
        :math:`\epsilon / k_b`, where :math:`\epsilon` is energy parameter
        for 12-6 LJ potential in K.
    m : float, [default=39.95]
        Atomic mass in AMU.

    Attributes
    -----------
    sigma : float
        :math:`\sigma` parameter for 12-6 LJ potential in m.
    epsilon_kb : float
        :math:`\epsilon / k_b`, where :math:`\epsilon` is energy parameter
        for 12-6 LJ potential in K.
    m : float
        Atomic mass in Kg.
    epsilon : float
        :math:`\epsilon` energy parameter for 12-6 LJ potential in J.

    Methods
    --------
    real_length(len_star)
        Convert non-dim LJ length to m.
    real_temp(temp_star)
        Convert non-dim LJ temperature to Kelvins.
    real_energy(nrg_star)
        Convert non-dim LJ energy to Joules.
    real_press(press_star)
        Convert non-dim LJ pressure to Pa.
    real_num_dens(nd_star)
        Convert non-dim LJ number density to [atoms/m^3]
    real_time(t_star)
        Convert non-dim LJ time to sec.
    real_velocity(v_star)
        Convert non-dim LJ velocity to [m/s].
    real_force(f_star)
        Convert non-dim LJ force to N.
    real_surface_tension(gamma_star)
        Convert non-dim LJ surface tension to [N/m].
    real_viscosity(mu_star)
        Convert non-dim LJ viscosity to [Pa.s].

    References
    ------------
    http://lammps.sandia.gov/doc/units.html

    """
    def __init__(self, sigma=3.405, epsilon_kb=119.8, m=39.95):
        self.sigma = sigma*1.e-10       # [m]
        self.epsilon_kb = epsilon_kb    # [K]
        self.m = m*AMU                  # [kg]

        self.epsilon = self.epsilon_kb * KBOLTZ['si']    # [J]

    def real_length(self, len_star):
        """
        Convert non-dim LJ length to m.

        Parameters
        -----------
        len_star : float
            Non-dimensional length.

        Returns
        --------
        len_si : float
            Length in m.

        """
        len_si = len_star*self.sigma
        return len_si

    def real_temp(self, temp_star):
        """
        Convert non-dim LJ temperature to Kelvins.

        Parameters
        -----------
        temp_star : float
            Non-dimensional temperature.

        Returns
        --------
        temp_si : float
            Temperature in Kelvins.

        """
        temp_si = temp_star*self.epsilon_kb
        return temp_si

    def real_energy(self, nrg_star):
        """
        Convert non-dim LJ energy to J.

        Parameters
        -----------
        nrg_star : float
            Non-dimensional energy.

        Returns
        --------
        nrg_si : real
            Energy in J.

        """
        nrg_si = nrg_star*self.epsilon
        return nrg_si

    def real_press(self, press_star):
        """
        Convert non-dim LJ pressure to Pa.

        Parameters
        -----------
        press_star : float
            Non-dimensional pressure.

        Returns
        --------
        press_si : float
            Pressure in Pa.

        """
        press_si = press_star / (self.sigma**3 / self.epsilon)
        return press_si

    def real_num_dens(self, nd_star):
        """
        Convert non-dim LJ number density to [atoms/m^3].

        Parameters
        -----------
        nd_star : float
            Non-dimensional number density.

        Returns
        --------
        nd_si : float
            Number density in [atoms/m^3].

        """
        nd_si = nd_star / self.sigma**3
        return nd_si

    def real_time(self, t_star):
        """
        Convert non-dim LJ time to sec.

        Parameters
        -----------
        t_star : float
            Non-dimensional time.

        Returns
        --------
        t_si : float
            Time in sec.

        """
        t_si = t_star / np.sqrt(self.epsilon / (self.m * (self.sigma)**2))
        return t_si

    def real_velocity(self, v_star):
        """
        Convert non-dim LJ velocity to m/s.

        Parameters
        -----------
        v_star : float
            Non-dimensional velocity.

        Returns
        --------
        v_si : float
            Velocity in m/s.

        """
        v_si = v_star / np.sqrt(self.m / self.epsilon)
        return v_si

    def real_force(self, f_star):
        """
        Convert non-dim LJ force to N.

        Parameters
        -----------
        f_star : float
            Non-dimensional force.

        Returns
        --------
        f_si : float
            Force in N.

        """
        f_si = f_star / (self.sigma/self.epsilon)
        return f_si

    def real_surface_tension(self, gamma_star):
        """
        Convert non-dim LJ surface tension to [N/m].

        Parameters
        -----------
        gamma_star : float
            Non-dimensional surface tension.

        Returns
        --------
        gamma_si : float
            Survace tension in [N/m].

        """
        gamma_si = gamma_star / (self.sigma**2 / self.epsilon)
        return gamma_si

    def real_viscosity(self, mu_star):
        """
        Convert non-dim LJ viscosity to [Pa.s].

        Parameters
        -----------
        mu_star : float
            Non-dimensional viscosity.

        Returns
        --------
        mu_si : float
            Viscosity in [Pa.s].

        """
        mu_si = mu_star / (self.sigma**2 / np.sqrt(self.m* self.epsilon))
        return mu_si


class MetalUnits(object):
    """
    Helper class used to convert metal units to non-dimensional LJ units.

    Each method returns a parameter converted to non-dimensional
    LJ units from LAMMPS metal units as defined by
    sigma, epsilon and mass.

    Parameters
    ----------
    sigma : float, [default=3.405]
        :math:`\sigma` parameter for 12-6 LJ potential in angstroms.
    epsilon_kb : float, [default=119.8]
        :math:`\epsilon / k_b`, where :math:`\epsilon` is energy parameter
        for 12-6 LJ potential. in K.
    m : float, [default=39.95]
        Atomic mass in AMU.

    Attributes
    ----------
    sigma : float
        :math:`\sigma` parameter for 12-6 LJ potential in angstrom.
    epsilon_kb : float
        :math:`\epsilon / k_b`, where :math:`\epsilon` is energy
        parameter for 12-6 LJ potential in K.
    epsilon : float
        epsilon parameter from 12-6 LJ potential function.
    m : float
        Atomic mass in AMU.

    Methods
    --------
    LJ_length(len_metal)
        Convert length to non-dim LJ units.
    LJ_temp(temp_metal)
        Convert temperature to non-dim LJ units.
    LJ_time(t_metal)
        Convert time to non-dim LJ units.
    LJ_energy(nrg_metal)
        Convert energy to non-dim LJ units.
    LJ_num_dens(nd_metal)
        Convert number density to non-dim LJ units.
    LJ_press(press_metal)
        Convert pressure to non-dim LJ units.
    LJ_velocity(v_metal)
        Convert velocity to non-dim LJ units.
    LJ_surface_tension(gamma_metal)
        Convert surface tension to non-dim LJ units.
    LJ_viscosity(mu_metal)
        Convert viscosity to non-dim LJ units.

    References
    ------------
    http://lammps.sandia.gov/doc/units.html

    """

    def __init__(self, sigma=3.405, epsilon_kb=119.8, m=39.95):
        """


        """
        self.sigma = sigma                              # [\angstrom]
        self.epsilon_kb = epsilon_kb                    # [K]
        self.m = m
        self.epsilon = self.epsilon_kb * KBOLTZ['metal']   # [J]

    def LJ_length(self, len_metal):
        """
        Convert metal length to non-dim LJ units.

        Parameters
        -----------
        len_metal : float
            Metal units length.

        Returns
        --------
        len_lj : float
             Non-dim LJ length.

        """
        len_lj = len_metal/self.sigma
        return len_lj

    def LJ_temp(self, temp_metal):
        """
        Convert metal temperature [K] to non-dim LJ units.

        Parameters
        -----------
        temp_metal : float
            Metal units temperature [K].

        Returns
        --------
        temp_lj : float
             Non-dim LJ temperature.

        """
        temp_lj = temp_metal/self.epsilon_kb
        return temp_lj

    def LJ_energy(self, nrg_metal):
        """
        Convert metal energy to non-dim LJ units.

        Parameters
        -----------
        nrg_metal : float
            Metal units energy.

        Returns
        --------
        nrg_lj : float
             Non-dim LJ energy.

        """
        nrg_lj = nrg_metal/self.epsilon
        return nrg_lj

    def LJ_num_dens(self, nd_metal):
        """
        Convert metal units number density to non-dim LJ units.

        Parameters
        -----------
        nd_metal : float
            Metal units number density.

        Returns
        --------
        nd_lj : float
             Non-dim LJ number density.

        """
        nd_lj = nd_metal * self.sigma**3
        return nd_lj

    def LJ_press(self, press_metal):
        """
        Convert metal units pressure to non-dim LJ units.

        Parameters
        -----------
        press_metal : float
            Metal units pressure.

        Returns
        --------
        press_lj : float
             Non-dim LJ pressure.

        """
        press_lj = press_metal*1e5 *((self.sigma*1e-10)**3 / (self.epsilon*eV))
        return press_lj

    def LJ_time(self, t_metal):
        """
        Convert metal units time to non-dim LJ units.

        Parameters
        -----------
        t_metal : float
            Metal units time.

        Returns
        --------
        t_lj : float
             Non-dim LJ time.

        """
        t_lj = t_metal*1e-12 * np.sqrt(self.epsilon*eV /
                                        (self.m*AMU * (self.sigma*1e-10)**2))
        return t_lj

    def LJ_velocity(self, v_metal):
        """
        Convert metal units velocity to non-dim LJ units.

        Parameters
        -----------
        v_metal : float
            Metal units velocity.
        Returns
        --------
        v_lj : float
             Non-dim LJ velocity.

        """
        v_lj = v_metal*(1e-10/1e-12)*np.sqrt(self.m*AMU / (eV*self.epsilon))
        return v_lj

    def LJ_force(self, f_metal):
        """
        Convert metal units force to non-dim LJ units.

        Parameters
        -----------
        f_metal : float
            Metal units force.
        Returns
        --------
        f_lj : float
             Non-dim LJ force.

        """
        f_lj = f_metal * (self.sigma/self.epsilon)
        return f_lj

    def LJ_surface_tension(self, gamma_metal):
        """
        Convert metal units surface tension to non-dim LJ units.

        Parameters
        -----------
        gamma_metal : float
            Metal units surface tension.

        Returns
        -------
        gamma_lj : float
            Non-dim LJ surface tension.

        """
        gamma_lj = gamma_metal * (self.sigma**2 / self.epsilon)



    def LJ_viscosity(self, mu_metal):
        """
        Convert metal units viscosity to non-dim LJ units.

        Parameters
        -----------
        mu_metal : float
            Metal units viscosity.
        Returns
        --------
        mu_lj : float
             Non-dim LJ viscosity.

        """
        mu_lj = mu_metal*0.1 / ((self.sigma*1e-10)**2 /
                                        np.sqrt(self.m*AMU* self.epsilon*eV))
        return mu_lj


 #=============================================================================


def main(argv=None):
    pass

if __name__ == "__main__":
    sys.exit(main())
