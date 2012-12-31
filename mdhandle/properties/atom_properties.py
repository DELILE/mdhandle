# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Calculates atomwise properties based on existing simulation results.

Returns atom-wise results for property calculation

**Credits**:

* Inspiration for virial and stress calc from LAMMPS source.

"""

import sys
import abc

import numpy as np

from mdhandle import units
from mdhandle.lammps import lammps
from mdhandle.logger import Logger
import mdhandle.vectors as vectors
import mdhandle.properties.cell_decomposition as cell_decomposition
from mdhandle.settings import PROP_SIZE, DTYPES

# TODO: Change AtomCalc.add_new_calc(...) so not manually monkeypatching
#       new calculation implemenation.
# TODO: Add additional class like NonLocal to handle other unit types
#       - or add unit factors like in mdhandle.properties.snap_properties


# -----------------------------------------------------------------------------

logger = Logger()

# -----------------------------------------------------------------------------


class AtomCalc(object):
    """
    Abstract base class for atomwise calculations.

    Parameters
    ----------
    snap : :class:`mdhandle.snap.Snap`
        Snap object to be processed.
    calcs : iterable
        List of calculations to perform.
        Each sub-class must supply its own calculation functions which
        are run either directly,  or via the :meth:`run` method.

        The list of available calculations is in ``_all_calcs``.
        If ``calcs is None``, all implemented calculations are done.
        [default=None]
    wrap : boolean
        If ``True``, atom coordinates are wrapped back to primary
        simulation domain, subject to PBC.
        [default=True]

    Attributes
    -----------
    snap : :class:`mdhandle.snap.Snap`
        Must be ``'atoms'`` dataset.
    xyz : :class:`numpy.ndarray` 
        Atomwise positions within :attr:`AtomCalc.snap`.
    scalars : dict
        Dictionary of scalars and associated metadata produced by calculations.
    vectors : dict
        Dictionary of vectors and associated metadata produced by calculations.
    symm_tensors : dict
        Dictionary of symmetric tensors and associated metadata 
        produced by calculations.
    tensors : dict
        Dictionary of tensors and associated metadata produced by calculations.
    calcs : dict
        Stores input for calculations to run.  Stages calculation names until
        processed by``setup_calculations()``.
    to_calc : dict
        Dictionary of calculations performed by :meth:`AtomCalc:run`.
    results : dict
        Dictionary used for storing result of calculations by name.

    _all_calcs : dict 
        Dictionary sotring metadata of all implemented calculations.

    Methods
    --------
    add_new_calc(name, dtype='float', collection_type='scalars')
        Adds metadata for new calculation.
        Calculation function itself must be monkey patched into
        :class:`AtomCalc` and have a name in style ``_calc_NAME()``
    setup_calculations()
        Runs configurations required before calling :meth:`AtomCalc.run`.
    run()
        Runs all calculations stored in :attr:`AtomCalc.to_calc`.
        Must be overridden in sub-classes of :class:`AtomCalc`.

    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, snap, calcs=None, wrap=True):
        self.snap = snap
        self.snap.gather_data()

        if snap.meta['grid_type'] != 'atoms':
            logger.error("Snap dataset type must be 'atoms' not %s"
                                                    % snap.meta['grid_type'])

        self.xyz = snap.get_vector(snap.vectors['xyz']['comps'])

        if wrap is True:
            self.xyz = snap.sim_cell.wrap_to_central(self.xyz)

        # Implemented calculations are added via add_new_calc()
        self._all_calcs = {}

        self.calcs = calcs

        self.scalars = {}
        self.vectors = {}
        self.tensors = {}
        self.symm_tensors = {}
        self.to_calc = {}
        
        self.types = {}
        self.dtype = {}
        self.locations = {}
        
        self.results = {}

        # setup_calculations(...) not called in __init__ for AtomCalc
        # because relies on subclass implementation

    def add_new_calc(self, name, dtype='float', collection_type='scalars',
                                                location='Node'):
        """
        Adds information about the output of a new calculation type beyond the
        default set.

        Implementation of new calculation is done by manually monkeypatching
        function with name in style "_calc_NAME" into existing object.

        Parameters
        -----------
        name : string
            Name of new calculation.
        dtype : string, {'int' | 'float'}, [default='float']
            Data type for calculation.
        collection_type : string, [default='scalars']
            Defines the shape of the output (i.e. vectors: Nx3)
            {'scalars' | 'vectors' | 'tensors' | 'symm_tensors'}
        location : string, [default='Node']
            For atomwise properties as calculated by subclasses of
            :class:`AtomCalc`, only `'Node'` is relevant.

        """
        dtype = dtype.lower()
        collection_type = collection_type.lower()

        if dtype not in DTYPES:
            logger.warning("dtype must be 'int' or 'float', %s given" % dtype)
            return

        if collection_type not in PROP_SIZE:
            logger.warning("collection_type %s does not exist."
                                                            % collection_type)
            return

        if location != 'Node':
            logger.warning("Location must be set to 'Node'")
            return

        self.calcs.append(name)
        self._all_calcs[name] = {'dtype': dtype, 'location': location,
                                  'type': collection_type}
        self.setup_calculations()

    def setup_calculations(self):
        """
        Perform initialization necessary for running calculations.

        Must also be called if any additional new calculations are
        added at runtime.

        Warning: New results dict is created each time this method is called.
        Information from any previous runs is lost.

        """
        # Logic for default calcs == None --> selecting all possible calcs
        if self.calcs is None:
            self.calcs = [i for i in self._all_calcs]
            to_loop = self._all_calcs
        else:
            to_loop = self.calcs

        # Assining items to be calculated.
        for i in to_loop:
            self.to_calc[i] = self._all_calcs[i]

            if i['type'] == 'scalars':
                self.scalars[i] = self._all_calcs[i]
            elif i['type'] == 'vectors':
                self.vectors[i] = self._all_calcs[i]
            elif i['type'] == 'tensors':
                self.tensors[i] = self._all_calcs[i]
            elif i['type'] == 'symm_tensors':
                self.symm_tensors[i] = self._all_calcs[i]

        # Human readable.
        self.dtypes = {}
        for i in self.to_calc:
            self.dtypes[i] = self.to_calc[i]['dtype']

        # Human readable.
        self.locations = {}
        for i in self.to_calc:
            self.locations = self.to_calc[i]['location']

        # Human readable.
        self.types = {}
        for i in self.to_calc:
            self.types = self.to_calc[i]['type']

        self.results = {}

    @abc.abstractmethod
    def run(self):
        """
        Abstract method for running all calculations in 
        :attr:`AtomCalc.to_calc`.
        
        Must be overridden in sub-classes.
        """
        pass


# TODO: Add ability to handle systems w more than one cutoff.
# TODO: Outer product method for getting virial and stress tensor.
class NonLocal(AtomCalc):
    """
    Class for calculating atomwise properties which are dependent on infomation
    gatherered from all atoms within the snap.
    (e.g. potential energy, coordination number.)

    Calculation of interatomic distance (``NonLocal._calc_diff``) assumes
    local effect so that uses only nearest neighbours.  The resulting
    interatomic distance is limited to the cutoff radius and is modified
    to conform with the minimum image convention.

    Due to high cost of calculations - should only be run using ``run()``
    rather than directly calling individual calculation functions
    (e.g. ``_calc_epair()``)

    .. note::
       Interatomic interaction is calculated using 12-6 Lennard-Jones
       potential and therefore must have units style 'lj' (non-dimensional
       Lennard-Jones units.).

    Parameters
    -----------
    snap : :class:`mdhandle.snap.Snap`
        Snap object to be processed.  Must have units style, 'lj'
    cutoff : float
        Cutoff radius for interatomic interaction.
    calcs : iterable, [default=None]
        List of calculations to perform.

        Implemented calculations available:
        ``['coord', 'epair', 'virial', 'force', 'stress']``.

        If ``calcs is None``, all implemented calculations are done.

    wrap : boolean, [default=True]
        If ``True``, atom coordinates are wrapped back to primary
        simulation domain, subject to PBC.

    Attributes
    ------------
    cutoff : float 
        Cutoff for LJ 12-6 interatomic interaction.
    binning : :class:`mdhandle.properties.cell_decomposition.SparseBinning`
        Cell decomposition of atoms for efficient interaction calculation.
    results : dict
        All results from calculations in :attr:`NonLocal.to_calc`
        are stored in `dict` organized by name.  Results also stored in
        attributes of :class:`NonLocal`.

    Methods
    --------
    setup_calculations()
        Runs configuration to setup for calculation.
    run()
        Runs all calculations in :attr:`NonLocal.to_calc`.
    _calc_epair(atom_id, diff, diff_mag, any_there)
        Calculates interatomic potential energy via LJ 12-6.
    _calc_diff(atom_id, diff, diff_mag, any_there)
        Calculates the distance (vector and magnitude) to an atom.
    _calc_coord(atom_id, diff, diff_mag, any_there)
        Calculates coordination number up to ``cutoff``.
    _calc_force(atom_id, diff, diff_mag, any_there)
        Calculates total interatomic force on an atom via LJ 12-6.
    _calc_virial(atom_id, diff, diff_mag, any_there)
        Calculates virial on an atom via LJ 12-6.
    _calc_stress(atom_id, diff, diff_mag, any_there)
        Calculates atomwise stress tensor for an atom via LJ 12-6.

    """

    def __init__(self, snap, cutoff, calcs=None, wrap=True):
        super(NonLocal, self).__init__(snap, calcs, wrap)

        if self.snap.meta['units'] != 'lj':
            logger.error("Units must be 'lj' for NonLocal, currently: %s"
                                                     % self.snap.meta['units'])

        self.cutoff = cutoff

        # Uses cell decomposition to speed up calculation
        self.binning = cell_decomposition.SparseBinning(self.snap,
                                        self.cutoff*np.ones(3, dtype=np.float))

        # Default set of calculations for NonLocal
        self.add_new_calc('coord', dtype='int')
        self.add_new_calc('epair')
        self.add_new_calc('virial', collection_type='symm_tensors')
        self.add_new_calc('force', collection_type='vectors')
        self.add_new_calc('stress', collection_type='symm_tensors')

        # Must be run before proceeding with run()
        self.setup_calculations()

    def setup_calculations(self):
        """
        Perform initialization necessary for running calculations.

        Must also be called if any additional new calculations are
        added at runtime.

        .. warning::
           New results arrays are created each time this method is called.
           Information from any previous runs is lost.

        """
        super(NonLocal, self).setup_calculations()

        # Different than other AtomCalc subclasses because
        # cannot calculate vectorially so need destination arrays
        for i in self.to_calc:
            self.results[i] = np.empty(
                                (self.xyz.shape[0], PROP_SIZE[self.types[i]]),
                                dtype=self.dtypes[i]
                                      )

    def _calc_epair(self, atom_id, diff, diff_mag, any_there):
        """
        Calculate the kinetic energy of atom atom_id.

        Implementation is limited to monotomic system.  Additional terms
        needed for bonded contributions.

        Also not applicable to long-range forces (e.g. electrostatics), which
        cannot be calculated by direct summation under the minimum image
        convention.

        Parameters
        ----------
        atom_id : int
            Atom index in snap.
        diff : :class:`numpy.ndarray` (Nx3, 'float')
            Vector distance between ``atom_id`` and all other atoms.
            Distance is adjusted to conform with the minimum image conv.
        diff_mag : :class:`numpy.ndarray` (Nx1, 'float' )
            Distance magnitude between ``atom_id`` and all other atoms.
            Distance is adjusted to conform with the minimum image conv.
        any_there : boolean
            If ``True``, at least one atom is within the cutoff distance.

        Returns
        -------
        epair : float
            Scalar value for total potential energy of atom ``atom_id``.

        """
        if any_there:
            epair = (2*(diff_mag**(-12) - diff_mag**(-6))).sum()
        else:
            epair = 0.0
        return epair

    def _calc_coord(self, atom_id, diff, diff_mag, any_there):
        """
        Calculate the coordination number of atom ``atom_id``.

        As implemented coordination number produced is the total coordination
        number due to all atoms in the system regardless of type.

        Parameters
        ----------
        atom_id : int
            Atom index in snap.
        diff : :class:`numpy.ndarray` (Nx3, 'float')
            Vector distance between ``atom_id`` and all other atoms.
            Distance is adjusted to conform with the minimum image conv.
        diff_mag : :class:`numpy.ndarray` (Nx1, 'float' )
            Distance magnitude between ``atom_id`` and all other atoms.
            Distance is adjusted to conform with the minimum image conv.
        any_there : boolean
            If ``True``, at least one atom is within the cutoff distance.

        Returns
        -------
        coord : float
            Scalar value for coordination number of atom ``atom_id``.

        """
        if any_there:
            coord = diff.shape[0]
        else:
            coord = 0
        return coord

    def _calc_force(self, atom_id, diff, diff_mag, any_there):
        """
        Calculate the net force due to 12-6 Lennard-Jones
        interaction on atom ``atom_id``.

        Implementation is limited to monotomic system.  Additional terms
        needed for bonded contributions.

        Calculation algorithm not applicable to long-range forces (e.g.
        electrostatics), which cannot be calculated by direct summation under
        the minimum image convention.

        Parameters
        ----------
        atom_id : int
            Atom index in snap.
        diff : :class:`numpy.ndarray` (Nx3, 'float')
            Vector distance between ``atom_id`` and all other atoms.
            Distance is adjusted to conform with the minimum image conv.
        diff_mag : :class:`numpy.ndarray` (Nx1, 'float' )
            Distance magnitude between ``atom_id`` and all other atoms.
            Distance is adjusted to conform with the minimum image conv.
        any_there : boolean
            If ``True``, at least one atom is within the cutoff distance.

        Returns
        -------
        force : :class:`numpy.ndarray`, (1x3)
             Net force on atom ``atom_id``.

        """
        if any_there:
            force = np.dot((2*diff_mag**(-14) - diff_mag**(-8)), -24*diff)
        else:
            force = np.array([0., 0., 0.])
        return force

    def _calc_virial(self, atom_id, diff, diff_mag, any_there):
        """
        Calculate the virial on atom ``atom_id`` based on LJ 12-6
        interaction.

        Implementation is limited to monotomic system.  Additional terms
        needed for bonded contributions.

        Calculation algorithm not applicable to long-range forces (e.g.
        electrostatics), which cannot be calculated by direct summation under
        the minimum image convention.

        Parameters
        ----------
        atom_id : int
            Atom index in snap.
        diff : :class:`numpy.ndarray` (Nx3, 'float')
            Vector distance between ``atom_id`` and all other atoms.
            Distance is adjusted to conform with the minimum image conv.
        diff_mag : :class:`numpy.ndarray` (Nx1, 'float' )
            Distance magnitude between ``atom_id`` and all other atoms.
            Distance is adjusted to conform with the minimum image conv.
        any_there : boolean
            If ``True``, at least one atom is within the cutoff distance.

        Returns
        -------
        virial : :class:`numpy.ndarray`, (1x6)
             Virial (symmetric tensor) of atom ``atom_id``.

        """
        if any_there:
            force_per = ( (2*diff_mag**(-14) - diff_mag**(-8))*(-24*diff).T ).T
            v0 = -0.5*np.dot( force_per[:,0], diff[:,0] )
            v1 = -0.5*np.dot( force_per[:,0], diff[:,1] )
            v2 = -0.5*np.dot( force_per[:,0], diff[:,2] )
            v3 = -0.5*np.dot( force_per[:,1], diff[:,1] )
            v4 = -0.5*np.dot( force_per[:,1], diff[:,2] )
            v5 = -0.5*np.dot( force_per[:,2], diff[:,2] )

            virial = np.array( (v0, v1, v2, v3, v4, v5))
        else:
            virial = np.array( (0., 0., 0., 0., 0., 0.))
        return virial

    def _calc_stress(self, atom_id, diff, diff_mag, any_there):
        """
        Calculates the per-atom stress tensor as defined in the LAMMPS
        MD software package based on 12-6 Lennard-Jones interaction.

        In this case the per-atom stress is ``stress*per_atom_volume``.
        Therefore, to calculate the scalar pressure from the trace the result
        should be divided by (``dim*volume``) to obtain the pressure, where
        ``dim`` is the dimensionality of the system.

        Implementation is limited to monotomic system.  Additional terms
        needed for bonded contributions or constraints on atom motion.

        Calculation algorithm not applicable to long-range forces (e.g.
        electrostatics), which cannot be calculated by direct summation under
        the minimum image convention.

        **See**: http://lammps.sandia.gov/doc/compute_stress_atom.html

        Parameters
        ----------
        atom_id : int
            Atom index in snap.
        diff : :class:`numpy.ndarray` (Nx3, 'float')
            Vector distance between ``atom_id`` and all other atoms.
            Distance is adjusted to conform with the minimum image conv.
        diff_mag : :class:`numpy.ndarray` (Nx1, 'float' )
            Distance magnitude between ``atom_id`` and all other atoms.
            Distance is adjusted to conform with the minimum image conv.
        any_there : boolean
            If ``True``, at least one atom is within the cutoff distance.

        Returns
        -------
        stress : :class:`numpy.ndarray`, (1x6)
            Atomwise stress (symmetric tensor) of atom ``atom_id``.

        """
        if any_there:

            m = self.mass[atom_id]
            v = self.velocity[atom_id]

            force_per = ( (2*diff_mag**(-14)-diff_mag**(-8))*(-24*diff).T ).T

            s0 = -1*( m*v[0]*v[0]- 0.5*np.dot(force_per[:,0], diff[:,0]) )
            s1 = -1*( m*v[0]*v[1]- 0.5*np.dot(force_per[:,0], diff[:,1]) )
            s2 = -1*( m*v[0]*v[2]- 0.5*np.dot(force_per[:,0], diff[:,2]) )
            s3 = -1*( m*v[1]*v[1]- 0.5*np.dot(force_per[:,1], diff[:,1]) )
            s4 = -1*( m*v[1]*v[2]- 0.5*np.dot(force_per[:,1], diff[:,2]) )
            s5 = -1*( m*v[2]*v[2]- 0.5*np.dot(force_per[:,2], diff[:,2]) )

            stress = np.array( (s0, s1, s2, s3, s4, s5))
        else:
            stress = np.array(( 0., 0., 0., 0., 0., 0.))
        return stress

    def _calc_diff(self, atom_id):
        """
        Calculates the distance (magnitude and vector) between ``atom_id`` and
        other atoms in the snap.

        Calculated distance is limited by :attr:`NonLocalLammps.cutoff`.

        Resuts are stored in :attr:`NonLocal.diff` and
        :attr:`NonLocal.diff_mag`.

        `NxN` difference arrays are not calculated or stored due to large
        memory requirements.

        Parameters
        -----------
        atom_id : int
            Index number of central atom in distance calculation.
            (i.e. all distances are relative to its position.)

        Returns
        -------
        diff : :class:`numpy.ndarray`, (Nx3)
            Vector distance betweeen atom ``atom_id`` and all other atoms
            within the cutoff distance.
        diff_mag : :class:`numpy.ndarray`, (N,1)
            Scalar distance betweeen atom ``atom_id`` and all other atoms
            within the cutoff distance.
        any_there : boolean
            If ``True``, there is at least one atom within the cutoff radius
            of atom ``atom_id``.

        """
        r0 = self.xyz[atom_id]

        bin_coll = self.binning.get_surrounding(r0)
        bin_coll.sort()

        # Removing itself from nearest neighbours.
        try:
            bin_coll.remove(atom_id)
        except ValueError:
            # Can be caused by unhandled edge cases at boundaries.
            logger.warning('Atom %s not in own surrounding bins' % atom_id )

        diff = self.xyz[bin_coll] - r0
        diff = self.snap.sim_cell.mic_distance(diff)

        diff_mag = vectors.magnitude(diff)

        # Eliminating atoms farther than cutoff.
        mask = diff_mag < self.cutoff
        diff = np.compress(mask, diff, axis=0)
        diff_mag = np.compress(mask, diff_mag)

        # There are no atoms within the cutoff distance.
        if diff.shape[0] == 0:
            any_there = False
        # There is at least one atom within the cutoff distance.
        else:
            any_there = True

        return diff, diff_mag, any_there

    def run(self, mass_name='mass', velocity_name='velocity' ):
        """
        Run all calculations as defined by :attr:`NonLocal.to_calc`.

        Makes use of cell decomposition to reduce work in looking for nearest
        neighbours.

        Parameters
        ----------
        mass_name : string, [default='mass']
            Name for atomwise mass vector.
        velocity_name : string, [default='velocity']
            Name for atomwise velocity vector collection.

        """
        if hasattr(self,'velocity') is not True:
            vel_comps = self.snap.vectors[velocity_name]['comps']
            self.velocity = self.snap.get_vector(vel_comps)
        if hasattr(self,'mass') is not True:
            self.mass = self.snap.get_scalar(mass_name)

        for atom_id in xrange(0, self.xyz.shape[0]):
            diff, diff_mag, any_there = self._calc_diff(atom_id)

            for i in self.to_calc:
                self.results[i][atom_id] = \
                  getattr(self, '_calc_'+i)(atom_id, diff, diff_mag, any_there)

        # Making each results array an attribute
        for i in self.to_calc:
            setattr(self, i, self.results[i])


# -----------------------------------------------------------------------------


# TODO: Use overt Python calls to LAMMPS library rather than input file.
#           - e.g. calculations defined by NonLocalLammps method not input file
# TODO: Add call to generate LAMMPS data file from snap object at run time.
# TODO: Extend to other types of properties fix, atom properties (velocity)


class NonLocalLammps(AtomCalc):
    """
    Class used for calculating non-local properties using a Python wrapped
    LAMMPS library to perform the calculation.

    Non-local properties are atomwise properties like stress and potential
    energy which are a function of all other atoms in the simulation.

    Calculation is restricted to properties which are produced by LAMMPS
    compute commands.

    LAMMPS is run for a single timestep to perform the calculation.  Therefore,
    this operation is intended for use if a calculation was not done at
    the original simulation runtime.

    Current implementation uses a modified serial version of LAMMPS built
    as a library.

    Parameters
    -----------
    snap : :class:`mdhandle.snap.Snap`
        Snap object used for calculation.
    input_fn : string
        LAMMPS input file for controlling LAMMPS library.
        LAMMPS input file is written separately as per format specified
        by the LAMMPS user manual.
        Input file should contain reference to LAMMPS data file which
        contains the atom positions and velocities.
    calcs : iterable, [default=None ]
        List of property names to calculate.
        If ``calcs is None``, then all implemented calculations are done.

        Implemented calculations: ``['coord', 'epair', 'force', 'stress']``
        which must match the name of LAMMPS compute.

        To add additional calculaion types must add additional item
        using :meth:`NonLocalLammps.add_new_calc` method.
    lammps_log : string, [default='lammps_mdhandle.log']
        Log file name for LAMMPS run.

    Attributes
    -----------
    input_fn : string
        LAMMPS input file name used to setup calculations.
    lammps_log : string
        LAMMPS log name.
    cutoff : float
        Cutoff associated with interatomic interaction.
    results : dict
        All results from calculations in :attr:`NonLocalLammps.to_calc`
        are stored in `dict` organized by name.  Results also stored in
        attributes of :class:`NonLocalLammps`.

    Methods
    --------
    add_new_calc(name, dtype='float', collection_type='scalars') 
        Add metadata for new calculation function.
        Overrides completely :meth:`AtomCalc.add_new_calc`.
    run()
        Runs all calculations stored in :attr:`NonLocalLammps.to_calc`.

    """

    def __init__(self, snap, input_fn, cutoff, calcs=None, wrap=True,
                                       lammps_log='lammps_mdhandle.log'):
        super(NonLocalLammps, self).__init__(snap, calcs, wrap)
        self.input_fn = input_fn
        self.lammps_log = lammps_log

        # Not currently used within NonLocalLammps obj
        self.cutoff = cutoff

        # Default calculations
        self.add_new_calc('coord', dtype='int')
        self.add_new_calc('epair')
        self.add_new_calc('force', collection_type='vectors')
        self.add_new_calc('stress', collection_type='symm_tensors')

        # must be called before NonLocalLammps.run(...)
        self.setup_calculations()

    def add_new_calc(self, name, dtype='float', collection_type='scalars',
                                                location='Node'):
        """
        Adds information about the output of a new calculation type beyond the
        default set: ``['coord', 'epair', 'force', 'stress']``

        For :class:`NonLocalLammps`, new calculations are defined within the
        LAMMPS data file rather than via a NonLocalLammps method.  Therefore,
        the name must match the name of the LAMMPS compute.

        Parameters
        -----------
        name : string
            Name of new calculation.  Must be the name of the LAMMPS compute.
        dtype : string, {'int' | 'float'}, [default='float']
            Data type for calculation.
        collection_type : string, [default='scalars']
            Defines the shape of the output (i.e. vectors: Nx3)
            {'scalars' | 'vectors' | 'tensors' | 'symm_tensors'}
        location : string, [default='Node']
            For atomwise properties as calculated by subclasses of
            :class:`AtomCalc`, only `'Node'` is relevant.

        """
        super(NonLocalLammps, self).add_new_calc(name, dtype, 
                                                  collection_type, location)
        logger.warning('Add new calculation %s to LAMMPS input file' % name)

    def run(self):
        """
        Launches LAMMPS based calculation defined by requested calculations and
        the LAMMPS input file (:attr:`NonLocalLammps.input_fn`).

        Result of calculations stored in :attr:`NonLocalLammps.results` and
        :class:`NonLocalLammps` attribute with name corresponding
        to the calculation name stored in :attr:`NonLocalLammps.calcs`.

        """

        # Launching LAMMPS
        lmp = lammps()
        lmp.command(self.lammps_log)

        # Reading LAMMPS input file.
        lines = open(self.input_fn, 'r').readlines()
        for line in lines:
            lmp.command(line)

        # HARDCODED: Running LAMMPS for a single timestep
        lmp.command('run 0')

        for name in self.calcs:
            if self.types[name] == 'scalars':
                output_type = 0
                tmp = lmp.lazy_extract_compute(name, output_type)
                self.results[name] = np.ctypeslib.as_array(tmp).\
                                                astype(self.dtypes[name])

            # Non-scalars required special handling of 2-dim shapes
            else:
                if self.types[name] == 'vectors':
                    output_type = 1
                elif self.types[name] == 'tensors':
                    output_type = 2
                elif self.types[name] == 'symm_tensors':
                    output_type = 3
                # Pathalogical case - cannot find matching type. Skip
                else:
                    logger.warning('Could not find type %s' % self.types[name])
                    continue

                tmp = lmp.lazy_extract_compute(name, output_type)

                # Saving output to NonLocalLammps.results.
                self.results[name] = np.ctypeslib.as_array(tmp).\
                                        reshape((self.snap.meta['num_atoms'], 
                                                            PROP_SIZE[name]))

            # Putting each result in an attribute as well
            setattr(self, name, self.results[name])

        # Closing LAMMPS session
        lmp.close()


# -----------------------------------------------------------------------------


class Local(AtomCalc):
    """
    Calculates atomwise quantities which are only dependent on the
    properties of each individual atom (e.g. atomwise kinetic energy).

    Parameters
    -----------
    snap : :class:`mdhandle.snap.Snap`
        Snap object to be processed.
    calcs : iterable
        List of local atomwise property calculations to perform.
        Implemented calculations: ``['ke', ]``
        If ``calcs is None``, all available calculations are done.
        [default=None]
    wrap : boolean
        If ``True``, atom coordinates are wrapped back into the primary
        simullation domain, subject to PBC.

    Attributes
    -----------
    results : dict
        Results from calculations stored in 
        :attr:`Local.results` and in attributes of :class:`Local`.

    Methods
    ---------
    run()
        Runs all calculations stored in :attr:`Local.to_calc`.

    _calc_ke()
        Calculates atomiwise kinetic energy.

    """

    def __init__(self, snap, calcs=None, wrap=True):
        super(Local, self).__init__(snap, calcs, wrap)

        self.add_new_calc('ke')

        # Must be run before calling Local.run(...)
        self.setup_calculations()

    def _calc_ke(self):
        """
        Calculates atomwise kinetic energy from atomwise mass and velocity.

        Returns
        -------
            Atomwise kinetic energy.
            Result is also stored in :attr:`Local.ke` and
            :attr:`Local.results['ke']`

        """
        u = self.snap.meta['units']
        units_factor = units.MASS[u]*units.VELOCITY[u]**2  / units.ENERGY[u]

        mass_name = self._all_calcs['ke']['name'][0]
        velocity_name = self._all_calcs['ke']['name'][1]

        mass = self.snap.get_scalar(mass_name)
        vel = self.snap.get_vector(velocity_name)
        self.ke = units_factor * (0.5 * mass * 
                       vectors.sq_magnitude(vel)).astype(self.to_calc['dtype'])
        self.results['ke'] = self.ke
        return self.ke

    def run(self):
        """
        Runs all atomwise property calculations as listed in
        :attr:`Local.to_calc`.

        """
        for i in self.to_calc:
            getattr('_calc_'+i)()


# =============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
