# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Calculate average quantities on grid from atomiwise values.
(e.g. local streaming velocity, local temperature)

**Credits**: Inspiration for binning algorithm from binning module in ``mol_mod`` module in ``MD-Tracks`` package


"""

import sys
import abc
from math import fsum
import inspect

import numpy as np

import mdhandle.utilities as util
import mdhandle.vectors as vectors
import mdhandle.units as units
from mdhandle.logger import Logger
from mdhandle import settings

# -----------------------------------------------------------------------------

# TODO: Cylindrical or spherical grids.
# TODO: Grid manipulations (refinement)
# TODO: Rectalinear grid with variation along coordinate directions.
# TODO: Grids that are not aligned to XYZ coordinate directions.
# TODO: Change Grid.add_new_calc(...) so not manually monkeypatching
#       new calculation implemenation.
# TODO: Edit docs to ensure that user knows how to sub-class Grid
#       (i.e. what is required of mapping functions), and how to implement
#       a calculation function.
# TODO: Refactor abstract class so more elegant.
#       - Currentl abstract method: mapping_func, and also require
#         definition of

# -----------------------------------------------------------------------------

logger = Logger()

# -----------------------------------------------------------------------------

def get_implemented_grids():
    """
    Get the names of grid types implemented (i.e. sub-classes of :class:`Grid`
    class).

    Taken from the :attr:`Grid.grid_type` properties in subclasses of
    :class:`Grid`.

    'atoms' is the default grid type for an atomwise dataset.

    Returns
    --------
    implemented_grids : list
        List of strings for implemented grid type.

    """
    # TODO: Handling new grids implemented later by other users.
    current_mod = sys.modules[__name__]
    implemented_grids = ['atoms']
    for name, obj in inspect.getmembers(current_mod, inspect.isclass):
        if issubclass(obj, Grid) and (obj is not Grid):
            implemented_grids.append(obj.grid_type)
    return implemented_grids

# -----------------------------------------------------------------------------
#                       Grid Objects
# -----------------------------------------------------------------------------

class Grid(object):
    """
    Abstract grid class.

    Intent is to provide a useful base class for other grid types to inherit
    as the property calculations are not dependent on the grid type or the
    mapping function.

    Parameters
    ----------
    snap : :class:`mdhandle.snap.Snap`
        Snap for calculating binned properties.
        Current data set in snap should be atomwise data
        (e.g. ``'rawSimResults'``)
    calcs : iterable, [default=None]
        List of calculations to perform on the grid.
        If ``calcs is None``, then all existing calculations are performed.
        Existing calculations :
            ``['count', 'num_density', 'mass', 'mass_density',``
            ``'epair_count', 'coord_count', 'velocity',``
            ``'ke_count', 'force_count', 'temp', 'stress']``
    wrap : boolean, [default=True]
        If ``True``, coordinates in active data set in snap are
        wrapped to the central simulation cell subject to PBC.
    selection : boolean, [default=True]
        If ``True``, selections are applied to snap object when getting
        atom coordinates.
    location : string, [default='Node']
        Grid location of calculated values.
        {'Node' } ``'Cell'`` location not currently implemented.

    Attributes
    -----------
    snap : :class:`mdhandle.snap.Snap`
        Snap object.
    use_selection : boolean
        If ``True``, any selections in :attr:`Grid.snap` are used.
    atom_xyz : :class:`numpy.ndarray`
        Atomic positions for atoms in :attr:`Grid.snap`.
    calcs : iterable
        List of calculation names (strings) to perform.
    to_calc : dict
        List of calculations to do via :meth:`Grid.run_calc` and metadata.
    scalars : dict
        List of scalars produced by grid calcs.
    vectors : dict
        List of vectors produced by grid calcs.
    tensors : dict
        List of tensors produced by grid calcs.
    symm_tensors : dict
        List of symmetric tensors produced by grid calcs.
    dtypes : dict
        dtype of data produced by grid calcs
    locations : dict
        Location of data produced by grid calcs
    types : dict
        Type of output produced by grid calcs
    bins : dict
        Dictionary storing result of sparse atoms-to-grid mapping.
    grid_to_atom : :class:`numpy.ndarray`
        Stores grid cells for each atom in the :attr:`Grid.snap`.
    elem_volume : float
        Abstract property - volume of grid element.
    num_elem : :class:`numpy.ndarray`
        Abstract property - number of grid elements.
    grid_type : string
        Abstract property - Name of grid type (e.g. '3DRegular')
    _KBOLTZ : float
        Boltzmann constant in units consistent with :attr:`Grid.snap`.
    _all_calcs : dict
        List of all implemented grid calculations.
    _run_once : dict
        If values are true - then grid calculation has been run.

    Notes
    ------

    * Each grid calculation (``_calc_NAME()``) stores its result in 
      :class:`Grid` attribute with name matching name in :attr:`Grid.to_calc`.
    * See :meth:`RectGrid._mapping_ngp` or :meth:`RectGrid._mapping_cic`
      for basic requirements of a mapping function.
    * Abstract methods and properties are dependent on the specific and 
      concrete grid implemenation.
     * Grid calculation shouldonly be run once per intance of :class:`Grid`.

    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, snap, calcs=None, wrap=True, use_selection=True,
                                                    location='Node'):

        logger.user_message('Using abstract ABC Grid __init__() function')
        self.snap = snap
        self.snap.gather_data()

        if self.snap.meta['grid_type'] != 'atoms':
            logger.error("Snap dataset must be of type 'atoms'")

        self.calcs = calcs
        self.use_sel = use_selection

        self._KBOLTZ = units.KBOLTZ[self.snap.meta['units']]
        self.atom_xyz = self.snap.get_vector( ('x', 'y', 'z'),
                                              raw_data=(not self.use_sel) )

        if wrap is True:
            self.atom_xyz = snap.sim_cell.wrap_to_central(self.atom_xyz)

        # Setting up for binning calculation
        self.bins = {}
        self.grid_to_atom = None

        # Existing calculations
        self._all_calcs = {}
        self.add_new_calc('count', location=location)
        self.add_new_calc('num_density', location=location)
        self.add_new_calc('mass', location=location,
                                  atomwise_name={'mass': 'mass'})
        self.add_new_calc('mass_density', location=location)
        self.add_new_calc('epair_count', location=location,
                                         atomwise_name={'epair': 'epair'})
        self.add_new_calc('coord_count', location=location,
                                         atomwise_name={'coord': 'coord'})
        self.add_new_calc('velocity', location=location,
                                      collection_type='vectors',
                                      atomwise_name={'velocity': 'velocity',
                                                     'mass': 'mass'})
        self.add_new_calc('ke_count', location=location,
                                      atomwise_name={'ke': 'ke'})
        self.add_new_calc('force_count', location=location,
                                         collection_type='vectors',
                                         atomwise_name={'force': 'force'})
        self.add_new_calc('temp', location=location,
                                  atomwise_name={'velocity': 'velocity',
                                                 'mass': 'mass'})
        self.add_new_calc('temp_uncorr', location=location,
                                         atomwise_name={'velocity': 'velocity',
                                                        'mass': 'mass'})
        self.add_new_calc('stress', location=location,
                                    collection_type='symm_tensors',
                                    atomwise_name={'stress': 'stress'})
        self.to_calc = {}
        self.scalars = {}
        self.vectors = {}
        self.tensors = {}
        self.symm_tensors = {}

        self.types = {}
        self.dtypes = {}
        self.locations = {}
        
        self.results = {}
        self.coverage = False
        self.num_atoms = 0.
        self._run_once = {'run_binning': False,
                          'run_calc': False}

        # Do not call setup_calculations in abstrct Grid class
        # because setup may rely on details of concrete class
        # implementation.  Should call within concrete class __init__().

    @abc.abstractmethod
    def mapping_func(self):
        pass
    
    @abc.abstractproperty
    def grid_type(self):
        pass
    
    @abc.abstractproperty
    def elem_volume(self):
        pass

    @abc.abstractproperty
    def num_elem(self):
        pass


    def set_mapping_function(self, name):
        """

        Parameters
        ----------
        name : string
            Name of mapping function.
            Must match method within grid class
            (e.g. ``GridClass._mapping_NAME()``):

        """
        if ('_mapping_' + name) in dir(self):
            self.mapping_func = getattr(self, '_mapping_' + name)
        else:
            logger.error('Mapping function %s does not exist' % name)

    def add_new_calc(self, calc_name, dtype='float', collection_type='scalars',
                                      location='Node', atomwise_name=None,\
                                      run_setup=False):
        """
        Adds metadata about the output of a new calculation type beyond the
        default set.

        Implementation of new calculation is done by manually monkeypatching
        function with name in style ``_calc_NAME(...)`` into existing object.

        Parameters
        -----------
        calc_name : string
            Name of calculation.
        dtype : string, {'int' | 'float'}, [default='float']
            Data type for calculation.
        collection_type : string, [default='scalars']
            Defines the shape of the output (i.e. vectors: Nx3)
            {'scalars' | 'vectors' | 'tensors' | 'symm_tensors'}
        location : string, [default='Node']
            Grid location of calculated values.
            {'Node'} ``'Cell'`` location not currently implemented.
        atomwise_name : dict, [default=None]
            Name of atomwise property associated with grid calculation to
            allow for variation in the name stored for propety in
            :attr:`Grid.snap`. Each ``dict`` entry is formatted as:
        
            ``'identifying_name':'name_in_snap'``
        run_setup : boolean, [default=False]
            Runs :meth:`Grid.setup_calculations` after new calculationi is
            added.  Not recommended for use during initial setup such as
            in grid __init__.
        
        """
        dtype = dtype.lower()
        collection_type = collection_type.lower()

        if dtype not in settings.DTYPES:
            logger.warning("dtype must be 'int' or 'float', %s given" % dtype)
            return

        if collection_type not in settings.PROP_SIZE:
            logger.warning('collection_type %s does not exist.'
                                                            % collection_type)
            return

        if location != 'Node':
            logger.warning("Location must be set to 'Node'")
            return

        self.calcs.append(calc_name)
        self._all_calcs[calc_name] = {'dtype': dtype, 'location': location,
                                      'type': collection_type,
                                      'atomwise_name': atomwise_name}

        if run_setup is True:
            self.setup_calculations()

    def setup_calculations(self):
        """
        Configures list of calculations to performs and list of scalars
        (:attr:`Grid.scalars`), vectors (:attr:`Grid.vectors`), tensors
        (:attr:`Grid.tensors`) and symmetric tensors
        (:attr:`Grid.symm_tensors`).

        Moves requested calculations from :attr:`Grid.calcs` to
        :attr:`Grid.to_calc` so can be used by :meth:`Grid.run_calc`.

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

            if self.to_calc[i]['type'] == 'scalars':
                self.scalars[i] = self._all_calcs[i]
            elif self.to_calc[i]['type'] == 'vectors':
                self.vectors[i] = self._all_calcs[i]
            elif self.to_calc[i]['type'] == 'tensors':
                self.tensors[i] = self._all_calcs[i]
            elif self.to_calc[i]['type'] == 'symm_tensors':
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

    def change_default_name(self, calc_name, 
                                  identifying_name, new_atomwise_name):
        """
        Changes the default name of atomwise property used in calculation.
        (e.g. ``'epair'`` is default for potential energy, but could be
        ``'c_epair'`` depending on LAMMPS configuration).

        See :attr:`Grid._all_calcs` for defaults.

        Parameters
        -----------
        calc_name : string
            Name of grid calculation.
        identifying_name : string
            Identifying name for atomwise property (i.e. 'velocity').
            This may not correspond with the name stored in 
            :attr:`Grid.snap`.
        new_atomwise_name : string
            New name for related atomwise property.

        """
        self._all_calcs[calc_name]['atomwise_name'][identifying_name] = \
                                                              new_atomwise_name

    # ---------------- RUNNING BINNING ----------------------------------------

    def run_mapping(self, *args, **kwargs):
        """
        Performs atom-to-grid mapping according to
        :meth:`Grid.mapping_func`.

        The result of mapping is stored in :attr:`Grid.bins`.

        .. note::
           This method must be called before :meth:`Grid.run_calc`.

        """
        if self._run_once['run_binning'] is True:
            logger.error('Grid calculation cannot be run more than once.')
        else:
            self._run_once['run_binning'] = True
            self.mapping_func(*args, **kwargs)

    def run_calc(self):
        """
        Systematically calls all grid calculation functions, ``calc_NAME()``
        as listed in :attr:`Grid.to_calc`.

        The result of each calculation is stored in :attr:`Grid.results`
        dictionary and in a :class:`Grid` attribute with with the name
        corresponding to the name in :attr:`Grid.to_calc`.

        """
        if self._run_once['run_calc'] is True:
            logger.error('Grid calculation cannot be run more than once.')
        
        else:
            if self.bins == {}:
                logger.error('Atom-to-grid mapping must be done first.')

            for calc_name in self.to_calc:
                # Skipping calculation if it already exists.
                if hasattr(self, calc_name):
                    continue
                else:
                    # Storing the result in self.results and self.calc_name.
                    self.results[calc_name] = getattr(self, 
                                                      '_calc_' + calc_name)()
                    setattr(self, calc_name, self.results[calc_name])

    # --------------- CALCULATION UTILITY FUNCTIONS --------------------------

    def _requirements(self, *list_reqs):
        """
        Determines if gridded property exists and if not calculates it.

        Used in grid calculations that are dependent on existence of other
        grid calculations.

        Parameters
        ----------
        *list_reqs: strings
            Names of required grid properties.  The name should match the
            grid properties attribute name (e.g. 'epair' for `Grid.epair`).
            Can supply many requirements as comma sparated arguments.
            (e.g. ``'count', 'epair'``).

        """
        for req in list_reqs:
            if not hasattr(self, req):
                self.results[req] = getattr(self, '_calc_' + req)()
                setattr(self, req, self.results[req])
            else:
                continue

    def _reverse_mapping(self, grid_value):
        """
        Performs grid-to-atom mapping of ``grid_value``.

        Uses mapping information stored in :attr:`Grid.grid_to_atoms` when
        atom-to-grid mapping is originally performed.

        Parameters
        -----------
        grid_value : array-like or float
            Gridded data array.  If a constant ``float`` is given, each atom is
            assigned the constant value.

        Returns
        -------
        atom_value : :class:`numpy.ndarray`
            Atomwise array resulting from grid-to-atom mapping of
            ``grid_value``.

        """
        # Defining helper function to get atom values
        # back from grid_value.
        def get_total_weighting(idx, weights):
            """
            Calculates the weighted average of grid values
            from all grid points that an individual atom contributes
            to.

            ``@np.vectorize`` decorator produces a ``numpy.ufunc`` which
            operates over array inputs without looping.

            Parameters
            -----------
            idx : array-like, (Nx1)
                Numpy array containing a list of index tuples
                (i.e. ``(1,1,1)``) for each atom.
            weights : array-like (Nx1)
                Numpy array containing list of weights (i.e. [0.9, 0.6])
                for each atom.

            Returns
            -------
            weighted_avg : :class:`numpy.ndarray` (Nx1)
                Weighted average of `grid_value` for each atom based
                on atomic contributions at all grid points given the
                :meth:`Grid.mapping_func`.

            """
            if idx is not None:
                total = fsum(weights)
                temp = [grid_value[i]*w for (i, w) in zip(idx, weights)]
                                
                prod = np.array(temp).sum(axis = 0)
                return prod / total
            else:
                return np.zeros(prod.shape)

        # end get_total_weighting(...)
        vec_get_total_weighting = np.vectorize(get_total_weighting, 
                                               otypes=[np.object,])
        
        
        tmp = vec_get_total_weighting(self.grid_to_atom[:,0], 
                                      self.grid_to_atom[:,1])

        atom_value = np.zeros((tmp.shape[0], tmp[0].shape[0]))
        for i in range(tmp.shape[0]):
            for j in range(tmp[0].shape[0]):
                atom_value[i,j] = tmp[i,j]

        return atom_value

    def _sum_from_atom_value(self, atom_value):
        """
        Perform summation of atomwise values on the grid for arbitrary
        atomwise data (``atom_value``).

        Utility function used in other grid calculations, ``_calc_mass()``.

        For example, to count the number of atoms at the gird points
        (or grid cells), :meth:`Grid._sum_from_atom_value(1.0)` or to
        calculate the mass of each grid point,
        :meth:`Grid._sum_from_atom_value(atom_masses)`.

        Parameters
        ----------
        atom_value : array-like or float
            Atomwise data object to be summed according to grid weighting at
            each grid point or cell.

        """
        # Add special case for temp, which is active when
        # Flag (to be added) is thrown asking for velocity correction to
        # be used.

        # atom_value is a scalar, constant across all atoms
        if not isinstance(atom_value, np.ndarray):

            gridded_value = np.zeros(self.num_elem)

            for bin0, contents in self.bins.iteritems():
                gridded_value[bin0[0], bin0[1], bin0[2]] = \
                                    atom_value*fsum(contents['weights'])
            return gridded_value

        # atom_value varies over atoms
        else:
            if atom_value.ndim == 1:
                gridded_value = np.zeros(self.num_elem)
            else:
                individual_size = atom_value.shape[1]
                gridded_value = np.zeros( self.num_elem,
                                dtype=(np.float, (individual_size,))  )

            for bin0, contents in self.bins.iteritems():
                gridded_value[bin0[0], bin0[1], bin0[2]] = \
                    (atom_value[list(contents['atoms'])].T *
                                     contents['weights']).T.sum(axis=0)
            return gridded_value

    # --------------- SPECIFIC GRID CALCULATIONS ------------------------------

    def _calc_count(self):
        """
        Calculate the atom count on grid.

        Also calls tests on grid coverage: :meth:`Grid._calc_coverage`.

        Returns
        -------
        count : :class:`numpy.ndarray`
            Number of atoms at each grid point.
            Also stored in :attr:`Grid.results['count']` and :attr:`Grid.count`
            in `Grid.run()`.

        """
        count = self._sum_from_atom_value(1.0)
        return count

    def _num_atoms(self):
        """
        Calculates the total number of atoms on grid.

        .. note::
           The result, stored in :attr:`Grid.num_atoms`, is not a grid array 
           and is not  stored in :attr:`Grid.coverage.results`.

        """
        self._requirements('count')
        self.num_atoms = self.count.sum()

    def _coverage(self):
        """
        Test whether atom count on the grid captures all atoms
        in the snap by comparing to number of atoms in
        :attr:`Grid.snap` (i.e. metadata stored on disk).

        The result is ``True`` if the number of atoms in grid is equal to the
        total number of atoms in the snap.

        .. note::
           The result, stored in :attr:`Grid.coverage`, is not a grid array and
           is not  stored in :attr:`Grid.results`.

        """
        self._requirements('count')
        self._num_atoms()

        # Comparing total number of atoms
        if self.num_atoms == self.snap.meta['num_atoms']:
            self.coverage = True
        else:
            self.coverage = False

    def _calc_num_density(self):
        """
        Calculates the gridded number density.

        Returns
        -------
        num_density : :class:`numpy.ndarray`
            Number density at each grid point.
            Also stored in :attr:`Grid.results['num_density']` and
            :attr:`Grid.num_density` in `Grid.run()`.

        """
        self._requirements('count')
        num_density = self.count / self.elem_volume
        return num_density

    def _calc_mass(self):
        """
        Calculates the mass of at each grid point (or grid cell).

        Returns
        -------
        mass : :class:`numpy.ndarray`
            Mass at each grid point or grid cell.
            Also stored in :attr:`Grid.results['mass']` and :attr:`Grid.mass`
            in `Grid.run()`.

        """
        mass_name = self._all_calcs['mass']['atomwise_name']['mass']
        mass = self._sum_from_atom_value( self.snap.get_scalar(mass_name,
                                                raw_data=(not self.use_sel)) )
        return mass

    def _calc_mass_density(self):
        """
        Calculates mass density at each grid point (or grid cell).

        Returns
        -------
        mass_density : :class:`numpy.ndarray`
            Mass density at on the grid.
            Also stored in :attr:`Grid.results['mass_density']` and
            :attr:`Grid.mass_density` in `Grid.run()`.

        """
        self._requirements('mass')
        mass_density = self.mass / self.elem_volume
        return mass_density

    def _calc_epair_count(self):
        """
        Calculates the potential energy normalized by atom count on the grid.

        Returns
        -------
        epair_count : :class:`numpy.ndarray`
            Specific potential energy normalized by atom count on the grid.
            Also stored in :attr:`Grid.results['epair_count']` and
            :attr:`Grid.epair_count` in `Grid.run()`.

        """
        self._requirements('count')
        epair_name = self._all_calcs['epair']['atomwise_name']['epair']

        epair_count = self._sum_from_atom_value(
                self.snap.get_scalar(epair_name, raw_data=(not self.use_sel))
                                                ) / self.count
        util.clean(epair_count, ~np.isfinite(epair_count), 0.)
        return epair_count

    def _calc_coord_count(self):
        """
        Calculates the coordintion number normalized by atom count on the grid.

        Returns
        -------
        coord_count : :class:`numpy.ndarray`
            Coordination number normalized by atom count on the grid.
            Also stored in :attr:`Grid.results['coord_count']` and
            :attr:`Grid.coord_count` in `Grid.run()`.

        """
        self._requirements('count')
        coord_name = self._all_calcs['coord']['atomwise_name']['coord']

        coord_count = self._sum_from_atom_value(
                  self.snap.get_scalar(coord_name, raw_data=(not self.use_sel))
                                               ) / self.count
        util.clean(coord_count, ~np.isfinite(coord_count), 0.)
        return coord_count

    def _calc_velocity(self):
        """
        Calculates the streaming velocity on the grid.

        Returns
        -------
        velocity : :class:`numpy.ndarray`
            Streaming velocity on the grid.
            Also stored in :attr:`Grid.results['velocity']` and
            :attr:`Grid.velocity` in `Grid.run()`.

        """
        self._requirements('mass')
        vel_name = self._all_calcs['velocity']['atomwise_name']['velocity']
        vel_comps = self.snap.vectors[vel_name]['comps']

        mass_name = self._all_calcs['mass']['atomwise_name']['mass']

        velocity = (self._sum_from_atom_value(
             (self.snap.get_vector(vel_comps, raw_data=(not self.use_sel)).T
             *self.snap.get_scalar(mass_name, raw_data=(not self.use_sel))).T
                                             ).T / self.mass.T).T
        util.clean(velocity, ~np.isfinite(velocity), 0.)
        return velocity

    def _calc_ke_count(self):
        """
        Calculates the gridded kinetic energy normalized by the atom count.

        Returns
        -------
        ke_count : :class:`numpy.ndarray`
            Gridded kinetic energy normalized by atom count.
            Also stored in :attr:`Grid.results['ke_count']` and
            :attr:`Grid.ke_count` in `Grid.run()`.

        """
        self._requirements('count')
        ke_name = self._all_calcs['ke_count']['atomwise_name']['ke']

        ke_count = self._sum_from_atom_value(
                     self.snap.get_scalar(ke_name, raw_data=(not self.use_sel))
                                            ) / self.count
        util.clean(ke_count, ~np.isfinite(ke_count), 0.)
        return ke_count

    def _calc_force_count(self):
        """
        Calculates gridded net force normalized by the atom count.

        Returns
        -------
        force_count : :class:`numpy.ndarray`
            Gridded net force normalized by atom count.
            Also stored in :attr:`Grid.results['force_count']` and
            :attr:`Grid.force_count` in `Grid.run()`.

        """
        self._requirements('count')
        force_name = self._all_calcs['force_count']['atomwise_name']['force']
        force_comps = self.snap.vectors[force_name]['comps']

        force_count = (self._sum_from_atom_value(
               self.snap.get_vector(force_comps, raw_data=(not self.use_sel))).T
               / self.count.T).T
        util.clean(force_count, ~np.isfinite(force_count), 0.)
        return force_count

    def _calc_temp(self):
        """
        Calculates gridded temperature, corrected for the local streaming
        velocity.

        Returns
        -------
        temp : :class:`numpy.ndarray`
            Temperature corrected by local streaming velocity.
            Also stored in :attr:`Grid.results['temp']` and
            :attr:`Grid.temp` in `Grid.run()`.

        """
        self._requirements('count', 'velocity')

        u = self.snap.meta['units']
        units_factor = units.MASS[u]*units.VELOCITY[u]**2 / units.ENERGY[u]

        # Getting effective streaming velocity for each atom.
        vel_grid_to_atom = self._reverse_mapping(self.velocity)

        mass_name = self._all_calcs['mass']['atomwise_name']['mass']
        vel_name = self._all_calcs['velocity']['atomwise_name']['velocity']
        vel_comps = self.snap.vectors[vel_name]['comps']

        vel_corr = self.snap.get_vector(vel_comps, raw_data=(not self.use_sel))\
                   - vel_grid_to_atom

        temp = units_factor * 2.0 * \
            self._sum_from_atom_value(
                        0.5*self.snap.get_scalar(mass_name, 
                                                 raw_data=(not self.use_sel))\
                                       *(vel_corr*vel_corr).sum(axis=1)
                                      ) / (3*self.count*self._KBOLTZ)

        util.clean(temp, ~np.isfinite(temp), 0.)
        return temp

    def _calc_temp_uncorr(self):
        """
        Calculates *uncorrected* gridded temperature.

        Returns
        -------
        temp_uncorr : :class:`numpy.ndarray`
            Gridded uncorrected temperature.
            Also stored in :attr:`Grid.results['temp_uncorr']` and
            :attr:`Grid.temp_uncorr` in :meth:`Grid.run()`.

        """
        self._requirements('count')

        u = self.snap.meta['units']
        units_factor = units.MASS[u]*units.VELOCITY[u]**2 / units.ENERGY[u]

        vel_name = self._all_calcs['velocity']['atomwise_name']['velocity']
        vel_comps = self.snap.vectors[vel_name]['comps']

        mass_name = self._all_calcs['mass']['atomwise_name']['mass']
        
        vel = self.snap.get_vector(vel_comps, raw_data=(not self.use_sel))

        temp_uncorr = units_factor * 2.0 * \
                      self._sum_from_atom_value( 0.5 *
                            self.snap.get_scalar(mass_name,
                                                 raw_data=(not self.use_sel))*\
                           (vel*vel).sum(axis=1)) / (3*self.count*self._KBOLTZ)

        util.clean(temp_uncorr, ~np.isfinite(temp_uncorr), 0.)
        return temp_uncorr

    def _calc_stress(self):
        """
        Calculates the gridded stress tensor.

        .. note::
           The gridded stress tensor is based on spatial averaging of existing
           atomic stress tensor.  Therefore, any streaming velocity correction
           must be applied in the original calculation of the atomic
           stress tensor.

        Returns
        -------
        stress : :class:`numpy.ndarray`
            Gridded stress tensor.
            Also stored in :attr:`Grid.results['stress']` and
            :attr:`Grid.stress` in `Grid.run()`.

        """
        stress_name = self._all_calcs['stress']['atomwise_name']['stress']
        stress_comps = self.snap.symm_tensors[stress_name]['comps']

        s = self.snap.get_symm_tensor(stress_comps, raw_data=(not self.use_sel))
        s = util.clean(s, ~np.isfinite(s), 0.)

        # atomwise stress from LAMMPS is multiplied by the volume
        stress = self._sum_from_atom_value(s) / (self.elem_volume)
        return stress


# -----------------------------------------------------------------------------
#       Grid Sub-classes
# -----------------------------------------------------------------------------


class RectGrid(Grid):
    """
    Regular grid over retangular region.

    Aligned to global coordinate system of the central simulation cell.

    Sub-class of :class:`Grid`.

    Parameters
    ----------

    snap : :class:`mdhandle.snap.Snap`
        Snap object.
    axes : iterable
        Total length of grid in each coordinate direction
        (e.g. [5., 5., 5.])
    grid_spacing : :class:`numpy.ndarray`
        Size of grid element in each coordinate direction.
    origin : :class:`numpy.ndarray`
        Coordinate of low corner of grid.
    offset : boolean
        If ``True``, the grid is offset so the first grid node
        is located at 0.5*:attr:`RectGrid.grid_spacing` from the
        :attr:`RectGrid.origin` rather than at the :attr:`RectGrid.origin`.
    mapping : string
        Name of atom-to-grid mapping function to match
        (e.g. ``'ngp'`` --? ``_mapping_ngp()`` ).
    calcs : iterable
        List of grid calculations to perform.  Names should match
        associated function ``_calc_NAME()``.
    wrap : boolean, [default=True]
        If ``True``, coordinates in active data set in snap are
        wrapped to the central simulation cell subject to PBC.
    selection : boolean, [default=True]
        If ``True``, selections are applied to snap object when getting
        atom coordinates.
    location : string, [default='Node']
        Grid location of calculated values.
        {'Node' } ``'Cell'`` location not currently implemented.

    Attributes
    -----------
    origin : :class:`numpy.ndarray`
        Coordinate of low corner of grid.
    grid_spacing : :class:`numpy.ndarray`
        Size of grid element in each coordinate direction.
    elem_volume : float
        Volume of grid element.
    num_elem : :clas:`numpy.ndarray`
        Number of grid elements in each coordinate direction.
    grid_type : string
        Name of type of grid (e.g. '3DRegular')
    snap : :class:`mdhandle.snap.Snap` object.
    use_selection :
        If ``True``, any selections in :attr:`Grid.snap` are used.
    atom_xyz :
        Atomic positions for atoms in :attr:`Grid.snap`.
    calcs :
        List of calculations to perform, which is processed by
        :meth:`Grid.setup_calculations` to generate :attr:`Grid.to_calc`.
    to_calc :
        List of calculations to do via :meth:`Grid.run_calc()`.
    scalars :
        List of scalars produced by grid calcs.
    vectors :
        List of vectors produced by grid calcs.
    tensors :
        List of tensors produced by grid calcs.
    symm_tensors :
        List of symmetric tensors produced by grid calcs.
    dtypes :
        dtype of data produced by grid calcs (i.e. 'int' or 'float')
    locations :
        Location of data produced by grid calcs (i.e. 'Node' or 'Cell').
    types :
        Type of output produced by grid calcs (i.e. 'scalars', 'vectors', ...)
    bins :
        Dictionary storing result of atoms-to-grid mapping.  Similar to
        :attr:`mdhandle.properties.cell_decomposition.SparseBinning.bins`.
        Must be supplied by the atom-to-grid mapping function.
    grid_to_atom : :class:`numpy.ndarray`
        Stores grid cells for each atom in the :attr:`Grid.snap`.
        Must be supplied by the atom-to-grid mapping function.
    _KBOLTZ :
        Boltzmann constant in units consistent with :attr:`Grid.snap`.
    _all_calcs :
        List of all implemented grid calculations.
    _matrix : :class:`numpy.ndarray` (3x3)
        The lenghts of the grid along each coordinate direction
        are the diagonal of :attr:`RectGrid._matrix`.
    _boxL : :class:`numpy.ndarray`, (3x1)
        Length of full grid in each coordinate direction.

    Notes
    ------

    * Each grid calculation (``_calc_NAME()``) stores its result in 
      :class:`Grid` attribute with name matching name in :attr:`Grid.to_calc`.
    * See :meth:`RectGrid._mapping_ngp` or :meth:`RectGrid._mapping_cic`
      for basic requirements of a mapping function.

    """
    grid_type = '3DRegular'

    def __init__(self, snap, axes, grid_spacing=np.ones(3, dtype=np.float),
                             origin=np.zeros(3, dtype=np.float),
                             offset=False,
                             mapping='ngp',
                             calcs=None,
                             wrap=True,
                             selection=True,
                             location='Node'):

        super(RectGrid,self).__init__(snap, calcs, wrap, selection, location)
        
        self.origin = origin.astype(np.float)
        self.grid_spacing = grid_spacing.astype(np.float)

        self._matrix = np.diag(axes).astype(np.float)

        # Testing that box defining RectGrid is orthogongal
        assert np.dot(self._matrix[:,0], self._matrix[:,1]) == 0, \
                    'RectGrid must be orthogonal.  Redefine _matrix.'
        assert np.dot(self._matrix[:,1], self._matrix[:,2]) == 0, \
                    'RectGrid must be orthogonal.  Redefine _matrix.'

        if offset is True:
            self.origin = origin + self.grid_spacing / 2.0

            new_matrix = np.zeros(self._matrix.shape)
            for i in range(self._matrix.shape[1]):
                col = self._matrix[:,i]
                new_len = np.linalg.norm(col) - self.grid_spacing
                new_matrix[:,i] = vectors.unit_vec(col)*new_len
            self._matrix = new_matrix

        self._boxL = np.array([ np.linalg.norm(i) for i in self._matrix.T])

        self.remainder = np.mod(self._boxL, self.grid_spacing)
        self._boxL -= self.remainder


        self.set_mapping_function(mapping)
        self.setup_calculations()
                      
    @property
    def num_elem(self):
        return (np.floor(1.0*self._boxL/self.grid_spacing) +
                                                      np.ones(3)).astype(np.int)

    @property
    def elem_volume(self):
        return np.prod(self.grid_spacing)

    def mapping_func(self):
        """
        Dummy function to override abstract function in :class:`Grid`
        base class so that class can be instantiated.
        
        The actual mapping function is assigned by 
        :meth:`RectGrid.set_mapping_function`.
        """
        pass

    def _mapping_ngp(self):
        """
        Peforms nearest-grid-point (NGP) zeroth order atom-to-grid mapping.

        The resulting mapping is stored in :attr:`RectGrid.bins`.

        The reverse grid-to-atom mapping is stored in
        :attr:`RectGrid.grid_to_atom`.

        References
        ----------
        * Hockney & Eastwood, Computer Simulation Using Particles, CRC Press,
          1988.

        """
        idx = ( (self.atom_xyz - self.origin) /
                                      self.grid_spacing).round().astype(np.int)

        self.grid_to_atom = np.empty((len(idx), 2), dtype=np.object)

        for (atom_id, idx_move) in enumerate(idx):
            idx_move = tuple(idx_move)

            # Avoiding atoms outside of the grid
            if np.any(idx_move >= self.num_elem) or\
               np.any(idx_move < np.zeros(3)):
                continue

            bin0 = self.bins.get(idx_move)
            if bin0 is None:
                bin0 = {'atoms':set(), 'weights':[]}
                # 1. hardcoded as this is the weighting
                bin0['atoms'].add(atom_id)
                bin0['weights'].append(1.0)
                self.bins[idx_move] = bin0
            else:
                # 1.0 hardcoded as this is the weighting
                bin0['atoms'].add(atom_id)
                bin0['weights'].append(1.0)

            # Storing reverse mapping.
            self.grid_to_atom[atom_id][0] = [idx_move]
            self.grid_to_atom[atom_id][1] = [1.0]

    def _mapping_cic(self):
        """
        Peforms cell-in-cloud (CIC) first order atom-to-grid mapping.

        The resulting mapping is stored in :attr:`RectGrid.bins`.

        The reverse grid-to-atom mapping is stored in
        :attr:`RectGrid.grid_to_atom`.

        References
        -----------
        * Hockney & Eastwood, Computer Simulation Using Particles, CRC Press,
          1988.

        """
        idx = ( self.atom_xyz-self.origin ) / self.grid_spacing
        idx_down = np.floor(idx).astype(np.int)

        self.grid_to_atom = np.empty((len(idx), 2), dtype=np.object)
        for movement in ( (0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0),
                          (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1) ):

            # Fractional distance to grid point defined by movement
            # which iterates through the NGP.
            diff_frac = 1 - np.abs(idx - (idx_down + movement))

            for (atom_id, idx_move) in enumerate(idx_down + movement):
                idx_move = tuple(idx_move)

                if np.any( idx_move >= self.num_elem) or\
                   np.any( idx_move < np.zeros(3)):
                    continue

                bin0 = self.bins.get(idx_move)
                if bin0 is None:
                    bin0 = {'atoms':set(), 'weights':[]}
                    # 1. hardcoded as this is the weighting
                    bin0['atoms'].add(atom_id)
                    bin0['weights'].append( np.prod( diff_frac[atom_id]) )
                    self.bins[idx_move] = bin0
                else:
                    bin0['atoms'].add(atom_id)
                    bin0['weights'].append( np.prod( diff_frac[atom_id]) )

                # Reverse mapping.
                if self.grid_to_atom[atom_id, 0] is None:
                    # First time through create lists at each element.
                    self.grid_to_atom[atom_id,0] = [idx_move]
                    self.grid_to_atom[atom_id,1] = [np.prod(diff_frac[atom_id])]
                else:
                    # Otherwise append to existing lists.
                    self.grid_to_atom[atom_id,0].append(idx_move)
                    self.grid_to_atom[atom_id,1].append(
                                                   np.prod(diff_frac[atom_id]))


# -----------------------------------------------------------------------------
#                   Setting module level constants
# -----------------------------------------------------------------------------

IMPLEMENTED_GRIDS = get_implemented_grids()

# =============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
