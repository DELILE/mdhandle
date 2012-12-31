# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Wrapper on set of simulation results.  :class:SimContainer objects are used
to hold and analyze a number of simulation snapshots.

"""

# TODO: Add metadata organizing functions
# TODO: Add method to delete and/or pop snaps from SimContainer collection
# TODO: Add check for SimCell object constant through snaps in SimContainer.
# TODO: Add ability to act only slice of slist within generators


#-----------------------------------------------------------------------------

import sys
import os
import resource

from mdhandle.snap import Snap
from mdhandle.readers.log import Log
from mdhandle.logger import Logger
from mdhandle import utilities as util

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

class SimContainer(object):
    """
    Wrapper on set of simulation results.  :class:`SimContainer` objects are 
    used to hold and analyze a number of simulation snapshots.

    Parameters
    -----------
    flist : iterable
        List of HDF files containing a snapshots of simulation restuls.
        Each file in ``flist`` becomes a :class:`Snap` object.
    dataset_name : string, [default='rawSimResulsts']
        Name of dataset to activate within each snapshot.
        Default ``'rawSimResulsts'`` is always present and contains
        the atomwise MD simulation output.
    verbose : boolean, [default=False]
        Boolean flag for the verbosity of user input/output.

    Attributes
    -----------
    slist : list
        List of :class:`mdhandle.snap.Snap` objects in :class:`SimContainer`.
    _max_open : int
        Maximum number of files that can be open in Python.
    nsnaps : int
        Number of snapshots within the :class:`SimContainer`.
    times : list
        List of time stamps for snaps within :class:`SimContainer`.
    selections : list
        List of selections active in snaps within :class:`SimContainer`.
    active_dataset_name : string
        Name of active dataset within snaps within :class:`SimContainer`.
    meta : dict
        Mapping of all simulation metadata.
    col_mapping : dict
        Metadata mapping for scalar collections.
    vectors : dict
        Dictionary of vector collections.
    tensors : dict
        Dictionary of tensor collections.
    symm_tensors : dict
        Dictionary of symmetric tensor collections.
    is_constant : dict
        Dictionary sorring boolean values for constant metadata in container.
    sim_cell : :class:`mdhandle.sim_cell.SimCell`
    LAMMPS_log : :class:`mdhandle.readers.log.Log`
        LAMMPS log file handle.
    logger : :class:`mdhandle.logger.Logger`
    verbose : boolean
        If ``True``, verbose user I/O.

    """

    def __init__(self, flist, dataset_name='rawSimResults', verbose=False):
        self.verbose = verbose
        self.logger = Logger(self.verbose)

        if self.verbose is True:
            if len(flist) == 0:
                self.logger.user_message('SimContainer: List of files is empty')
            else:
                self.logger.user_message('List of files in simContainer:')
                util.print_file_list(flist)

        self._max_open = resource.getrlimit(resource.RLIMIT_NOFILE)[0]

        # Adding Snap objects to SimContainer
        self.slist = []
        self.nsnaps = 0
        self.times = []
        self.selections = []

        # metadata  - Should be constant within SimContainer
        self.meta = {}
        self.col_mapping = {}
        self.vectors = {}
        self.tensors = {}
        self.symm_tensors = {}

        self.is_constant = {}
        self.sim_cell = None

        self.LAMMPS_log = None

        self.active_dataset_name = ''
        if isinstance(dataset_name, basestring):
            self.set_active_dataset(dataset_name)
        else:
            self.logger.error('Dataset name must be a string')

        flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists,
                                                                   self.verbose)
        for fn in flist:
            self.add_snap(fn)

        # In case any snaps didn't work.
        self.flist = [snap.filename for snap in self.slist]

        self.gather_data()
        # Testing for constant metadata across SimContainer.
        self.constant_props_test()

    def __repr__(self):
        """
        Prints single line unique text representation of
        :class:`SimContainer` object.
        
        Returns
        -------
        string_rep : string
            Single line string representation of :class:`simContainer`.

        """
        string_rep = '%s - %s - id: %s' % \
                    (self.meta['sim_name'], self.active_dataset_name, id(self))
        return string_rep

    def __str__(self):
        """
        Prints multiple line textual representation of 
        :class:`~mdhandle.simcontainer.SimContainer` object.

        Returns
        -------
        string_rep : string
            Multiple line string representation of :class:`SimContainer`.

        """
        string_rep = """
======  SimContainer Object ======
Sim. Name: \t\t %s
Dataset Name: \t\t %s
Snaps: \t\t %s
Units: \t\t %s
Units Constant? \t\t %s
Num. Atoms: \t\t %s
Num. Atoms. Constant? \t\t %s
Box Constant? \t\t %s
""" % (self.meta['sim_name'], self.active_dataset_name, str(self.flist)[:60],
       self.meta['units'], self.is_constant['units'], 
       self.meta['num_atoms'], self.is_constant['num_atoms'],
       self.is_constant['box'])

        return string_rep

    def set_active_dataset(self, name, reset_meta=False):
        """
        Sets active dataset in :class:`SimContainer` and attached snapshots.

        The active dataset is the one available via other methods.
        Dataset is selected by name at root ('/') of HDF file.

        Parameters
        ----------
        name : string
            Name of active dataset.
        reset_meta : boolean
            If ``True``, metadata and selections in snap are reset.

        """
        self.active_dataset_name = name

        for snap in self.slist:
            snap.gather_data()
            snap.set_active_dataset(name, reset_meta=reset_meta)

        if reset_meta is True:
            if self.selections != []:
                self.logger.user_message('Resetting metadata and selections.')
                self.reset_selections()
            if len(self.slist) != 0:
                self.gather_data()

    def add_snap(self, new_name, test_constant=True):
        """
        Adds a new snap to SimContainer object.

        Parameters
        -----------
        new_name : string
           File name of :class:`Snap` file (HDF).
        test_constant : boolean
           If True, new snap is tested against existing snaps to see if
           has consistent properties (simulation domain dimensions, number
           of atoms, boundary conditions, unit type).

        """
        # Active data set is not set - cannot proceed with adding snap.
        if (self.active_dataset_name is None
                                        or self.active_dataset_name == ''):
            self.logger.warning('Cannot add snapshot without dataset.')
            return

        try:
            snap_handle = Snap(new_name, self.active_dataset_name)
        # IOError thrown py pyTables on bad read.
        except IOError:
            self.logger.warning('Skipping %s: Not HDF file, or does not exist.'
                                                                % new_name)
        # Adding any existing selections to new snap
        for sel in self.selections:
            snap_handle.add_selection(sel)

        # Close snap file to avoid having too many open at one time.
        snap_handle.writer.cleanup()
        self.slist.append(snap_handle)
        # Incrementing the number of snaps.
        self.nsnaps += 1

        # Testing if new snap changes constant meta data tests
        if test_constant is True:
            self.constant_props_test()

    def gather_data(self, test_constant=True):
        """
        Gathers metadata based on data from first snapshot
        (i.e. at ``SimContainer.slist[0]``)

        Parameters
        -----------
        test_constant : boolean
            If True, new snap is tested against existing snaps to see if
            has consistent properties (simulation domain dimensions, number
            of atoms, boundary conditions, unit type)

        """
        # Perform gather_data(...) on each of attached snaps.
        self.slice_action_do('gather_data')
        
        self.meta = self.slist[0].meta
        self.sim_cell = self.slist[0].sim_cell
        self.col_mapping = self.slist[0].col_mapping
        self.vectors = self.slist[0].vectors
        self.tensors = self.slist[0].tensors
        self.symm_tensors = self.slist[0].symm_tensors

        self.times = [snap.meta['time'] for snap in self.slist]

        # Testing if metadata is constant across attached snaps
        if test_constant is True:
            self.constant_props_test()

    # ----------------- GETTING DATA FROM SNAPS ------------------------------

    def slice_action_do(self, action, *args, **kwargs):
        """
        Calls function named ``action`` on all attached snaps.
        
        The snap method ``action`` is called immediately.  This
        function does not produce a generator.

        .. note::
           This function can only act on the snaps in :attr:`slist`.
           Nothing is returned.
        
        Parameters
        ----------
        action : string
            Name of :class:`~mdhandle.snap.Snap` function.
        *args, *kwargs :
            Arguments and keyword arguments to ``action`` function.

        """
        for snap in self.slist:
            snap.gather_data()

            # Special case - gather data is already called.
            if action == 'gather_data':
                pass
            else:
                func = getattr(snap, action, None)
                if func is not None:
                    func(*args, **kwargs)
                else:
                    self.logger.warning('Function %s not available' % action)
            snap.writer.cleanup()

    def slice_action_return(self, action, *args, **kwargs):
        """
        Returns generator containing snap method ``action``.
        
        To get the result of the snap method the generator must be used.
        
        This method is required in cases where the return value is required,
        especially in cases where the returned value is large
        (i.e. atomwise or gridded datad)

        """
        for snap in self.slist:
            snap.gather_data()

            func = getattr(snap, action, None)
            if func is None:
                self.logger.warning('Snap function %s not available' % action)

            yield func(*args, **kwargs)

            # Closing if necessary to prevent accumulation of open files.
            if snap.writer.is_open() is True:
                snap.writer.cleanup()

    def slice_generator(self):
        """
        Returns generator used for iterating over 
        :class:`~mdhandle.snap.Snap` functions.
        
        .. note::
           Users must be careful to close the snap files as they are used
           in cases where the number of snaps within :attr:`slist` is
           larger than the maximum number of open files in Python.
        
        """
        for snap in self.slist:
            snap.gather_data()
            yield snap

    def slice_metadata(self, prop_name):
        """
        Get item of metadata from across entire :class:`SimContainer`.

        Parameters
        -----------
        prop_name : string
            Name of :class:`~mdhandle.snap.Snap` metadata property to get.

        Returns
        --------
            Ordered list of metadata value for each snapshot.

        """ 
        return [snap.meta[prop_name] for snap in self.slist]

    def slice_heavy_scalar(self, scalar_name, raw_data=False):
        """
        Get scalar collection (Nx1) from each attached 
        :class:`~mdhandle.snap.Snap`.

        Parameters
        ----------
        scalar_name : string
            Name of scalar attribute.
        raw_data : boolean, [default=False]
            If ``True``, avoids the use of any existing selections.

        Returns
        --------
        scalar_gen : generator
            Generator for iterating through snap files.  Each call to the
            generator returns a scalar collection.

        """
        scalar_gen = self.slice_action_return('get_scalar', 
                                              scalar_name, raw_data)
        return scalar_gen

    def slice_heavy_vector(self, comp_names, raw_data=False):
        """
        Get vector collection (Nx3 numpy.ndarray) from each attached 
        :class:`~mdhandle.snap.Snap`.

        Parameters
        -----------
        comp_names : iterate
            Names of vector components for vector collection.
        raw_data : boolean, [default=False]
            If ``True``, avoids the use of any existing selections.

        Returns
        --------
        vec_gen : generator
            Generator for iterating through snap files.  Each call to the
            generator returns a vector collection.

        """
        vec_gen = self.slice_action_return('get_vector', comp_names, raw_data)
        return vec_gen

    def slice_heavy_tensor(self, comp_names, raw_data=False):
        """
        Get tensor collection (Nx9 :class:`numpy.ndarray`) from each 
        attached :class:`~mdhandle.snap.Snap`.

        Parameters
        -----------
        comp_names : iterable
            Names of tensor components for tensor collection.
        raw_data : boolean, [default=False]
            If ``True``, avoids the use of any existing selections.

        Returns
        --------
        tens_gen : generator
            Generator for iterating through snap files.  Each call to the
            generator returns a tensor collection.

        """
        tens_gen = self.slice_action_return('get_tensor', comp_names, raw_data)
        return tens_gen

    def slice_heavy_symm_tensor(self, comp_names, raw_data=False):
        """
        Get symmetic tensor collection (:class:`numpy.ndarray`, (Nx6)) 
        from each attached :class:`~mdhandle.snap.Snap`.

        Parameters
        ----------
        comp_names : iterable
            Names of symmetric tensor components for symmetric
            tensor collection.
        raw_data : boolean, [default=False]
            If ``True``, avoids the use of any existing selections.

        Returns
        --------
        symm_gen : generator
            Generator for iterating through snap files.  Each call to the
            generator returns a symmetric tensor collection.

        """
        symm_gen = self.slice_action_return('get_symm_tensor', 
                                            comp_names, raw_data)
        return symm_gen

    # -------------------------- SELECTIONS ----------------------------------

    def add_selection(self, new_selection):
        """
        Adds a :class:`~mdhandle.properties.selection.Selection` object to
        :class:`SimContainer` and contained :class:`~mdhandle.snap.Snap` 
        objects.

        Selections are stored in a list at :attr:`SimContainer.selections` and
        act to sequentially select sub-sets of atoms.

        Parameters
        ----------
        new_selection : :class:`mdhandle.properties.selection.Selection`
            Selection object object used to choose a sub-set of
            atoms within the snapshot.

        """
        self.selections.append(new_selection)
        self.slice_action_do('add_selection', new_selection)

    def reset_selections(self):
        """
        Removes all existing selections stored in self.selections from
        simContainer and contained :class:`mdhandle.snap.Snap` objects.

        See Also
        ----------
        :meth:`SimContainer.add_selection` for more information.

        """
        self.selections = []
        self.slice_action_do('reset_selections')

    # ----------------- TESTING SNAP DATA ------------------------------------

    def constant_props_test(self):
        """
        Runs collection of metadata tests.

        The result of each test is stored in related 
        :class:`SimContainer` attribute.

        **:class:`SimContainer`  items tested**:
        
        * Simulation cell dimensions.
        * Number of atoms.
        * Unit system.
        * Attributes (scalars, vectors, tensors, symmetric tensors).

        """
        self.box_constant_test()
        self.meta_constant_test('col_mapping')
        self.meta_constant_test('num_atoms')
        self.meta_constant_test('units')

    def box_constant_test(self):
        """
        Tests if simulation box dimensions are constant for snaps within
        the :class:`SimContainer` object.

        **Assumptions**:
        
        * Cuboid simulation domain.

        """
        num_fail = 0

        self.logger.user_message('Running: box_constant_test')

        for coord in ['x', 'y', 'z']:
            for direc in ['lo', 'hi']:
                test_value = True
                test_value = self.meta_constant_test(coord+direc)
                if test_value is False:
                    num_fail += 1
        if num_fail > 0:
            self.logger.user_message('Result: False')
            self.is_constant['box'] = False
        else:
            self.logger.user_message('Result: True')
            self.is_constant['box'] = True

        return self.is_constant['box']

    def meta_constant_test(self, prop_name):
        """
        Test if a piece of metadata is constant across set of snapshots.

        Parameters
        ----------
        prop_name : string
            Name of metadata property for testing.

        Returns
        --------
        test_value : boolean
            Boolean flag for the value of the test.

        """
        self.logger.user_message('Running: %s constant test' % prop_name)

        test_value = True
        meta_list = self.slice_metadata(prop_name)

        if len(meta_list) == 0:
            self.logger.warning('SimContainer does not contain any data.')
            test_value = True

        elif len(meta_list) == 1:
            self.logger.warning('SimContainer contains only single\
                                            snapshot.\n\tMetadata is constant')
            test_value = True

        # More than one snap - iterate over gathered metadata
        elif len(meta_list) > 1:
            first_item = meta_list[0]
            for item in meta_list:
                if item != first_item:
                    test_value = False
                    break
        # Only if len(meta_list) < 0
        else:
            self.logger.error('simContainder metadata slice has \
                                                            negative length.')

        # Setting attribute and returning value
        self.logger.user_message('Result: %s' % test_value)
        self.is_constant[prop_name] = test_value
        return test_value

    # ----------------- MISC FUNCTIONS ----------------------------------------

    def get_log(self, log_name):
        """
        Attach LAMMPS log file to :class:`SimContainer`.

        Does not contain any logic to ensure that the log file is logically
        related to the snaps within the SimContainer.

        Uses :class:`mdhandle.readers.log.Log`.

        Parameters
        ----------
        log_name : string
            Absolute path to log file.

        """
        self.LAMMPS_log = Log(log_name)


# ============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
