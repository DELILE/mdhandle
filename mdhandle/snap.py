# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Continer for holding information on a single simulation snapshot.

"""


# TODO: Abstract out direct use of XDMF writer.
# TODO: Consider adding scalars list like vectors to metadata
# TODO: Add dedicated getter for binned data arrays.  Currently using
#       Snap.get_scalar().
# TODO: Use reader/writer rather than operating on the HDF file itself.
#       or make Snap the only place that files are opened and managed.
# TODO: Change get_scalar, get_vector, etc. so can be used by gridded datasets
#       - Currently data in gridded datasets stored as a single bloc
#         so best to use get_scalar.
#       - Need to use col_mapping meta to understand what to do.


#---------------------------------------------------------------------------

import os

import tables as pyT
import numpy as np

from mdhandle import units
from mdhandle.writers.xdmf import XDMFWriter
from mdhandle.logger import Logger
from mdhandle.sim_cell import SimCell
from mdhandle import vectors
from mdhandle import tensors
from mdhandle import settings

#---------------------------------------------------------------------------


class Snap(object):
    """
    A class to hold a single snapshot (i.e. heavy atomwise data and light
    metadata).

    Parameters
    -----------
    fn : string
        Filename of HDF file containing simulation snapshot data.
    dataset_name : string, [default='rawSimResults']
        Name of data set within HDF file.  Datasets are stored
        at the root of the HDF file.  Default dataest ``'rawSimResults'``
        should always be present and contains raw atomwise output.
    is_new_data : boolean, [default=False]
        Boolean flag. True if snap is house new data rather than
        reading existing dataset from disk.
    verbose : boolean, [default=False]
        Boolean flag for verbosity of user input/output.

    Attributes
    ----------
    filename : string
        Absolute filename of snapshot.
    fn_base : string 
        ``filename`` without absolute path.
    directory : string
        Absoulte path of containing directory.
    active_dataset_name : string
        Name of active dataset within snapshot file.
    active_dataset : :class:`pytables.tables.group.Group`.
        Active dataset object.
    selections : list of :class:`mdhandle.selection.Selection`
        List of selection objects.
    meta : dict
        Dictionary with all metadata in active dataset.
    sim_cell :  :class:`mdhandle.sim_cell.SimCell` object.
        Contains geometric information about the simulation domain.
    col_mapping: dict
        Dictionary of all atom/grid scalars and metadata.
    vectors : dict
        Dictionary of vectors and metadata.
    tensors : dict
        Dictionary of tensors and metadata.
    symm_tensors : dict
        Dictionary of symmetric tensors and metadata.
    verbose : boolean
        Flag for user I/O verbosity.
    _mask : boolean :class:`numpy.ndarray` 
        ``True`` entries mark atoms included in active selection.
    writer : :class:`mdhandle.writers.xdmf.XDMFWriter`
        XDMF writer object.

    Methods
    --------
    add_selection(new_selection)
        Adds new selection object.
    reset_selections()
        Clears the list of selections objects.
    get_datasets()
        Returns list of datasets within HDF file.
    set_active_dataset(name, create=False)
        Sets active dataset from within HDF file.
    gather_data()
        Retrieves metadata from active dataset.
    dict_metadata()
        Returns dictionary of all metadata in active dataset.
    get_metadata(name)
        Returns single metadata item.
    get_scalar(name, raw_data=False)
        Returns scalar collection (i.e. scalar for each atom or grid point).
    get_vector(comp_names, raw_data=False)
        Returns vector collection.
    get_tensor(comp_names, raw_data=False)
        Returns tensor collection.
    get_symm_tensor(comp_names, raw_data=False)
        Returns symmetric tensor collection.
    _use_selections()
        Applies selections to atomwise (or grided) data.

    See Also
    --------
    :class:`mdhandle.simcontainer.SimContainer`

    """

    def __init__(self, fn, dataset_name='rawSimResults',
                                            is_new_data=False, verbose=False):
        self.verbose = verbose
        self.logger = Logger(verbose)

        self.filename = fn

        self.fn_base = os.path.basename(fn)
        self.directory = os.path.dirname(fn)
        
        self.meta = {}

        # List of Selection object used to extract only part of snap.
        self.selections = []
        self._mask = None

        # Simulation output metadata information.
        self.col_mapping = {}
        self.vectors = {}
        self.tensors = {}
        self.symm_tensors = {}

        # Geometry meta-data.
        self.sim_cell = None
        self.active_dataset = None
        self.active_dataset_name = ''


        self.writer = XDMFWriter(self)
        self.f_hdf = self.writer.open_file()

        # Must be done after setting self.writer
        self.set_active_dataset(dataset_name, create_new=is_new_data)

    def __repr__(self):
        """
        Short and identifiable string represenation of :class:`Snap` object.

        Returns
        -------
        string_rep : string
            String representation of object (< 80 chars).

        """
        string_rep = 'Snap %s - %s - id: %s' % (self.fn_base,
                                            self.active_dataset_name, id(self))
        return string_rep

    def __str__(self):
        """
        Prints multi-line textual representation of :class:`Snap` object.

        Returns
        -------
        string_rep : string
            Multiple line string represenation of object.

        """
        string_rep = """
===  Snap Object ====
Sim Name: \t %s
Dataset Name: \t %s
Units: \t %s
Num. Atoms: \t %s
""" % (self.meta['sim_name'], self.active_dataset_name, 
       self.meta['units'], self.meta['num_atoms'])
        return string_rep

    def set_active_dataset(self, name, reset_meta=False, create_new=False):
        """
        Sets the active dataset within HDF file.
        The active dataset is the one available via other Snap methods.
        Dataset is selected by name at root ('/') of HDF file.

        Parameters
        ----------
        name : string
            Name of active dataset.
        create : boolean
            Boolean flag used for creating a new dataset.
        reset_meta : boolean
            If ``True``, metadata and selections in snap are reset.

        """
        try:
            self.active_dataset = self.f_hdf.getNode('/', name)
            self.active_dataset_name = name
            
            # Resetting metadata
            if reset_meta is True:
                if self.selections != []:
                    self.logger.user_message('Resetting metadata and\
                                                                selections.')
                    self.reset_selections()

                self.gather_data()

        except pyT.NoSuchNodeError:
            if create_new is True:
                self.writer.hdf_create_dataset(name)
                self.active_dataset = self.f_hdf.getNode('/', name)
                self.active_dataset_name = name
                self.logger.user_message('Created new dataset %s' % name)
                
                # Reinitializing attributes for newly created datasets
                self.reset_selections()
                self.meta = {}
                self.col_mapping = {}
                self.vectors = {}
                self.tensors = {}
                self.symm_tensors = {}
            else:
                self.logger.error('Requested dataset %s does not exist in %s' % 
                                        (name, self.filename))

    def get_datasets(self):
        """
        Walks the root of related HDF file and return list of datasets.

        Returns
        -------
        datasets : list
            List of datasets contained within :class:`Snap` file.

        """
        datasets = [i._v_name for i in self.f_hdf.walkGroups() 
                                                        if i._v_name != '/']
        return datasets

    def gather_data(self):
        """
        Gather and calculate metadata from active dataset.
        
        Also used to ensure that the HDF file associated with the snapshot
        is active.

        """
        if self.writer.is_open() is False:
            self.f_hdf = self.writer.open_file()
            self.set_active_dataset(self.active_dataset_name)

        self.meta = self.dict_metadata()
        self.col_mapping = self.meta['col_mapping']
        self.vectors = self.meta['vectors']
        self.tensors = self.meta['tensors']
        self.symm_tensors = self.meta['symm_tensors']

        # Creating SimCell object to store information re simulation cell
        box_len = np.array((self.meta['xhi'] - self.meta['xlo'],
                         self.meta['yhi'] - self.meta['ylo'],
                         self.meta['zhi'] - self.meta['zlo'] ))
        self.sim_cell = SimCell(box_len, pbc=self.meta['pbc'])

        if self.meta['units'] not in units.UNITS_ALLOWED:
            self.logger.error('Snap units %s not allowed' % self.meta['units'])


    def get_metadata(self, name):
        """
        Returns value for matadata item associated with  ``name``.
        Data is retreived from the current active dataset.

        Parameters
        ----------
        name : string
            Name of metadata item to retrieve.

        Returns
        -------
            Metadata item associated with ``name``.

        """
        return getattr(self.active_dataset._v_attrs, str(name))

    def dict_metadata(self):
        """
        Returns dictionary indexed by property name, for all matadata
        stored in HDF file for snapshot.

        Returns
        --------
        meta_dict : dict
            Dictionary of ``string``-`value`` pairs for all
            metadata associated with the active snapshot.

        """
        meta_dict = {}
        for item in self.active_dataset._v_attrs._f_list():
            meta_dict[item] = self.get_metadata(item)
        return meta_dict

    def add_selection(self, new_selection):
        """
        Adds a :class:`mdhandle.properties.selection.Selection` object to
        snapshot.

        Selections are stored in a list at :attr:`Snap.selections` and act to
        sequentially select sub-sets of atoms within the snapshot.

        Parameters
        ----------
        new_selection : :class:`mdhandle.properties.selection.Selection`
            Selection object object used to choose a sub-set of
            atoms within the snapshot.

        """
        self.selections.append(new_selection)

    def reset_selections(self):
        """
        Removes all existing selections stored in :attr:`Snap.selections`.

        See Also
        ----------
        :meth:`Snap.add_selection`.

        """
        self.selections = []

    def _use_selections(self):
        """
        Applies list of selections stored in self.selections.

        Selections are applied sequentially from ``Snap.selections[0]``,
        ``Snap.selections[1]``, ... to filter atoms within the snapshot.

        The result of all selection is stored in :attr:`Snap._mask`, which is an
        array with a value of either ``0`` (unselected) or ``1`` (selected),
        which is used to filter any atomwise results (i.e. scalar, vector,
        tensor collections)

        """
        # Setting default mask: All atoms are selected by np.ones.
        if self.meta['grid_type'] == 'atoms':
            self._mask = np.ones(self.meta['num_atoms'],  dtype=bool)
        elif self.meta['grid_type'] == '3DRegular':
            # TODO: change the way to get the shape of the grid arrays
            self._mask = np.ones(self.get_scalar('count').shape, dtype=bool)

        for sel in self.selections:
            # The snap object (i.e. self is passed through to the masking)
            self._mask *= sel.calculate_mask(self)

    def good_grid(self):
        """
        If the active dataset in HDF file is a known type as stored
        in :attr:`mdhandle.settings.IMPLEMENTED_GRIDS`, returns ``True``.

        """
        if self.meta['grid_type'] in settings.IMPLEMENTED_GRIDS:
            return True
        else:
            self.logger.warning('Grid type %s unknown' % self.meta['grid_type'])
            return False


    def get_scalar(self, name, raw_data=False):
        """
        Gets and returns atwomwise scalar (:class:`numpy.ndarray`, (Nx1)).

        If ``raw_data`` is ``False`` and :attr:`~Snap.selections` contains
        active selections, the data associated with the selected atoms
        is returned.

        Parameters
        ----------
        name : string
            Name of scalar in HDF file.
        raw_data : boolean, [default=False]
            If ``True``, avoids the use of any existing selections.

        Returns
        --------
        scalar : :class:`numpy.ndarray`, (Nx1)

        """
        arr = self.active_dataset._f_getChild(name)

        # No selections are applied.
        if (self.selections == [] or raw_data ) is True:
            scalar = arr.read()
        # Applying selections and extracting data for selected atoms.
        else:
            self._use_selections()
            scalar = np.extract(self._mask, arr.read())
        return scalar

    def get_vector(self, comp_names, raw_data=False):
        """
        Gets and returns atwomwise vector collection
        (:class:`numpy.ndarray`, (Nx3)).

        If ``raw_data`` is ``False`` and :attr:`~Snap.selections` contains
        active selections, the data associated with the selected atoms
        is returned.

        Parameters
        ----------
        comp_names : iterable
            List of names of scalar vector components in HDF file.
        raw_data : boolean, [default=False]
            If ``True``, avoids the use of any existing selections.

        Returns
        --------
        vec : :class:`numpy.ndarray`, (Nx3)

        """
        assert len(comp_names)==3

        return vectors.build_from_scalars((
                    self.get_scalar(comp_names[0], raw_data),
                    self.get_scalar(comp_names[1], raw_data),
                    self.get_scalar(comp_names[2], raw_data)  ))

    def get_tensor(self, comp_names, raw_data=False):
        """
        Gets and returns atwomwise tensor collection
        (:class:`numpy.ndarray`, (Nx9)).

        If ``raw_data is False`` and :attr:`~Snap.selections` contains
        active selections, the data associated with the selected atoms
        is returned.

        Parameters
        ----------
        comp_names : iterable
            List of names of scalar tensor components in HDF file.
        raw_data : boolean, [default=False]
            Boolean flag, which if `True` avoids the use of any existing
            selections.

        Returns
        --------
        tens : :class:`numpy.ndarray`, (Nx9)

        """
        assert len(comp_names) == 9

        tens = tensors.build_from_scalars((
                    self.get_scalar(comp_names[0], raw_data),
                    self.get_scalar(comp_names[1], raw_data),
                    self.get_scalar(comp_names[2], raw_data),
                    self.get_scalar(comp_names[3], raw_data),
                    self.get_scalar(comp_names[4], raw_data),
                    self.get_scalar(comp_names[5], raw_data),
                    self.get_scalar(comp_names[6], raw_data),
                    self.get_scalar(comp_names[7], raw_data),
                    self.get_scalar(comp_names[8], raw_data)  ))
        return tens

    def get_symm_tensor(self, comp_names, raw_data=False):
        """
        Gets and returns atwomwise symmetric tensor collection
        (:class:`numpy.ndarray`, (Nx6)).

        If ``raw_data is False`` and :attr:`~Snap.selections` contains
        active selections, the data associated with the selected atoms
        is returned.

        Parameters
        ----------
        comp_names : iterable
            List of names of scalar tensor components in HDF file.
        raw_data : boolean, [default=False]
            Boolean flag, which if ``True`` avoids the use of any existing
            selections.

        Returns
        --------
        symm_tens : :class:`numpy.ndarray`, (Nx6)

        """
        assert len(comp_names) == 6

        symm_tens = tensors.build_from_scalars((
                    self.get_scalar(comp_names[0], raw_data),
                    self.get_scalar(comp_names[1], raw_data),
                    self.get_scalar(comp_names[2], raw_data),
                    self.get_scalar(comp_names[3], raw_data),
                    self.get_scalar(comp_names[4], raw_data),
                    self.get_scalar(comp_names[5], raw_data)  ))
        return symm_tens

# =============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    main()
