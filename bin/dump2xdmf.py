    #!/usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Convert LAMMPS dump files into set of XDMF files (HDF and .xml).  A pair of
.xmf and HDF file is produced for each LAMMPS dump file and timestep.

Input LAMMPS dump files should have original file header intact.

Returns ``0`` if successful, ``-1`` otherwise.

**Credits**:

* dump file in Pizza package, (http://www.sandia.gov/~sjplimp/pizza.html).

.. note::
   Conversion is non-destructive. Original LAMMPS dump files are left
  unchanged and can be deleted once XDMF files are verified to be ok.

"""

# TODO: Extend capabilities beyond simple atomistic datasets (bonds, rotations)
# TODO: Add flexibility to handle more than one snapshot per dump file.
#          - Require that the file is scanned for the header
# TODO: Add support for compressed files (.gz, .tar.gz)

import os
import sys
from optparse import OptionParser
from multiprocessing import Process

import numpy as np

from mdhandle.readers import column_mapping as col
from mdhandle.readers.lammps_dump import LAMMPS_dump
from mdhandle import utilities as util
from mdhandle import settings
from mdhandle.units import UNITS_ALLOWED
from mdhandle.snap import Snap
from mdhandle.logger import Logger

# ----------------------------------------------------------------------------

logger = Logger()

# ----------------------------------------------------------------------------


def process_single(fn, col_mapping, vectors, tensors,
                       symm_tensors, sim_name, units, timestep, mass, pbc):
    """
    Processes individual LAMMPS dump file into XDMF format
    (i.e. HDF and .xmf files).

    Data is put into default group ``'/rawSimResults'`` within HDF file.

    Parameters
    -----------
    fn : string
        Filename of LAMMPS dump file.
    col_mapping : dict
        Dictionary with mapping of column name, number, data type and location.
    vectors : dict
        Dictionary with makeup of vectors.
    tensors : dict
        Dictionary with makeup of tensors.
    symm_tensors : dict
        Dictionary with makeup of symmetric tensors
    sim_name : string
        Name of simulation.
    units : string, {'lj' | 'metal' | 'si'}
        Units system identifier.
    timestep : float
        Integration timestep used in simulation.
    mass : iterable
        List of atomic masses found in simulation in AMU.
        Used to create an atomwise mass vector to be stored in HDF file.
    pbc : list or tuple
        Flag for periodic boundary conditions along each coordinate direction.
        (e.g. ``(1, 0, 1) == (True, False, True)`` is periodic along
        the x and z direction)

    """

    reader = LAMMPS_dump(fn, col_mapping)

    if reader.meta['ncol'] != len(col_mapping):
        logger.error(
        """
Column mapping does not agree with columns in LAMMPS dump.
LAMMPS dump --> %s
Column Mapping --> %s
        """ % (reader.meta['ncol'], len(col_mapping)))

    fn_hdf = os.path.splitext(fn)[0] + '.h5'
    snap = Snap(fn_hdf, dataset_name='rawSimResults', is_new_data=True)

    # Writing atomwise data:
    for name in col_mapping:
        snap.writer.hdf_create_array(name,
                             reader.atom_data[:, (col_mapping[name]['col'])],
                             dtype=col_mapping[name]['dtype'],
                             location=col_mapping[name]['location'],
                             update_meta=False)

    col_mapping_no_num = col.convert_from_num(col_mapping)

    # Writing atomwise mass vector if doesn't exist.
    # Same condition as 'mass' not in col_mapping
    if mass is not None:
        type_array = reader.atom_data[:, col_mapping['type']['col']]
        unique_types = np.unique(type_array)
        mass_array = np.zeros(type_array.size)

        if len(mass) != len(unique_types):
            logger.error('Number of masses and atom types incongruent')

        for (indiv_type, m) in zip(unique_types, mass):
            mass_array += (type_array == indiv_type)*m

        col_mapping_no_num['mass'] = {'dtype': 'float', 'location': 'Node'}
        snap.writer.hdf_create_array('mass', mass_array,
                                            dtype='float',
                                            location='Node',
                                            update_meta=False)
    else:
        logger.user_message('Mass is already in atomwise dump results')

    # Writing metadata
    snap.writer.hdf_add_meta('time', reader.meta['time'])
    snap.writer.hdf_add_meta('timestep', timestep)
    snap.writer.hdf_add_meta('num_atoms', reader.meta['num_atoms'])
    snap.writer.hdf_add_meta('xlo', reader.meta['xlo'])
    snap.writer.hdf_add_meta('xhi', reader.meta['xhi'])
    snap.writer.hdf_add_meta('ylo', reader.meta['ylo'])
    snap.writer.hdf_add_meta('yhi', reader.meta['yhi'])
    snap.writer.hdf_add_meta('zlo', reader.meta['zlo'])
    snap.writer.hdf_add_meta('zhi', reader.meta['zhi'])
    snap.writer.hdf_add_meta('col_mapping', col_mapping_no_num)
    snap.writer.hdf_add_meta('vectors', vectors)
    snap.writer.hdf_add_meta('tensors', tensors)
    snap.writer.hdf_add_meta('symm_tensors', symm_tensors)
    snap.writer.hdf_add_meta('sim_name', sim_name)

    # TODO: HARDCODED: 'atoms' 'grid_type'
    snap.writer.hdf_add_meta('grid_type', 'atoms')
    snap.writer.hdf_add_meta('units', units)
    snap.writer.hdf_add_meta('comments', '')
    snap.writer.hdf_add_meta('pbc', pbc)

    snap.writer.flush()
    snap.gather_data()

    # writing default light .xmf data
    snap.writer.xmf_snap_file()

    snap.writer.cleanup()


def main(argv=None):
    """
    Parses user input and converts LAMMPS ASCII dump files to
    XDMF format.

    Returns
    --------
    0
        :func:`main` has completed without problems.
    -1
        :func:`main` encountered a problem and unable to complete.

    """
    if argv is None:
        argv = sys.argv[1:]

    # setting up basic arguments.
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                                help='[default=%default]')
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-f", "--file", action="store", dest="fn_base",
                       help="Search string for files to be processed.",
                       metavar="FNBASE",
                       default=None)
    parser.add_option("-s", "--sim", action="store", dest="sim_name",
                       help="Simulation name, e.g. gridXX.", metavar="SIMNAME",
                       default=None)
    parser.add_option("-c", "--col", action="store", dest="fn_col",
                       help="Name of column definition file. Do not use\
                             for new self-describing LAMMPS dump files.",
                       metavar="COLFILE",
                       default=None)
    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger.procedure_banner('Running dump2xdmf')

    if options.fn_base is None:
        logger.warning('File search string must be supplied.')
        return -1

    flist = util.file_list(options.fn_base, verbose=options.verbose)
    flist = [fn for fn in flist if '.xmf' not in fn]
    flist = [fn for fn in flist if '.h5' not in fn]
    flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists,
                                                          options.verbose)

    logger.user_message('Files to be processed:')
    util.print_file_list(flist)

    if len(flist) < 1:
        logger.completion_msg('No ASCII dump files to process.')
        return -1

    if util.user_approve('Is this list of files correct?') is not True:
        logger.completion_msg('Exiting script.  Please try again.')
        return -1

    col_mapping = {}
    vectors = {}
    tensors = {}
    symm_tensors = {}

    if options.fn_col is None or os.path.exists(options.fn_col) is not True:
        if util.user_approve('New self-describing LAMMPS dump?') is True:
            col_mapping = col.lammps_mapping(flist[0])
        # TODO: Eliminate or require config file
        else:
            good_input = False
            while not good_input:
                fn_col = logger.request('Column mapping file [blank if none]: ',
                                         input_type='raw')
                if fn_col == '':
                    good_input = True
                    col_mapping = col.user_mapping(settings.REQUIRED_DATA)
                elif os.path.exists(fn_col) is True:
                    good_input = True
                    col_mapping = col.file_mapping(fn_col,
                                               required=settings.REQUIRED_DATA)
                else:
                    logger.bad_input('File %s does not exist.' % fn_col)
    # Using options.fn_col
    else:
        col_mapping = col.file_mapping(options.fn_col,
                                               required=settings.REQUIRED_DATA)

    # TODO: Replace with config file based mapping function.
    vectors = col.vector_mapping(col_mapping)
    tensors = col.tensor_mapping(col_mapping)
    symm_tensors = col.symm_tensor_mapping(col_mapping)

    good_input = False
    while not good_input:
        units = logger.request('Enter units system %s: ' % UNITS_ALLOWED,
                                input_type='raw')
        if units in UNITS_ALLOWED:
            good_input = True
        else:
            logger.bad_input('Units given not allowed: %s' % UNITS_ALLOWED)

    good_input = False
    while not good_input:
        try:
            timestep = float(logger.request("Enter the timestep: ",
                                                            input_type='raw'))
        except ValueError as exc:
            print('Timestep - ValueError caught: %s.' % exc)
            continue
        else:
            good_input = True

    if options.sim_name is None:
        good_input = False
        while not good_input:
            sim_name = logger.request("Enter simulation name: ",
                                                            input_type='raw')
            if (' ' not in sim_name or '_' not in sim_name) is True:
                logger.bad_input('Cannot have spaces or underscores.')
            else:
                good_input = True
    else:
        # Removing ' ' and '_' from options.sim_name
        sim_name = ''.join(options.sim_name.split())
        sim_name = ''.join(options.sim_name.split('_'))

    good_input = False
    while not good_input:
        pbc = logger.request('Enter PBC [list of bool]: ', input_type='input')
        pbc = np.array(pbc, dtype=np.bool)
        if pbc.shape != (3,):
            logger.bad_input('List must be shape (3x1).')
        else:
            good_input = True

    # Mass is not one of columns in dump file.
    if 'mass' not in col_mapping:
        good_input = False
        while not good_input:
            mass = logger.request('Enter of masses of atom types [list]: ',
                                  input_type='input')
            if (isinstance(mass, list) or isinstance(mass, tuple)) is not True:
                logger.bad_input('Input must be a Python list or tuple.')
            else:
                try:
                    mass = [float(i) for i in mass]
                except ValueError as exc:
                    logger.bad_input('Masses - ValueError caught: %s.' % exc)
                    continue
                good_input = True
    else:
        mass = None

    for fn in flist:
        logger.user_message('Processing: %s' % fn)

        # spawning secondary process done to aid in memory management.
        p = Process(target=process_single,
                    args=(fn, col_mapping, vectors, tensors, symm_tensors,
                          sim_name, units, timestep, mass, pbc)
                    )
        p.start()
        p.join()

    logger.completion_msg('Finished processing dumps: %s' % sim_name)
    return 0

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
