#!/usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Calculates spatially averaged binned quantities from atomwise snap data.
Calculation is done on a regular rectilinear grid.

Grid: :class:`mdhandle.properties.binned_properties.RectGrid`

Setup with default configurations for gridding applied to homogeneous
fluidic system.

Returns ``0`` if successful, ``-1`` otherwise.
"""

import sys
from optparse import OptionParser

import numpy as np

from mdhandle.snap import Snap
from mdhandle.properties import binned_properties
from mdhandle import utilities as util
from mdhandle.logger import Logger
from mdhandle.writers.xdmf import valid_file

#------------------------------------------------------------------------------


def main(argv=None):

    # User Interaction and Command Line options
    if argv is None:
        argv = sys.argv[1:]

    # setting up basic arguments.
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
            help='[default=%default]')
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-f", "--file", action="store", dest="fn_base",
                       help="Search string for siles to be processed.",
                       metavar="FILES",
                       default=None)
    parser.add_option("-c", "--calcs", action="append", dest="calcs",
                      help='Name of calculations to perform.  Repeat if req.',
                      metavar="CALCS",
                      default=None)
    parser.add_option("-s", "--spacing", action="store", dest="spacing",
                      help="Grid spacing: dx dy dz", metavar="SPACE",
                      type='float', nargs=3, default=None)
    parser.add_option("-n", "--name", action="store", dest="new_dataset",
                      help="Name of new binned dataset.", metavar="NAME",
                      default=None)

    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger()
    logger.procedure_banner('Running Binning Calculation')
    
    logger.warning('This version of binned_calc.py only for fluid only systems')

    if options.fn_base is None:
        logger.warning('File search string needed.')
        return -1

    flist = util.file_list(options.fn_base, verbose=options.verbose)
    flist = [fn for fn in flist if '.h5' in fn]
    flist, fn_culled = util.cull_flist_by_function(flist, valid_file,
                                                          options.verbose)

    util.print_file_list(flist, 'the list')

    if util.user_approve('Is this list of files correct?') is not True:
        logger.completion_msg('Exiting script.  Please try again.')
        return -1


    if options.spacing is None:
        good_input = False
        while good_input is False:
            spacing = logger.request('Enter Python list for grid spacing: ',
                                     input_type='input')
            if isinstance(spacing, list) is True:
                spacing = np.array(spacing, dtype=np.float)
                good_input = True
            else:
                logger.bad_input('Spacing must be in Python list of floats.')
    else:
        spacing = np.array(options.spacing, dtype=np.float)

    if options.new_dataset is None:
        good_input = False
        while good_input is False:
            new_dataset = logger.request('Name of new dataset for grid data: ',
                                     input_type='raw')
            if (' ' not in new_dataset) and ('_' not in new_dataset):
                good_input = True
            else:
                logger.bad_input("' ' and '_' cannot be present in \
                                                                dataset name")
    else:
        # Eliminating spaces and underscores in new_dataset from cmd line
        new_dataset = ''.join(options.new_dataset.split())
        new_dataset = ''.join(new_dataset.split('_'))

    if options.calcs is None:
        good_input = False
        while good_input is False:
            calcs = logger.request('Enter Python list of grid calc names: ',
                                     input_type='input')
            if isinstance(calcs, list) is True:
                good_input = True
            else:
                logger.bad_input('Submitted calc names must be contained\
                                  in a list')
    else:
        calcs = options.calcs


    for fn in flist:
        snap = Snap(fn, dataset_name='rawSimResults')
        snap.gather_data()

        grid = binned_properties.RectGrid(snap, snap.sim_cell._boxL,
                                                spacing,
                                                origin=snap.sim_cell.origin,
                                                offset=True,
                                                mapping='cic',
                                                calcs=calcs,
                                                wrap=True,
                                                selection=True,
                                                location='Node')
        grid.setup_calculations()
        grid.run_mapping()
        grid.run_calc()

        # Storing metadata from atomwise snap as blanked out when change
        # active dataset
        time = snap.meta['time']
        timestep = snap.meta['timestep']
        sim_name = snap.meta['sim_name']
        units = snap.meta['sim_name']

        snap.set_active_dataset(new_dataset, create_new=True)

        # adding metadata
        snap.writer.hdf_add_meta('time', time)
        snap.writer.hdf_add_meta('timestep', timestep)
        snap.writer.hdf_add_meta('xlo', grid.origin[0])
        snap.writer.hdf_add_meta('xhi', (grid.origin + grid._boxL)[0])
        snap.writer.hdf_add_meta('ylo', grid.origin[1])
        snap.writer.hdf_add_meta('yhi', (grid._boxL + grid.origin)[1])
        snap.writer.hdf_add_meta('zlo', grid.origin[2])
        snap.writer.hdf_add_meta('zhi', (grid._boxL + grid.origin)[2])
        snap.writer.hdf_add_meta('sim_name', sim_name)

        # HARDCODED: grid_type as '3DRegular'
        snap.writer.hdf_add_meta('grid_type', '3DRegular')
        snap.writer.hdf_add_meta('units', units)
        snap.writer.hdf_add_meta('origin', grid.origin)
        snap.writer.hdf_add_meta('boxL', grid._boxL)
        snap.writer.hdf_add_meta('grid_spacing', spacing)
        snap.writer.hdf_add_meta('num_elem', grid.num_elem)
        snap.writer.hdf_add_meta('pbc', snap.sim_cell.pbc)
        snap.writer.hdf_add_meta('num_atoms', grid.num_atoms)
        snap.writer.hdf_add_meta('comments', '')

        snap.writer.hdf_add_meta('coverage', grid.coverage)
        snap.writer.hdf_add_meta('col_mapping', grid.to_calc)
        snap.writer.hdf_add_meta('vectors', grid.vectors)
        snap.writer.hdf_add_meta('tensors', grid.tensors)
        snap.writer.hdf_add_meta('symm_tensors', grid.symm_tensors)

        for name in grid.to_calc:
            snap.writer.hdf_create_array(name, grid.results[name],
                                              dtype=grid.dtypes[name],
                                              location=grid.locations[name],
                                              update_meta=False)
        snap.writer.flush()
        snap.gather_data()

        # writing light .xmf data
        snap.writer.xmf_snap_file()
        snap.writer.cleanup()

    logger.completion_msg('Done gridding: %s' % options.fn_base)
    return 0

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
