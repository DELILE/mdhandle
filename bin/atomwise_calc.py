#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Calculating atomwise quantities from existing simulation results.

Returns ``0`` if successful, ``-1`` otherwise.

.. note::
   Some accuracy will be lost due to floating point rounding, so cannot expect
   exact agreement with data produced via simulation it self.

"""

#----------------------------------------------------------------------------

import os
import sys
from optparse import OptionParser

from mdhandle.snap import Snap
import mdhandle.utilities as util
import mdhandle.properties.atom_properties as atom_props
from mdhandle.logger import Logger

# TODO: handle case where target array already exists.
#       eg. tables.exceptions.NodeError: group ``/rawSimResults`` already has a
#       child node named ``epair``

# ----------------------------------------------------------------------------


def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # setting up basic arguments.
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
            help='[default=%default]')
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-f", "--files", action="store", dest="fn_base",
                       help="Search string for files to be processed.",
                       metavar="FNBASE",
                       default=None)
    parser.add_option("-c", "--calcs", action="append", dest="calcs",
                      help='Atomwise Calculations to perform, space separated.',
                      metavar="CALCS",
                      default=None)
    parser.add_option("-d", "--dataset", action="store", dest="dataset_name",
                      help="Name of dataset.", metavar="DSET",
                      default='rawSimResults')

    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    # posix reliant call
    logger = Logger()
    logger.procedure_banner('Recovering: Epair, Force, Coord, Stress')

    if options.fn_base is None:
        logger.warning('Search string for files must be given.')
        return -1

    if options.calcs is None:
        good_input = False
        while good_input is False:
            calcs = logger.request('Enter Python list of calculation names: ',
                                    input_type='input')
            if isinstance(calcs, list) is True:
                good_input = True
            else:
                logger.bad_input('Calculations must be a list of strings.')
    else:
        calcs = options.calcs

    flist = util.file_list(options.fn_base, verbose=options.verbose)
    flist = util.cull_flist_by_function(flist, os.path.exists,
                                                    verbose=options.verbose )
    flist = [fn for fn in flist if '.h5' in fn]

    if options.verbose:
        util.print_file_list(flist)

        if not util.user_approve('Is this list of files ok?'):
            logger.completion_msg('Ok - Exiting recover_atomwise.py')
            return -1

    for fn in flist:
        logger.user_message('Processing: %s' % fn)
        snap = Snap(fn, dataset_name='rawSimResults')
        snap.gather_data()

        # NOTE: xyz positions automatically wrapped to central image right now
        # TODO: take cutoff from the snap rather than hardcoded (i.e. 4.0).
        calc = atom_props.NonLocal(snap, 4.0, calcs, wrap=True)
        calc.setup_calculations()
        calc.run()
        
        for name in calc.to_calc:
            snap.writer.hdf_create_array(name, calc.results[name],
                                               dtype=calc.dtypes[name],
                                               location=calc.locations[name],
                                               update_meta=True)

            if calc.to_calc['type'] == 'scalars':
                pass
            elif calc.to_calc['type'] in ['vectors', 'tensors', 'symm_tensors']:
                snap.writer.hdf_update_composite_mapping(name, 'vectors', 
                                                 dtype=calc.dtype[name],
                                                 location=calc.locations[name],
                                                 comps=None)

        snap.writer.flush()
        snap.gather_data()

        # writing light .xmf data
        snap.writer.xmf_snap_file()
        snap.writer.cleanup()

    logger.completion_msg('Finished recover_atomwise.py for %s' %
                                                        snap.meta['sim_name'])

#-----------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
