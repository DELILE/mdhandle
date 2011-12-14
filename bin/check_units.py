#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Check consistency of units metadata is constant in set of snapshot files.

Return ``0`` if successful, otherwise ``-1`` for error and ``1`` if units not
consistent.


"""

import sys
import os
from optparse import OptionParser


import mdhandle
import mdhandle.utilities as util
from mdhandle.snap import Snap
from mdhandle.logger import Logger

# -----------------------------------------------------------------------------


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    # setting up basic arguments.
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage, version='%prog '+mdhandle.__version__)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                                         help='[default=%default]')
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-f", "--file", action="store", dest="fn_base",
                       help="Search string for files to be processed.",
                       metavar="FNBASE",
                       default=None)
    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger(options.verbose)
    logger.procedure_banner('Checking PBC consistency')

    if options.fn_base is None:
        logger.warning('File search string must be supplied.')
        return -1

    flist = util.file_list(options.fn_base, verbose=options.verbose)
    flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists,
                                                          options.verbose)
    logger.user_message('List of files to be processed: ')
    util.print_file_list(flist)

    if util.user_approve('Is this list of files correct?') is not True:
        logger.completion_msg('Exiting script.  Please try again.')
        return -1

    units_in_datasets = {}
    first_time = True
    units_constant = True
    for fn in flist:
        logger.user_message('Processing --> %s' % fn)

        snap = Snap(fn, dataset_name='rawSimResults', verbose=options.verbose)
        dataset_list = snap.get_datasets()

        for dataset_name in dataset_list:
            snap.set_active_dataset(dataset_name, reset_meta=True)
            snap.gather_data()
            logger.user_message('Dataset: %s --> %s' %
                                          (dataset_name, snap.meta['units']))

            if first_time is True:
                units_in_datasets[dataset_name] = snap.meta['units']
            else:
                if snap.meta['units'] != units_in_datasets[dataset_name]:
                    logger.warning('Units change %s.%s' % (fn, dataset_name))
                    units_constant = False

        first_time = False
        snap.writer.cleanup()

    logger.completion_msg('Done checking units')
    if units_constant is True:
        return 0
    else:
        return 1

#------------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
