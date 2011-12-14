#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Regenerate `.xmf` files for dataset within provided snapshot files.

"""

import sys
import os
from optparse import OptionParser

import mdhandle
from mdhandle import utilities as util
from mdhandle.snap import Snap
from mdhandle.logger import Logger

#------------------------------------------------------------------------------


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
                       help="Search string for HDF files to be processed.",
                       metavar="FNBASE",
                       default=None)
    parser.add_option("-d", "--dataset", action="store", dest="dataset_name",
                       help="Name of dataset.",
                       metavar="DSET",
                       default='rawSimResults')
    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger(options.verbose)
    logger.procedure_banner('Regenerating .xmf files: %s' % options.fn_base)

    if options.fn_base is None:
        logger.warning('Search string for files must be given')
        return -1

    flist = util.file_list(options.fn_base, verbose=options.verbose)
    flist = [fn for fn in flist if '.h5' in fn]
    flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists,
                                                          options.verbose)

    logger.user_message('Files to be processed:')
    util.print_file_list(flist)

    if util.user_approve('Is this list of files correct?') is not True:
        logger.completion_msg('Exiting script.  Please try again.')
        return -1

    # -------------------------------------------------------------------
    # Add any further modifications to snaps as necessary here or in loop.
    # -------------------------------------------------------------------

    for fn in flist:
        snap = Snap(fn, dataset_name=options.dataset_name)

        # Do further actions on individual snap file.
        snap.writer.xmf_snap_file()
        snap.writer.cleanup()

    logger.completion_msg('Done rewriting .xmf files: %s' % options.fn_base)
    return 0


if __name__ == "__main__":
    sys.exit(main())
