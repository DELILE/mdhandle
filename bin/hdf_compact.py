#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Compacts HDF files.

Important when HDF files have been repeatedly modified.

Returns ``0`` if successful, otherwise ``-1``.

"""

import sys
import os
from optparse import OptionParser

import mdhandle.utilities as util
from mdhandle.writers.xdmf import hdf_compact
from mdhandle.logger import Logger

#------------------------------------------------------------------------------


def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # setting up command line argument parsing
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                                         help='[default=%default]')
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-f", "--files", action="store", dest="fn_base",
                       help="File base name search string.",
                       metavar="FNBASE",
                       default=None)
    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger(options.verbose)

    logger.procedure_banner('Compacting HDF file script')

    if options.fn_base is None:
        logger.warning('Search string for file to be process is needed.')
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

    hdf_compact(flist, verbose=options.verbose)

    logger.completion_msg('ALL DONE HDF COMPACTION')
    return 0

#------------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
