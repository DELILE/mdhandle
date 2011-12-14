#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Reads a LAMMPS dump file, takes off first 9 lines of file and saves the
remaining contents of the file in sub-directory, ``output`` of the
current working directory'

Returns ``0`` if successful, otherwise ``-1``.

"""

import sys
import os
from optparse import OptionParser

from mdhandle.settings import LAMMPSHEADER
from mdhandle.logger import Logger
from mdhandle import utilities as util

# -----------------------------------------------------------------------------


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
                       help="Search string for HDF files to be processed.",
                       metavar="FNBASE",
                       default=None)
    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger()
    logger.procedure_banner('Fixing Metadata - PBC')

    if options.fn_base is None:
        logger.warning('File search string must be given.')
        return -1

    fn_base = options.fn_base

    flist = util.file_list(fn_base, verbose=options.verbose)
    flist = [fn for fn in flist if '.h5' in fn]
    flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists,
                                                          options.verbose)

    if len(fn_culled) > 0:
        logger.user_message('Files removed because of improper format:')
        util.print_file_list(fn_culled)

    util.print_file_list(flist, 'the list')
    if util.user_approve('Is this list of files correct?') is not True:
        logger.completion_msg('Exiting script.  Please try again.')
        return -1

    # creating output directory if it doesn't exist

    if (os.path.exists('output') is not True and
        os.path.isdir('output') is not True):
        os.mkdir('ouput')
    else:
        logger.warning("./output directory already exists")
        logger.completion_msg('Exiting')
        return -1

    for fn in flist:
        fn_new = os.path.join('output', os.path.basename(fn))
        logger.user_message("Processing: %s" % fn)
        cmd = "sed '1,%dd' %s > %s" % (fn, len(LAMMPSHEADER), fn_new)
        os.system(cmd)

    logger.completion_msg('Done removing LAMMPS dump header')
    return 0

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
