#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Make units metadata the same thoughout set of snapshots.

Returns ``0`` if successful, otherwise ``-1``.

"""

import sys
import os
from optparse import OptionParser


import mdhandle.utilities as util
from mdhandle.snap import Snap
from mdhandle import units
from mdhandle.logger import Logger

# ----------------------------------------------------------------------------


def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # setting up command line argument parsing
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                                         help='[default=%default]')
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-f", "--file", action="store", dest="fn_base",
                       help="Search string for files to be processed.",
                       metavar="FILE",
                       type='string', default=None)
    parser.add_option("-u", "--units", action="store", dest="units",
                        help="Units: {'lj', 'metal', 'si'}",
                        type="string", metavar="UNITS", default=None)
    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger(options.verbose)
    logger.procedure_banner('Setting units')

    if options.fn_base is None:
        logger.warning('File search string must be supplied.')
        return -1

    flist = util.file_list(options.fn_base, verbose=options.verbose)
    flist = [fn for fn in flist if '.xmf' in fn]
    flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists,
                                                          options.verbose)

    logger.user_message('Files to be processed:')
    util.print_file_list(flist)

    if util.user_approve('Is this list of files correct?') is not True:
        logger.completion_msg('Exiting script.  Please try again.')
        return -1

    if options.units not in units.UNITS_ALLOWED or options.units is None:
        logger.completion_msg('Units not in allowed: %s' % units.UNITS_ALLOWED)
        return -1

    for fn in flist:
        snap = Snap(fn, 'rawSimResults')
        snap.gather_data()

        logger.user_message('Changing %s: %s --> %s' % (fn, snap.meta['units'],
                                                            options.units))

        snap.writer.hdf_modify_meta('units', units)

        snap.writer.cleanup()

    logger.completion_msg("Done modifying 'units' metadata")
    return 0

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(main())
