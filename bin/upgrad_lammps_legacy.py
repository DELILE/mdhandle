#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Conerting legacy (pre ~2007) LAMMPS ASCII dump files to self-describing format.

The identity of each column in the dump file is specified in line 9 of the
file header.

1   ITEM: TIMESTEP
2   6000000
3   ITEM: NUMBER OF ATOMS
4   208852
5   ITEM: BOX BOUNDS
6   0 392
7   0 795.68
8   0 392
9   ITEM: ATOMS id type x y z ix iy iz xu yu zu vx vy vz fx fy fz

.. note::
   Self-describing format is needed to be read by recent versions of VMD.

"""

import sys
import os

from optparse import OptionParser

import mdhandle
import mdhandle.utilities as util
from mdhandle.logger import Logger
from mdhandle.readers.lammps_dump import LAMMPS_dump

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
                       help="Search string for files to be processed",
                       metavar="FNBASE",
                       default=None)
    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger(options.verbose)
    logger.procedure_banner('Legacy LAMMPS dumps to Self-Describing Format')

    if options.fn_base is None:
        logger.warning('File searchs string must be provided')
        return -1

    flist = util.file_list(options.fn_base, verbose=options.verbose)
    flist = [fn for fn in flist if '.xmf' not in fn]
    flist = [fn for fn in flist if '.h5' not in fn]
    flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists,
                                                          options.verbose)

    logger.user_message('Files to be processed: ')
    util.print_file_list(flist)

    if util.user_approve('Is this list of files correct?') is not True:
        logger.completion_msg('Exiting script.  Please try again.')
        return -1

    cols = logger.request('Enter ordered list of columns [space separated]: ',
                          input_type='raw')

    for fn in flist:
        dump_reader = LAMMPS_dump(fn)
        dump_reader.read_header()

        if dump_reader.is_self_describing is True:
            logger.warning('LAMMPS dump %s is already self-describing' % fn)
            continue

        if len(cols) != dump_reader.meta['ncol']:
            logger.warning('Number of columns does not agree number of\
                            columns in ASCII dump file')
            return -1

        logger.user_message('Processing: %s' % fn)
        fn_tmp = os.path.join(os.path.dirname(fn), 'tmp_dump_file.tmp')

        ret_code = os.system("sed -e 's,^ITEM: ATOMS,ITEM: ATOMS %s,' %s > %s"
                                                        % (cols, fn, fn_tmp) )
        if ret_code == 0:
            os.rename(fn_tmp, fn)
        else:
            if options.verbose is True:
                logger.warning('Failed: %s' % os.path.basename(fn))
            os.remove(fn_tmp)

    logger.completion_msg('Done upgrading LAMMPS dumps to self-describing\
                           format')
    return 0


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
