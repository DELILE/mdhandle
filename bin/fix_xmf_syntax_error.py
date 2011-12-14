#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Fix small syntax error in ``.xmf`` file syntax for set of snapshots.

Fixing is done by simple replacement of offending string.

Returns ``0`` if successful, otherwise ``-1``.

.. note::
   Alternative is to simply ask for new ``.xmf`` file to be generated
   based on the data contained in the snap file.

"""
#-----------------------------------------------------------------------------

import sys
import os
from optparse import OptionParser

import mdhandle.utilities as util
from mdhandle.logger import Logger

# -----------------------------------------------------------------------------


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
    parser.add_option("-t", "--text", action="store", dest="fixes",
                        help="bad_string new_string",
                        type="string", nargs=2, metavar="TEXT", default=None)
    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger(options.verbose)
    logger.procedure_banner('Fixing xml syntax script: %s --> %s' %
                                        (options.fixes[0], options.fixes[1]))
    logger.warning('Alternative to regenerating .xmf from snap file data')

    if options.fn_base is None:
        logger.warning('Search string for files is needed.')
        return -1

    if options.fixes is None:
        logger.warning('Nothing to do - offending and replacement text\
                                                        not given')

    flist = util.file_list(options.fn_base, verbose=options.verbose)
    flist = [fn for fn in flist if '.xmf' in fn]
    flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists,
                                                          options.verbose)

    logger.user_message('Files to be processed:')
    util.print_file_list(flist)

    if util.user_approve('Is this list of files correct?') is not True:
        logger.completion_msg('Exiting script.  Please try again.')
        return -1

    for fn in flist:
        logger.user_message('Processing --> %s' % fn)
        f_old = open(fn, 'r')

        # files are small enough to read into long string.
        contents = f_old.read()
        f_old.close()

        fixes_strings = [str(i) for i in options.fixes]
        contents = contents.replace(fixes_strings[0], fixes_strings[1])

        # Replacing old file
        f_new = open(fn, 'w')
        f_new.write(contents)
        f_new.close()

    logger.completion_msg('Done fixing XML syntax errors')
    return 0

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(main())
