#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Looks through set of HDF files and reports if set of properties are found
in the snapshot. 

Default properties are : 'epair', 'force', 'coord' and 'stress'.

Results reported to screen for those files which are missing desired
atomwise quantities.

"""

import sys
import os
from optparse import OptionParser

import mdhandle
import mdhandle.utilities as util
from mdhandle.snap import Snap
from mdhandle.logger import Logger

# ----------------------------------------------------------------------------


def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # setting up command line arguments
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage, version='%prog '+mdhandle.__version__)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                                         help='[default=%default]')
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-f", "--fnBase", action="store", dest="fn_base",
                       help="File basename to be processed.",
                       metavar="FNBASE",
                       default=None)
    parser.add_option('-d', '--dataset', action='store', dest='dataset_name',
                            metavar='DSET', default='rawSimResults',
                            help='Snapshot dataset.')
    parser.add_option('-n', '--names', action='append', dest='properties',
                            metavar='PROPS',
                            default=['epair', 'force', 'stress', 'coord'],
                            help='Properties to look for in HDF file')
    (options, args) = parser.parse_args(argv)

    logger = Logger(options.verbose)
    logger.procedure_banner('Looking through snapshots for scalar, vector \
                                                        or tensor quantities')

    if options.fn_base is None:
        logger.warning('Search string for file must be given.')
        return -1


    flist = util.file_list(options.fn_base)
    flist = [fn for fn in flist if '.h5' in fn]
    flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists, 
                                                    options.verbose)

    for fn in flist:
        logger.user_message('Processing: %s.%s' % (fn, options.dataset_name))

        snap = Snap(fn, options.dataset_name, verbose=options.verbose)
        snap.gather_data()

        for prop in options.properties:
            if ( prop not in snap.col_mapping and
                 prop not in snap.vectors    and
                 prop not in snap.tensors    and
                 prop not in snap.symm_tensors):
                logger.warning('%s not in %s' % (prop, fn))

        snap.writer.cleanup()

    logger.completion_msg('Finished testing HDF files for %s' %
                                                            options.properties)
    return 0

if __name__ == '__main__':
    sys.exit(main())
