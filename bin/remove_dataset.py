#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Removes dataset from snapshot file.

Cannot be applied to default dataset :attr:`settings.DEFAULT_DATASET`.

Returns ``0`` if successful, otherwise ``-1``.

.. warning::
   Destructively removes dataset and cannot be undone.

"""

            # TODO: removing XMF file associated with dataset.

import sys
import os
from optparse import OptionParser

from mdhandle import utilities as util
from mdhandle.snap import Snap
from mdhandle.logger import Logger
from mdhandle import settings

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
    parser.add_option("-d", "--dataset", action="store", dest="to_delete",
                        help="Name dataset to delete",
                        metavar="DSET",
                        default=None)

    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger()
    logger.procedure_banner('Fixing Metadata - Velocity')

    to_delete = options.to_delete
    if to_delete is None:
        logger.completion_msg('Dataset to delete not given.')
        return -1
    
    if to_delete == settings.DEFAULT_DATASET:
        logger.warning("Cannot delete %s dataset." % settings.DEFAULT_DATASET)
        return -1

    if options.fn_base is None:
        logger.warning('Search string for files must be given.')

    flist = util.file_list(options.fn_base, verbose=options.verbose)
    flist = [fn for fn in flist if '.h5' in fn]
    flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists,
                                                          options.verbose)
    logger.user_message('Files culled from original flist: ')
    util.print_file_list(fn_culled)

    util.print_file_list(flist, 'the list')
    if util.user_approve('Is this list of files correct?') is not True:
        logger.completion_msg('Exiting script.  Please try again.')
        return -1

    for fn in flist:
        logger.user_message('Processing: %s' % fn)
        snap = Snap(fn, dataset_name='rawSimResults')
        snap.gather_data()
        if to_delete not in snap.get_datasets():
            logger.completion_msg('%s not dataset in snapshot' % to_delete)
            return -1

        else:
            snap.writer.hdf_delete_dataset(to_delete)
            snap.writer.flush()
            snap.writer.cleanup()

    logger.user_message('Checking all files are ok after dataset deletion')
    for fn in flist:
        snap = Snap(fn, dataset_name='rawSimResults')
        snap.gather_data()
        logger.user_message('Datasets in snap: %s' % snap.get_datasets())
        snap.writer.cleanup()

    logger.completion_msg('All Done')
    return 0

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
