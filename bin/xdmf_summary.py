#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Creates temporal summary file for a collection of XDMF files.

The temporal summary file is a XDMF formatted XML file (``.xmf``)
which adds individual snapshot ``.xmf`` files into a XDMF temporal
collection.

"""

import sys
import os
from optparse import OptionParser

import mdhandle.utilities as util
from mdhandle.writers.xdmf import xmf_write_temporal_summary
from mdhandle.logger import Logger
from mdhandle.simcontainer import SimContainer

#===============================================================================

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    # setting up basic arguments.
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
            help='[default=%default]')
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-f", "--file", action="store", dest="fn_base",
                       help="File basename to be processed.",
                       metavar="FNBASE",
                       default=None)
    parser.add_option("-d", "--dataset", action="store", dest="data_set",
                       help="Dataset name, [default=rawSimResults]",
                       metavar="DATASET",
                       default='rawSimResults')
    parser.add_option("-t", "--time", action="store", dest="time",
                       help="time_lo time_hi t_interval - inclusive.\n\
                             set t_interval to 1 to select all.",
                       metavar="TIME",
                       type="int", nargs=3, default=None)
    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger(options.verbose)
    logger.procedure_banner('XDMF - Temporal Summary Creator')

    if options.fn_base is None:
        logger.warning('File search string must be given')
        return -1

    flist = util.file_list(options.fn_base, verbose=options.verbose)
    flist = [fn for fn in flist if '.h5' in fn]
    flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists,
                                                          options.verbose)

    if options.time is not None:
        if len(options.time) != 3:
            logger.warning('Time values provide must be t_lo t_hi t_interval')
            logger.warning('Set t_interval to 1 to select all on [t_lo, t_hi]')
            logger.completion_msg('Exiting')
            return -1

        if not (options.time[1] > options.time[0]):
            logger.warning('Time values must be in order t_lo t_hi t_interval')
            logger.completion_msg('Exiting')
            return -1

        if options.time[2] == 1:
            # Using all timesteps that matched the glob
            flist_culled = []
        else:
            # Produced list of times - inclusive of time_hi
            time_list = range(options.time[0], options.time[1]+options.time[2],
                                               options.time[2])
            flist, flist_culled = util.cull_flist_by_time(flist, time_list, 
                                                               options.verbose)

    logger.user_message('Files to be processed for dataset=%s:' % 
                                                            options.data_set)

    util.print_file_list(flist)
    
    if options.verbose is True:
        logger.user_message('Number of files culled from glob: %s' 
                            % len(flist_culled) )

    if util.user_approve('Is this list of files correct?') is not True:
        logger.completion_msg('Exiting script.  Please try again.')
        return -1

    sim_cont = SimContainer(flist, dataset_name=options.data_set,
                                   verbose=options.verbose)
    xmf_write_temporal_summary(sim_cont)

    logger.completion_msg('All done - XMDF generated: %s' %
                                                    sim_cont.meta['sim_name'])
    return 0

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
