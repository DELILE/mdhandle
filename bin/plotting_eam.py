#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Plotting EAM function components using :class:`mdhandle.readers.eam.EAM`
object.

"""

import sys
from optparse import OptionParser

import matplotlib.pyplot as plt

from mdhandle.readers.eam import EAM
from mdhandle.logger import Logger
from mdhandle import interactive

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
    parser.add_option("-f", "--file", action="store", dest="fn",
                       help="EAM filename.", metavar="FILE", default=None)
    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    fn = options.fn
    verbose = options.verbose

    logger = Logger(verbose)
    logger.procedure_banner('Plotting EAM Potential file: %s' % fn)

    if options.fn is None:
        logger.warning('EAM potential file must be given.')
        return -1

    plt.close('all')

    eam = EAM(fn, verbose)
    eam.plot()

    logger.completion_msg('Done Plotting - Going to Interactive Prompt')
    return eam

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    interactive.launch_ipython()

    eam = main()
