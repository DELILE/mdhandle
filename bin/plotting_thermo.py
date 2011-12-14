#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Loads log file supplied at command line for thermo plotting via
:mod:`matplotlib`.

Log file handle is stored in ``thermo`` in top namespace.

Returns log file object (``thermo``), or ``-1`` if problem arises.

"""

import sys
from optparse import OptionParser

import numpy as np
import matplotlib.pyplot as plt

from mdhandle.readers.log import Log
from mdhandle.logger import Logger
import mdhandle.interactive as interactive

# -----------------------------------------------------------------------------

# If not in IPython - lauch IPython.
interactive.launch_ipython()

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    # setting up basic arguments.
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
            help='[default=%default]')
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-f", "--file", action="store", dest="fn_log",
                       help="Log file name. Quote string if using wildcards.",
                       metavar="FNLOG",
                       default=None)
    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger(options.verbose)
    logger.procedure_banner('Loading LAMMPS log file')
    
    if options.fn_log is None:
        logger.warning('Name of log file must be given.')
        return -1

    thermo = Log(options.fn_log)
    return thermo


#------------------------------------------------------------------------------


if __name__ == "__main__":
    thermo = main()
