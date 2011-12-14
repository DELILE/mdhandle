# !/usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2

"""
Loads LAMMPS log file and plots temperature, kinetic energy, potential energy
and total energy versus intrinsic array length.

Useful as diagnositc script for LAMMPS simulations.

After initial plots are complete - user can continue plotting in 
interactive namespace.

"""

#----------------------------------------------------------------------------

import sys
from optparse import OptionParser

import matplotlib.pyplot as plt

from mdhandle.readers.log import Log
import mdhandle.logger as logger
import mdhandle.utilities as util
import mdhandle.interactive as interactive

# ----------------------------------------------------------------------------

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
    parser.add_option("-f", "--file", action="store" ,dest="fn_log",
                       help='Log file name. Quote string if using wildcards.',
                       metavar="FNLOG",
                       default=None)

    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    log = logger.Logger()
    log.procedure_banner('Testing Thermodynamic Output - %s' % options.fn_log)

    if options.fn_log is None:
        logger.warning('Log file name must be given.')
        return -1

    thermo = Log(options.fn_log)

    # ASSUMPTION: Assume that LAMMPS log columns have default names.
    temp = thermo.get_vec('Temp')
    epair = thermo.get_vec('PotEng')
    ke = thermo.get_vec('KinEng')
    toteng = thermo.get_vec('TotEng')

    for (data, title) in [(temp, 'temp'), (epair, 'epair'), 
                          (ke,'ke'), (toteng, 'toteng')]:
        plt.figure()
        plt.plot(data)
        plt.title(title)

# ----------------------------------------------------------------------------

if __name__ == "__main__":
    return_code = main()
    
    if util.user_approve('Quit? ') is True:
        sys.exit(return_code)
    else:
        logger.user_message('Continuing with interactive plotting.')
