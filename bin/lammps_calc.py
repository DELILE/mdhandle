#! /usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Run a single timestep of LAMMPS solid simulations to recover coord and stress.

The name of the data file for LAMMPS initialization 
is 'lammps_data_stress_calc.data'.  This datafile name should be used in the
LAMMPS input file.

Return ``0`` if successful, ``-1`` otherwise.

"""

import os
import sys
from optparse import OptionParser

from mdhandle.snap import Snap
import mdhandle.utilities as util
import mdhandle.writers.lammps_data as lammps_data
import mdhandle.properties.atom_properties as atom_props
from mdhandle.logger import Logger
from mdhandle.settings import TMP_LOCATION

# ----------------------------------------------------------------------------

tmp_data_file_name = os.join(TMP_LOCATION, 'lammps_data_stress_calc.data')
default_calcs = ['stress', 'coord', 'epair', 'force']

def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # setting up command line arguments.
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
            help='[default=%default]')
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
    parser.add_option("-f", "--filebase", action="store", dest="fn_base",
                       help="File basename to be processed.",
                       metavar="FNBASE",
                       default=None)
    parser.add_option("-t", "--time", action="store", dest="time",
                       help="time_lo time_hi time_interval", 
                       metavar="TIME",
                       type="int", nargs=3, default=None)
    parser.add_option("-p", "--props", action="store",dest="prop",
                       help="Space separated atomwise properties.", 
                       metavar="PROPS",
                       default=default_calcs)
    parser.add_option("-l", "--lammps", action="store" ,dest="lammps_in",
                       help="LAMMPS input file",
                       metavar="LAMMPS",
                       default=None)

    parser.set_defaults(verbose=False)
    (options, args) = parser.parse_args(argv)

    logger = Logger()
    logger.procedure_banner('Recovering: coord and stress and epair')

    if options.fn_base is None:
        logger.warning('Search string for files must be given.')
        return -1

    if options.lammps_in is None:
        logger.warning('LAMMPS input file must be given.')
        return -1

    # --------------------- Processing CMD LINE ARGS --------------------------

    if len(options.time) != 3:
        logger.warning('Time values provide must be t_lo t_hi')
        logger.completion_msg('Exiting')
        return -1

    time_list = range(options.time[0], options.time[1], options.time[2])

    if not (options.time[1] > options.time[0]):
        logger.warning('Time values must be in order t_lo t_hi t_interval')
        logger.completion_msg('Exiting')
        return -1

    flist = util.file_list(options.fn_base, verbose=options.verbose)
    flist = [fn for fn in flist if '.h5' in fn]
    flist, fn_culled = util.cull_flist_by_function(flist, os.path.exists,
                                                          options.verbose)
    flist, fn_culled2 = util.cull_flist_by_time(flist, time_list, 
                                                       options.verbose)

    logger.user_message('Files to be processed:')
    util.print_file_list(flist)

    if util.user_approve('Is this list of files correct?') is not True:
        logger.completion_msg('Exiting script.  Please try again.')
        return -1

    if options.prop == default_calcs:
        props_list = options.prop
    else:
        props_list = options.prop.split()
    logger.user_message('Calculating: %s' % props_list) 

    if flist == []:
        logger.completion_msg('No files to process. Exiting.')
        return 0

    # ----------------------------- PROCESSING FILES --------------------------

    for fn in flist:
        if options.verbose is True:
            logger.user_message('Processing: %s' % fn)
        
        snap = Snap(fn,'rawSimResults', verbose=options.verbose)
        snap.gather_data()

        making_data = lammps_data.LAMMPS_Data(snap)
        making_data.new(tmp_data_file_name)
        
        # Padding is required in non-periodic directions because LAMMPS 
        # produced error if atom sitting right at the boundary
        padding = snap.sim_cell.pbc*(0.000001)
        making_data.header(non_pbc_pad=padding)
        making_data.body_atoms(wrap=False)
        making_data.body_masses()
        making_data.body_velocities()
        making_data.finish()

        calc = atom_props.NonLocalLammps(snap, options.lammps_in, props_list)
        calc.run()

        for calc_name in props_list:
            if options.verbose is True:
                logger.user_message('Writing %s' % calc_name)

            snap.writer.hdf_create_array(calc_name,
                                    calc.result[calc_name],
                                    dtype=calc.dtypes[calc_name],
                                    location=calc.locations[calc_name],
                                    update_meta=True)

            if calc.types[calc_name] != 'scalars':
                snap.writer.update_composite_mapping(calc_name,
                                            calc.types[calc_name],
                                            dtype=calc.dtypes[calc_name],
                                            location=calc.locations[calc_name])

        snap.writer.xmf_snap_file()
        snap.writer.cleanup()
        os.remove(tmp_data_file_name)

    logger.completion_msg('Finished recover_atomwise.py: %s' % snap.sim_name)
    return 0

#-----------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
