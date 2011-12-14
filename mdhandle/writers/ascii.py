# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Takes a :class:`mdhandle.snap.Snap` and outputs atom data to a LAMMPS
formatted dump file.  Columns in the LAMMPS dump file are not ordered according
to original LAMMPS output file beyond atom ``id``, which remains the first
column.

"""

# TODO: Ability to handle simulations that are not just atoms
#       (e.g. polyatomic molecules)
# TODO: Ability to output only a subset of columns, while ensuring that
#       important columns are included (e.g. id, type, x, y, z)

# TODO: Add ability to handle case where only getting subset of atoms.
#       - Currently hard coded self.snap.meta['num_atoms']
#       - Extraction of atomwise data already able to handle use of selection.

# ------------------------------------------------------------------------------

import sys
import os

import numpy as np

from mdhandle.logger import Logger

# ------------------------------------------------------------------------------

class AsciiWriter(object):
    """
    Class for writing ASCII dump files for a simulation snapshot in the
    format of LAMMPS ASCII output.

    Parameters
    ----------
    snap : :class:`mdhandle.snap.Snap`
        Current simulation snapshot.
        Active dataset should be of type ``'atoms'``.
    verbose : boolean,  [default=False]
        Flag indicating verbosity of user input and output.

    Attributes
    -----------
    snap : :class:`mdhandle.snap.Snap`
    verbose : boolean
        If ``True`` user I/O is verbose
    logger : :class:`mdhandle.logger.Logger` 
        Object for user I/O and logging.
    ascii_fn : string
        Name of output ASCII LAMMPS dump file.
    output_dir : string
        Absolute path to output directory for LAMMPS dump file.
    column_names : list of string
        Names of atomwise scalars (i.e. ``'x'``, ``'fx'``, ``'epair'``)
    header : string
        LAMMPS header - See LAMMPS user guide for more information.
    _f : file
        File handle to open ``ascii_fn``.

    Methods
    -------
    new: Creates ASCII file and writes LAMMPS dump header.
    write: Writes body of LAMMPS dump (i.e. atomwise data) to existing file.
    finish: Flushes buffers and closes ASCII file.

    References
    -----------
    * LAMMPS ``dump`` comand, http://lammps.sandia.gov/doc/dump.html

    """

    def __init__(self, snap, verbose=False):
        self.snap = snap
        self.verbose = verbose

        self.logger = Logger(self.verbose)

        self.snap.gather_data()
        
        if self.snap.meta['grid_type'] != 'atoms':
            self.logger.error("Active dataset is of type %s not 'atoms'." %
                                                   self.snap.meta['grid_type'])

        self.ascii_fn = 'ASCII_%s+.dat' % os.path.splitext(self.snap.fn_base)[0]
        self._f = None
        self.output_dir = self.snap.directory

        # Order is preserved by getting it in list form
        self.column_names = self.snap.col_mapping.keys()

        # Col for atom id is always the first column
        self.column_names.remove('id')
        # Col for atom type is always the second column
        self.column_names.remove('type')

        # TODO: meld with header constant in mdhandle.settings
        self.header = '''ITEM: TIMESTEP
        %d
        ITEM: NUMBER OF ATOMS
        %d
        ITEM: BOX BOUNDS
        %f %f
        %f %f
        %f %f
        ITEM: ATOMS id type %s
        ''' % (self.snap.meta['time'], self.snap.meta['num_atoms'],
               self.snap.meta['xlo'], self.snap.meta['xhi'],
               self.snap.meta['ylo'], self.snap.meta['yhi'],
               self.snap.meta['zlo'], self.snap.meta['zhi'],
               ' '.join(self.column_names) )

    def new(self):
        """
        Creates new ASCII file and writes necessary LAMMPS header.

        Header is new-style (post-2009) with column names in the final
        ``'ITEM:ATOMS ... '`` line.

        """
        # np.savetxt(..) will have optional header parameter for numpy>=2.0

        self._f = open(self.ascii_fn, 'w')
        self._f.write(self.header)

    def write(self):
        """
        Write atomwise data to ASCII dump file..

        """
        # Checking that new file has been opened.
        if self._f is None:
            self.logger.warning('Ouput file %s does not exist. Use new().' %
                                                self.ascii_fn)
            return
        elif self._f.closed is False:
            self.logger.warning('Output file %s is closed.' % self.ascii_fn)
            return

        # Empty array
        data = np.zeros((self.snap.meta['num_atoms'],
                         len(self.snap.col_mapping) ))

        data[:,0] = self.snap.get_scalar('id')
        data[:,1] = self.snap.get_scalar('type')

        line_format = ''
        for i,col_name in enumerate(self.column_names):
            data[:,2+i] = self.snap.get_scalar(col_name)
            if self.snap.col_mapping[col_name]['dtype'].lower() == 'int':
                line_format += '%d '
            elif self.snap.col_mapping[col_name]['dtype'].lower() == 'float':
                line_format += '%f '
            else:
                raise UserWarning('Column %s dtype is not int or float.')

        np.savetxt(self._f, data, fmt=line_format, delimiter=' ')

    def finish(self):
        """
        Closes ASCII file and performs any other finalizing operations on new
        ASCII file.

        """
        self._f.close()

        if self.verbose is True:
            self.logger.completion_msg('ASCII file creaed %s' % self.ascii_fn)


#=============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
