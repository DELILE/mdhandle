# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Writes out atom information for ASCII LAMMPS data file, used by LAMMPS
``read_data`` command.  Useful for restarting when do not have restart or
for starting more complex initial molecular geometries.

"""

#TODO: Handle other types of LAMMPS data (i.e. molecular information like
#      angle, dihedral, etc.)
#TODO: Add controls to manage precision of output, maybe be in the global
#      mdhandle package settings
#TODO: Add functionality so could write a subset of a complete snap?
# TODO: Add check to ensure that newly created LAMMPS data file has required
#       parts (i.e. header, etc.)

# -----------------------------------------------------------------------------

import os
import sys

import numpy as np

from mdhandle.logger import Logger

# -----------------------------------------------------------------------------

class LAMMPS_Data(object):
    """
    Writes LAMMPS data file for snap object.

    Parameters
    ----------
    snap : :class:`mdhandle.snap.Snap`
        Current simulation snapshot.
    verbose : boolean,  [default=False]
        Flag for verbosity of user input and output.

    Attributes
    ----------
    snap : :class:`mdhandle.snap.Snap`
    verbose : boolean
        If ``True``, verbose user I/O.
    logger : :class:`mdhandle.logger.Logger`
        Object for user I/O and logging.
    atom_types : :class:`numpy.ndarray`, (Nx1) 
        Atomwise vector of atom types.
    atom_ids : :class:`numpy.ndarray`, (Nx1) 
        Atomwise vector of unique atom ids.
    data_fn : string
        Name of LAMMPS data file.
    _f : file
        File handle for LAMMPS data file.

    Methods
    --------
    new(data_fn=None)
        Create new ASCII LAMMPS data file.
    header(non_pbc_pad=(0., 0., 0.))
        Write LAMMPS data file header. Manditory section.
    body_atoms(wrap=False)
        Write section for position of atoms. Manditory section.
    body_velocities()
        Write section for velocities of atoms.
    body_masses()
        Write section for atomic masses.
    finish()
        Flush buffers and close ASCII LAMMPS data file.

    References
    ------------
    * LAMMPS ``read_data`` command - http://lammps.sandia.gov/doc/read_data.html

    """

    def __init__(self, snap, verbose=False):
        self.snap = snap
        self.verbose = verbose

        self.logger = Logger(self.verbose)

        self.atom_types = self.snap.get_scalar('type')
        self.atom_ids = self.snap.get_scalar('id')
        
        self.data_fn = ''
        self._f = None

    def new(self, data_fn=None):
        """
        Creates new file for LAMMPS data file.

        Parameters
        ----------
        data_fn : string, [default=None]
            Optional parameter to provide new name for data file.
            If not given, data file name is ``Snap`` file name with
            ``'.data'`` extension.

        """
        if data_fn is None:
            self.data_fn = os.path.join(self.snap.directory,
                            os.path.splitext(self.snap.fn_base)[0] + '.data')
        else:
            self.data_fn = data_fn

        self._f = open(self.data_fn, 'w')

    # TODO: Add additional optional header sections for bonds, angles,
    # dihedrals, impropers, bond types, angle types, dihedral types,
    # improper types, extra bond per atom (leaves space in memory for new
    # bonds).

    def header(self, non_pbc_pad=(0., 0., 0.)):
        """
        Writes header for LAMMPS data file (> December 2009).  Structure
        of data file is discussed in LAMMPS ``read_data`` command
        documentation.

        .. note:: The header is a manditory section.

        Parameters
        ----------
        pbc_pad : iterable, [default=(0., 0., 0.)]
            Vector (``len(pbc_pad) == 3``) giving size of padding added
            around edges of simulation cell to avoid undesirable overlaps.

        """
        non_pbc_pad = np.array(non_pbc_pad)
        # ASSERT: Padding for non-pbc >= 0.
        assert np.all(non_pbc_pad >= 0.)

        self._f.write('\n')
        self._f.write('# Data file for: %s\n' % self.snap.filename)
        self._f.write('\n')
        self._f.write('%d atoms \n' % self.snap.meta['num_atoms'])

        self._f.write('%d atom types\n' % np.unique(self.atom_types).size)

        self._f.write('%f %f xlo xhi\n' % (self.snap.meta['xlo'] -
                                                                non_pbc_pad[0],
                                          self.snap.meta['xhi'] +
                                                               non_pbc_pad[0]))
        self._f.write('%f %f ylo yhi\n' % (self.snap.meta['ylo'] -
                                                                non_pbc_pad[1],
                                          self.snap.meta['yhi'] +
                                                               non_pbc_pad[1]))
        self._f.write('%f %f zlo zhi\n' % (self.snap.meta['zlo'] -
                                                                non_pbc_pad[2],
                                          self.snap.meta['zhi'] +
                                                               non_pbc_pad[2]))

        self._f.write('\n')

    def body_atoms(self, wrap=False):
        """
        Write ``Atoms`` section describing the position of each atom.

        .. note:: Manditory section  if number of atoms  > 0

        Parameters
        ----------
        wrap:   boolean,  [default=False]
            If ``True``, atom coordinates are wrapped back into central
            simulaton cell from their current periodic image.

        """
        self._f.write('\n')
        self._f.write('   Atoms\n')
        self._f.write('\n')

        xyz = self.snap.get_vector(('x', 'y', 'z'))

        if wrap is True:
            if self.verbose is True:
                self.logger.user_message('Wrapping coordinates back to central\
                                                                simulation box')
            xyz = self.snap.sim_cell.wrap_to_central(xyz)

        np.savetxt(self._f, np.column_stack(
                                        (self.atom_ids, self.atom_types, xyz)),
                                        fmt='%d %d %f %f %f')
        self._f.write('\n')

    def body_velocities(self):
        """
        Writes ``Velocities`` section for atomwise velocities.
        Can alternatively be specified via the LAMMPS ``velocity`` command.

        """

        self._f.write('\n')
        self._f.write('   Velocities\n')
        self._f.write('\n')

        # TODO: This format does not apply to atom_type dipole, ellipsoid,
        #       granular
        velocity = self.snap.get_vector(self.snap.vectors['velocity']['comps'])
        np.savetxt(self._f, np.column_stack((self.atom_ids, velocity)),
                                   fmt='%d %f %f %f')
        self._f.write('\n')

    def body_masses(self):
        """
        Writes ``Masses`` section for atomwise masses.
        Can also be done via LAMMPS ``mass`` commands

        """

        self._f.write('\n')
        self._f.write('   Masses\n')
        self._f.write('\n')

        for i in range(np.unique(self.atom_types).size):
            self._f.write('%d %f\n' % (np.unique(self.atom_types)[i],
                                    np.unique(self.snap.get_scalar('mass'))[i]))
        self._f.write('\n')

    # TODO: body_potential and body_bonds functions for writing parameters for
    #       atomwise interatomic potential and bonding.

    def finish(self):
        """
        Closes data file and performs any finalizing operations.

        """
        self._f.close()

        if self.verbose is True:
            self.logger.completion_msg('LAMMPS data file creaed %s' %
                                                                self.data_fn)


# =============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
