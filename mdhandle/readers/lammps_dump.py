#!/usr/bin/env python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Reads LAMMPS ASCII dump files and extracts metadata.

"""
import sys
import os

import numpy as np
import tables as pyT

from mdhandle import settings
from mdhandle.logger import Logger

# -----------------------------------------------------------------------------

logger = Logger()

# -----------------------------------------------------------------------------


# TODO: Refactor to remove use of os.system, sed, temporary files 
#       (at least the HDF file)
class LAMMPS_dump(object):
    """
    Reads LAMMPS ASCII dump files and extracts metadata.

    Parameters
    -----------
    fn : string
        Filename of LAMMPS dump file.

    Attributes
    -----------
    fn : string
        LAMMPS ASCII dump file.
    atom_data : :class:`numpy.ndarray` (N x ncol)
        2-dimensional array of atomwise data.
    meta : dict
        Dictionary of LAMMPS metadata, read from file header.
    is_self_describing : boolean
        If ``True``, LAMMPS dump is self-describing.
    _inc : int
        Number of rows of dump file, which are processed in
        a single step.

    """

    def __init__(self, fn):
        assert os.path.exists(fn) is True

        self.fn = fn
        self._inc = 5000000
        self.meta = {}
        self.is_self_describing = False
        self.atom_data = None

    def read_header(self):
        """
        Reads the LAMMPS dump header.

        """
        f = open(self.fn)

        # TODO: use LAMMPS_HEADER info from mdhandle.settings to parse
        # readline() statements configured for LAMMPS ASCII format.
        f.readline()
        self.meta['time'] = int(f.readline())
        f.readline()
        self.meta['num_atoms'] = int(f.readline())

        f.readline()
        words = f.readline().split()
        self.meta['xlo'], self.meta['xhi'] = float(words[0]), float(words[1])
        words = f.readline().split()
        self.meta['ylo'], self.meta['yhi'] = float(words[0]), float(words[1])
        words = f.readline().split()
        self.meta['zlo'], self.meta['zhi'] = float(words[0]), float(words[1])

        # First two words of line are 'ITEM:' 'ATOMS'.
        # If old-style LAMMPS dump, self.meta['description'] == []
        self.meta['description'] = f.readline().split[2:]

        if len(self.meta['description']) > 0:
            self.is_self_describing = True
        else:
            self.is_self_describing = False

        # Getting actual number of columns from first line in file.
        self.meta['ncol'] = len(f.readline().split())

        f.close()

    def read_body(self, col_mapping=None):
        """
        Reads body of LAMMPS dump file containing atomwise data.

        Parameters
        -----------
        col_mapping : dict, [default=None]
            Mapping of column name to column number.  Required for handling
            pre-2007 LAMMPS output, which is not self-describing.
            Leave as ``None`` if have self-describing LAMMPS dump.

        """
        # TODO: Refactor to clarify the execution flow.

        tmp_list = []
        range_lo = 1 + len(settings.LAMMPSHEADER)
        range_hi = self._inc + len(settings.LAMMPSHEADER)

        fn_hdf = os.path.join(settings.TMP_LOCATION,
                              os.path.basename(self.fn) + '_EArray_tmp.h5')

        # Splitting ASCII dump into smaller chunks for processing w sed
        for i in range(0, int(round(self.meta['num_atoms']/self._inc) + 1)):

            fn_tmp = os.path.join(settings.TMP_LOCATION,
                                  os.path.basename(self.fn)+'_'+str(i)+".tmp")
            tmp_list.append(fn_tmp)
            os.system("sed -n '" + str(range_hi+1) + "q;" + str(range_lo) + ","
                            + str(range_hi) + "p' " + self.fn + " > " + fn_tmp)
            range_lo = range_hi + 1
            range_hi = range_lo + self._inc - 1

        if self.meta['num_atoms'] > self._inc:
            f_hdf = pyT.openFile(fn_hdf, 'w')
            # 0 arg indicates that EArray can add rows later.
            f_hdf.createEArray(f_hdf.root, 'atom_data', pyT.FloatAtom(),
                               (0, self.meta['ncol']),
                               expectedrows=self.meta['num_atoms'])
            f_hdf.close()

            for fn_tmp in tmp_list:
                _append_data_hdf(fn_tmp, fn_hdf, self.meta['ncol'])

            f_hdf = pyT.openFile(fn_hdf, 'r')
            self.atom_data = f_hdf.root.atom_data.read()
            f_hdf.close()
            os.remove(fn_hdf)

        # self.meta['num_atoms'] <= increment
        else:
            # TODO: Use col_mapping to get dtypes right from the start
            # TODO: Use names to get array with fields.
            raw = np.fromfile(tmp_list[0], dtype=np.float, sep=' ')
            self.atom_data = raw.reshape( (raw.shape[0]/self.meta['ncol'],
                                            self.meta['ncol']))

        for fn_tmp in tmp_list:
            os.remove(fn_tmp)

        # Sorting rows - may be out of order due to parallel I/O.

        if self.is_self_describing is True:
            if 'id' not in self.meta['description']:
                logger.error("LAMMPS dump must have 'id' column to sort")
            id_idx = self.meta['description'].index('id')
        elif col_mapping is not None:
            id_idx = col_mapping['id']['col']
        else:
            logger.error('LAMMPS dump is not self-describing.\
                          col_mapping is required to sort rows.')

        # TODO: Alternative exists for using np.sort() with labeled arrays.
        #       see SO for example.
        self.atom_data = self.atom_data[ self.atom_data[:, id_idx].argsort() ]

    def read_column(self, name, col_mapping=None):
        """
        Read a single column from LAMMPS ASCII dump.

        Raises
        -------
        NotImplementedError : :class:`Exception`

        """
        logger.error('LAMMPS_dump.read_column(), Not Implemented.')


def _append_data_hdf(fn, fn_hdf, ncol):
    """
    Helper function used as target by :mod:`multiprocessing.Process` object.

    Reads ASCII file with only table of floating point numbers
    (space separated) and appends it to atomData array in h5TmpFile.

    Parameters
    ----------
    fn : string
        Name of ASCII file containing large table of floats.
    fn_hdf : string
        HDF file with atom_data EArray (extendible) hanging from root.
    ncol : int
        Number of columns in LAMMPS atomwise data.

    """
    raw = np.fromfile(fn, dtype=np.float, sep=' ')
    data = np.reshape(raw, (raw.shape[0]/ncol, ncol))
    f_hdf = pyT.openFile(fn_hdf, 'a')
    f_hdf.root.atom_data.append(data)
    f_hdf.close()


# -----------------------------------------------------------------------------


def main():
    pass

if __name__ == "__main__":
    sys.exit(main())
