# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Reads LAMMPS log files.

Provides wrapper to ``log.py`` in Pizza post-processing package.

**Credits** :

* Uses ``log.py`` from Pizza v19Apr10 for reading LAMMPS log files
  (http://www.sandia.gov/~sjplimp/pizza.html).  Within :mod:`mdhandle`
  it is renamed ``log_pizza.py``.

"""
# TODO: Embed the log data within one or all of HDF files.
# TODO: Add ability to use a list of log files rather than string to search
# TODO: Store log data w/i Log object as attribute (setattr)
# TODO: Write log data to HDF file so that it is more useful at later time.

import sys

import numpy as np

# LAMMPS log class.
from mdhandle.readers.log_pizza import log
from mdhandle.logger import Logger

#------------------------------------------------------------------------------

logger = Logger()

#------------------------------------------------------------------------------


class Log(object):
    """

    Attributes
    -----------
    log_fn : string
        Search string for LAMMPS log file.  May include wildcards.
    flist : list
        List of LAMMPS log files captured by ``log_fn``
    num : int
        Number of columns in LAMMPS log file.
    len : int
        Length of columns in LAMMPS log file.
    names : list
        Name of columns in LAMMPS log file.  Act as input to :meth:`get_vec`.
    _log_pizza : file
        Handle to Pizza log file object.

    Methods
    --------
    get_vec(name, dtype='float')
        Returns a vector for single log file column (i.e. kinetic energy).
    get_all_data()
        Returns array with complete log file data.

    See Also
    ---------
    :mod:`mdhandle.readers.log`

    References
    -----------
    * Pizza ``log.py``, http://www.sandia.gov/~sjplimp/pizza/doc/log.html

    """

    def __init__(self, log_fn):
        """
        Parameters
        -----------
        log_fn : string
            Search string for LAMMPS log files.
            Shell wildcards may be used to pull in many affiliated log
            files from the same simulation.

        """
        self.log_fn = log_fn
        self._log_pizza = log.log(self.log_fn)

        # Pulling metadata up into Log(...) object
        self.flist = self._log_pizza.flist
        self.num = self._log_pizza.nvec
        self.len = self._log_pizza.nlen
        self.names = self._log_pizza.names

    def get_vec(self, name, dtype='float'):
        """
        Retrieves data column from log file.

        Parameters
        ----------
        name : string
            Name of column from log file.
        dtype : string
            Desired data type of log data.
            {'int' | 'float'}
            [default='float']

        Returns
        -------
        vec : :class:`numpy.ndarray`
            Vector for column of log data.

        """
        vec = np.array(self._log_pizza.get(name), dtype=dtype)
        return vec

    def get_all_data(self):
        """
        Retrieves all data from log file as a large two-dimensional
        :class:`numpy.ndarray`.

        Returns
        --------
        block : :class:`numpy.ndarray`
            All log data.

        """
        block = np.transpose(np.vstack([self.get_vec(i) for i in self.names]))
        return block


# ============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
