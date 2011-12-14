# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Parses LAMMPS RDF files to provide analysis and graphing functionality.

"""

# TODO: Add ability to handle multiple atom types.
# TODO: Add ability to select timesteps from within RDF data file.
# TODO: Add compute for standard deviation
# TODO: Get nbins from rdf data file rather than as input to RDF __init__(...)
# TODO: Use existing LAMMPS averages which exist after any run command ends
# TODO: Add subset argument to plot_rdf_movie to allow for use of only some of
#       total data.

# -----------------------------------------------------------------------------

import sys
from StringIO import StringIO

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from mdhandle.logger import Logger

# -----------------------------------------------------------------------------

logger = Logger()

# -----------------------------------------------------------------------------


class RDF(object):
    """
    Class for parsing and processing LAMMPS RDF data produced by
    ``compute rdf`` command.

    LAMMPS RDF data:
    
    * Always ASCII.
    * File broken down into blocks based on timestep.
    * At each timestep RDF data given as a function of interatomic
      distance (``r``)
    * At each timestep three columns : ``r rdf coord``
    * ``coord`` is the sum of rdf up to that point.

    Parameters
    ----------
    fn : string
        Name of rdf file.
    nbins : int
        Number of bins used in LAMMPS compute command to produce data.
    comment_char : string,  [default='#']
        Comment character - skip lines with leading comment character.
    verbose : boolean,  [default=False]
        If ``True``, verbose user I/O.

    Attributes
    ----------
    fn : string
        RDF filename.
    verbose : boolean
        If ``True``, verbose user I/O.
    nbins : int
        Number of bins used in calculation of RDF file.
    comment : string
        Comment character in RDF file.
    times : list of int
        List of timestep values at which RDF function was stored.
    is_lammps_averaged : boolean
        If ``True`` average RDF has been computed by LAMMPS.
    computed_avg : :class:`numpy.ndarray`
        Average RDF computed from RDF file read from disk.
    data : :class:`numpy.ndarray`
        3-dimensional array storing tables of (``r rdf coord``).
    r : :class:`numpy.ndarray`.
        Discrete sample of radius used as indep variable in RDF function.
        Cutoff used in LAMMPS calculation can be inverred from ``r.max()``

    Methods
    --------
    average()
        Calculate average of RDF functions over time.
    plot_rdf(step_num, clear=True, smoothing=False)
        Plot single rdf file.
    plot_rdf_movie(clear=True, save=True, fmt='.png')
        Plot *movie* of RDF function over time.

    """

    def __init__(self, fn, nbins, comment_char='#', verbose=False):
        self.fn = fn
        self.nbins = nbins
        self.verbose = verbose

        f_rdf = open(self.fn, 'r')

        self.times = []
        self.is_lammps_averaged = False
        self.computed_avg = None

        num_steps = 0
        for line in f_rdf:
            if ('%s Timestep' % comment_char) in line:
                num_steps += 1
                self.times.append(line.split()[-2])
            elif ('%s Run Average' % comment_char) in line:
                self.is_lammps_averaged = True

        f_rdf.seek(0)

        # Filter off Header files so just data is left.
        good_data = ''
        f_iter = iter(f_rdf)
        for line in f_iter:
            if '# Timestep' in line:
                # Also skip the next line: r, g(1,1,r), ncoord(1,1,r)
                f_iter.next()
            elif '# Run Average' in line:
                # Also skip the next line: r, g(1,1,r), ncoord(1,1,r)
                f_iter.next()
                # Skipping the averaged RDF results under "# Run Average"
                for j in range(nbins):
                    f_iter.next()
            else:
                good_data += line

        # Making numpy.ndarray from good lines of RDF data file
        raw_data = np.loadtxt(StringIO(good_data), dtype=np.float)

        # HARDCODED: number of columns.
        num_cols = 3
        self.data = raw_data.reshape(num_steps, self.nbins, num_cols)
        self.r = self.data[0,:,0]

        f_iter.close()

    def average(self):
        """
        Computing average RDF in time.

        """
        self.computed_avg = self.data.mean(axis=0)

    def plot_rdf(self, step_num, clear=True, smoothing=False):
        """
        Plots RDF as a function of the interparticle distance (``r``).

        Parameters
        ----------
        step_num : int
            Index within :attr:`~mdhandle.rdf.RDF.data` to select
            one set of RDF data.
        clear : boolean, [default=True]
            If ``True``, creates new figure for RDF plot.
            If ``False``, RDF plot is created on existing matplotlib window
            if it exists.
        smoothing : boolean,  [default=False]
            If ``True``, use :meth:`scipy.interp1d` to do cubic interpolation
            to RDF data to produce a smoother line.
            If ``False``, plots raw RDF data.

        Returns
        --------
            Pointer to :meth:`matplotlib` figure containing RDF data.

        """
        rdf_to_plot = self.data[step_num,:,:]
        logger.user_message('Plotting RDF from timestep %d' %
                                                        self.times[step_num])

        if clear is True:
            plt.figure()
        if smoothing is False:
            return plt.plot(self.r, rdf_to_plot[:,1])
        elif smoothing is True:
            interp_res = interp1d(self.r, rdf_to_plot[:,1], kind='cubic')
            xnew = np.linspace(self.r.min(), self.r.max(), 10*self.nbins)
            return plt.plot(xnew, interp_res(xnew))

    def plot_rdf_movie(self, clear=True, save=True, fmt='.png'):
        """
        Plot sequence of RDF output for each timestep in the RDF data file.

        Intention is for diagnostic use only as no plotting controls and
        file name control.

        Parameters
        -----------
        clear : boolean, [default=True]
            If ``True``, creates new figure for RDF plot.
            If ``False``, RDF plot is created on existing matplotlib window
            if it exists.
        save : boolean,  [default=True]
            If ``True`` - plotting output is saved to disk.
        fmt : string,  [default='.png']
            String defining output file format for matplotlib.

        """
        plt.clf()
        for i in range(self.data.shape[0]):
            plt.plot(self.r, self.data[i,:,1])
            if save is True:
                plt.savefig(str(i)+fmt)
            else:
                plt.draw()
            if clear is True:
                plt.clf()


# ============================================================================


def main(argv=None):
    pass


if __name__ == "__main__":
    sys.exit(main())
