#!/usr/local/bin/python
# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2

import os
import sys
import glob
from pkg_resources import require, DistributionNotFound

# Trying to use distribute but defaulting to setuptools
try:
    from distribute import setup, Extension, find_packages
    print("Using distribute....")
except ImportError:
    from setuptools import setup, Extension, find_packages
    print("Using setuptools....")

# ---------------------------------------------------------------------------
#               Checking Python
# ---------------------------------------------------------------------------

# Require Python>=2.6
if sys.version_info[:2] < (2, 6):
    print('mdhandle requires Python>=2.6.  Currently Python %d.%d' %
                        sys.version_info[:2])
    print('************************    Exiting ******************************')
    sys.exit(-1)

# ----------------------------------------------------------------------------
#               General Options
# ---------------------------------------------------------------------------

PACKAGE = "mdhandle"
NAME = "mdhandle"
DESCRIPTION = "Post-processing package for molecular dynamics simulation based on NumPy/SciPy for fluid mechanics, and material science simulation studies."
AUTHOR = "Dan Lussier, FBG Group, University of Oxford"
AUTHOR_EMAIL = "dtlussier@gmail.com"
URL = "https://github.com/dtlussier/mdhandle"
VERSION = __import__(PACKAGE).__version__
LICENSE = 'GPLv2'

# Currently removed in favour of using setup.py from LAMMPS itself.
# TODO: Add option for optionally compiling LAMMPS: 
    #http://www.sqlalchemy.org/trac/browser/setup.py
#http://stackoverflow.com/questions/2709278/setup-py-adding-options-aka-setup-py-enable-feature)
# ---------------------------------------------------------------------------
#               LAMMPS Extension
# ---------------------------------------------------------------------------

# TODO: See HACK below for building LAMMPS library.
#       - Problem with 'expected flat namespace' errors when building
#         in more conventional wy via setuptools.Extension

## list of src files for LAMMPS
#lmp_path = os.path.join(os.getcwd(), 'lammps-25Sep11-modified')
#
#libfiles = glob.glob( os.path.join(lmp_path, 'src', '*.cpp'))+\
#           glob.glob( os.path.join(lmp_path, 'src', 'STUBS', '*.cpp'))
#
## Building Serial LAMMPS library - see LAMMPS Python wrapper for parallel
#lammps_library = Extension("_lammps_serial",
#                           sources = libfiles,
#                           define_macros = [("MPICH_IGNORE_CXX_SEEK",1),
#                                            ("LAMMPS_GZIP",1),
#                                            ("FFT_NONE",1), ],
#                           # src files for LAMMPS and MPI STUBS
#                           include_dirs = [ os.path.join('..', 'src'),
#                                            os.path.join('..', 'src', 'STUBS')
#                                          ],
#                           extra_compile_args = []
#                           )
#
# ---------------------------------------------------------------------------
#               Running setup
# ---------------------------------------------------------------------------

script_files = glob.glob(os.path.join(os.getcwd(), 'bin', '*.py') )

# Looking for pytables - name in PyPi is tables and pytables in EPD.
try:
    require('pytables')
    pytables_name = 'pytables>=2.2' 
except DistributionNotFound:
    pytables_name = 'tables>=2.2'


setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=open("README.txt").read(),
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    url=URL,
    packages=find_packages(),

    # dictionary mapping package names to lists of wildcard patterns
    package_data={},
    scripts=script_files,
#    ext_modules = [lammps_library],
    # Ordered in revrse order of dependancies
    install_requires=[
            pytables_name,
            "numpydoc>=0.4",
            "Sphinx>=1.0",
            "scipy>=0.8",
            "numpy>=1.5",
            "iPython>=0.11"
            ],
    classifiers=[
        "Development Status :: 4 - Beta ",
        "Environment :: Console",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Researchers",
        "License ::  OSI Approved :: GPLv2",
        "Operating System :: POSIX",
        'Operating System :: MacOS :: MacOS X',
        "Programming Language :: Python",
                ],
    zip_safe=False,
    )

# ---------------------------------------------------------------------------
#               LAMMPS Extension - using LAMMPS setup.py
# ---------------------------------------------------------------------------
# HACK - building LAMMPS if doingo build or installing but not
if ('dist' not in sys.argv[1] and sys.argv[1] != 'register'):
    print('')
    print('------------------ Using LAMMPS setup.py script ------------------')
    print('')
    os.chdir('lammps-25Sep11-modified/python')
    os.system('python setup_serial.py %s' % sys.argv[1])
