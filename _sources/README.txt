.. mdhandle,  http://github.com/dtlussier/mdhandle
   Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
   Released under the GNU General Public License, v2

``mdhandle`` - README - v0.9
================================

Introduction
-------------

:mod:`mdhandle` is a post-processing framework for molecular dynamics 
simulation, particularly fluid mechanics, materials science and other
non-macromolecule simulation.

Currently, :mod:`mdhandle` is focused on interacting with result produced by the 
LAMMPS (http://lammps.sandia.gov) software package, but could be extended to 
handle the results from other software packages.

System Requirements
---------------------

1. **Linux**:  

* Requires ``Python 2.6`` or higher in the ``2.x`` series. Not compatible with ``Python 3.x`` due to dependancy with ``NumPy``.

* ``Python`` package dependencies:
	* ``NumPy``
	* ``SciPy``
	* ``matplotlib``
	* ``PyTables``
	* ``IPython``
	* ``sphinx`` and ``numpydoc`` for documentation.
	* Python dependancies are downloaded automatically during install
	  by ``easy_install`` or ``pip``.
        
2. **OS X**:  See Linux.

3. **Windows**: Windows support is unknown.  As shipped, :mod:`mdhandle` is developed for POSIX systems and may operate using the ``cygwin`` environment.
    

Install
---------

To install :mod:`mdhandle`, run ::

>>> python setup.py install

This command will install Python modules, and build the LAMMPS library
extension.  It is advisable to install the module rather than using it in place within the source tree.  To avoid polluting the system ``site-packages`` directory, ``virtualenv`` should be used.

Source and binary distributions can also be created by ::

>>> python setup.py sdist

or ::

>>> python setup.py bdist

Configuration
^^^^^^^^^^^^^^^

Systemwide settings are contained within :mod:`mdhandle.settings`.  They can
be modified within the file, or at runtime during program/script execution.


Scripts
--------

Example scripts can be found in the ``bin`` directory, which is installed along with the :mod:`mdhandle` module by the ``python setup.py install``.

Additional scripts can be found at
http://github.com/dtlussier/mdhandle_more_scripts

Documentation
----------------

After installation is complete, documentation can be found within the ``doc`` folder.

Sphinx (http://sphinx.pocoo.org/) is required to build the documentation 
via: ::

% cd doc
% make html

The resulting docs are then at ``doc/_build/html/index.html``.  Other available
build formats can be seen by invoking: ::

% cd doc
% make

Building the documentation at the local level reflects any changes made at the
local level.  :mod:`mdhandle` must be in your ``PYTHONPATH`` to build docs using ``sphinx``.  

A persistent version of the documentation reflecting the 
github version is at http://dtlussier.github.com/mdhandle/ and
on mirrored at http://mdhandle.readthedocs.org

Authors
---------

* Dan Lussier, Fluidics and Biocomplexity Group, Oxford University
  (dtlussier@gmail.com, http://github.com/dtlussier)


Acknowledgements
------------------

**Related Packages**:

* ``MMTK``, http://dirac.cnrs-orleans.fr/MMTK/ 
* ``MD-Tracks``, http://molmod.ugent.be/code/wiki/MD-Tracks
* ``mdanalysis``, http://code.google.com/p/mdanalysis/
* ``Pizza``, http://www.sandia.gov/~sjplimp/pizza.html

**Thanks**:

* Yiannis Ventikos, Fluidics and Biocomplexity Group, Oxford University
