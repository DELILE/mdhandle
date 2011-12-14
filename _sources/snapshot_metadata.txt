.. mdhandle,  http://github.com/dtlussier/mdhandle
   Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
   Released under the GNU General Public License, v2

snap metadata
-------------

* ``col_mapping``
* ``vectors``
* ``tensors``
* ``symm_tensors``
* ``sim_cell``

Within the :attr:`Snap.meta` dictionary the following metadata is required :

* ``grid_type``
* ``xlo``
* ``xhi``
* ``ylo``
* ``yhi``
* ``zlo``
* ``zhi``
* ``num_atoms``
* ``col_mapping`` : 
    - List of all scalar values and their associated metadata.
* ``vectors`` :
    - List of all vector values and their associated metadata.
* ``tensors`` : 
    List of all tensor values and their associated metadata.
* ``symm_tensors`` : 
    List of all symmetric tensor values and their associated metadata.
* ``comments``
* ``time``
* ``units``
* ``sim_name``
* ``timestep``

Gridded
++++++++

For gridded datasets the following items are also added to the
:attr:`Snap.meta` :

* ``origin``
* ``num_elem``
* ``boxL``
* ``grid_spacing``
