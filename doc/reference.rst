.. General reference to all modules and classes within the :mod:`mdhandle`
   module.  New users should start with the Tutorial and the User Guide
   before consulting this reference.
   
Reference Documentation
========================

Modules within the :mod:`mdhandle` package can be sub-divided into
categories based on their functional role.  In the reference guide that follows
each module is listed within its functional group.

A similar organizational idea is also used within the source code.


Container Objects
-------------------

:mod:`~mdhandle.snap`
^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.snap

.. autoclass:: mdhandle.snap.Snap
   :members:

   .. document private functions
   .. automethod:: mdhandle.snap.Snap._use_selections
   

:mod:`~mdhandle.simcontainer`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.simcontainer
   :members:

:mod:`~mdhandle.sim_cell`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.sim_cell
   :members:

:mod:`~mdhandle.material`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.material
   :members:

Data Manipulators
-------------------

:mod:`~mdhandle.vectors`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.vectors
   :members:

:mod:`~mdhandle.tensors`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.tensors
   :members:


Property Calculations
------------------------

:mod:`~mdhandle.properties.selection`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.properties.selection
   :members: Selection

   Functions within :mod:`selection` are function factories that produce
   masking functions that are used to define a selection object.  They are
   passed as the ``masking_function`` argument when creating
   :class:`Selection`.
   
   Users can add additional new selection functions to :mod:`selection`,
   or within their own code which can be used as the ``masking_function``
   for a :class:`Selection` so long as they conform to the style seen in
   the native selection functions defined in :mod:`selection`.
   
.. automethod:: mdhandle.properties.selection.box
.. automethod:: mdhandle.properties.selection.sphere
.. automethod:: mdhandle.properties.selection.liquid_only
.. automethod:: mdhandle.properties.selection.fluid_only
.. automethod:: mdhandle.properties.selection.coordinate_less_than
.. automethod:: mdhandle.properties.selection.coordinate_greater_than
.. automethod:: mdhandle.properties.selection.opposite_selection

:mod:`~mdhandle.properties.snap_properties`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.properties.snap_properties
   :members:

:mod:`~mdhandle.properties.atom_properties`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.properties.atom_properties

.. autoclass:: mdhandle.properties.atom_properties.AtomCalc
   :members:

.. autoclass:: mdhandle.properties.atom_properties.NonLocal
   :members:
   :inherited-members:
   :show-inheritance:

   .. document private functions
   .. automethod:: _calc_epair
   .. automethod:: _calc_coord
   .. automethod:: _calc_force
   .. automethod:: _calc_virial 
   .. automethod:: _calc_stress 
   .. automethod:: _calc_diff

.. autoclass:: mdhandle.properties.atom_properties.NonLammps
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: mdhandle.properties.atom_properties.Local
   :members:
   :inherited-members:
   :show-inheritance:

   .. document private functions
   .. automethod:: _calc_ke


:mod:`~mdhandle.properties.binned_properties`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.properties.binned_properties

.. autoclass:: mdhandle.properties.binned_properties.Grid
   :members:

   .. document private functions
   .. automethod:: _requirements
   .. automethod:: _reverse_mapping
   .. automethod:: _sum_from_atom_value
   .. automethod:: _calc_count
   .. automethod:: _num_atoms
   .. automethod:: _coverage
   .. automethod:: _calc_num_density
   .. automethod:: _calc_mass
   .. automethod:: _calc_mass_density
   .. automethod:: _calc_epair_count
   .. automethod:: _calc_coord_count
   .. automethod:: _calc_velocity
   .. automethod:: _calc_ke_count
   .. automethod:: _calc_force_count
   .. automethod:: _calc_temp
   .. automethod:: _calc_temp_uncorr
   .. automethod:: _calc_stress

.. autoclass:: mdhandle.properties.binned_properties.RectGrid
   :members:
   :inherited-members:
   :show-inheritance:

   .. document private functions
   .. automethod:: _mapping_ngp
   .. automethod:: _mapping_cic
   .. automethod:: _requirements
   .. automethod:: _reverse_mapping
   .. automethod:: _sum_from_atom_value
   .. automethod:: _calc_count
   .. automethod:: _num_atoms
   .. automethod:: _coverage
   .. automethod:: _calc_num_density
   .. automethod:: _calc_mass
   .. automethod:: _calc_mass_density
   .. automethod:: _calc_epair_count
   .. automethod:: _calc_coord_count
   .. automethod:: _calc_velocity
   .. automethod:: _calc_ke_count
   .. automethod:: _calc_force_count
   .. automethod:: _calc_temp
   .. automethod:: _calc_temp_uncorr
   .. automethod:: _calc_stress

:mod:`~mdhandle.properties.cell_decomposition`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.properties.cell_decomposition

.. autoclass:: mdhandle.properties.cell_decomposition.SparseBinning
  :members:

  .. document private functions
  .. automethod:: _ngp

Statistics
------------

:mod:`~mdhandle.statistics.blockavg`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.statistics.blockavg
   :members:

:mod:`~mdhandle.statistics.fluctuations`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.statistics.fluctuations
   :members:

:mod:`~mdhandle.statistics.velocity_distribution`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.statistics.velocity_distribution
   :members:

Readers and Writers
---------------------

:mod:`~mdhandle.writers.xdmf`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.writers.xdmf
   :members:


:mod:`~mdhandle.writers.ascii`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.writers.ascii
   :members:


:mod:`~mdhandle.readers.column_mapping`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.readers.column_mapping
   :members:
   
   .. document private functions
   .. automethod:: mdhandle.readers.column_mapping._identify_collections

:mod:`~mdhandle.readers.lammps_dump.LAMMPS_dump`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.readers.lammps_dump
   :members:

:mod:`~mdhandle.writers.lammps_data`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.writers.lammps_data
   :members:


:mod:`~mdhandle.readers.log`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.readers.log
   :members:

:mod:`~mdhandle.readers.rdf`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.readers.rdf
   :members:

:mod:`~mdhandle.readers.eam`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.readers.eam
   :members:

LAMMPS Wrapping Module
--------------------------

.. automodule:: mdhandle.lammps

   .. autoclass:: mdhandle.lammps.lammps
      :members:

Helper Modules
-----------------

:mod:`~mdhandle.settings`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.settings
   :members:

:mod:`~mdhandle.units`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.units
   :members:


:mod:`~mdhandle.logger`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.logger
   :members:

:mod:`~mdhandle.utilities`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.utilities
   :members:
   
..  :mod:`~mdhandle.exceptions`
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    .. automodule:: mdhandle.exceptions
       :members:

:mod:`~mdhandle.interactive`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.interactive
   :members:

:mod:`~mdhandle.natsort`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: mdhandle.natsort
   :members:
