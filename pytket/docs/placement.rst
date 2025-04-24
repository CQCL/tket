pytket.placement
==================================
In order for the constraints of a :py:class:`Backend` to be solved we must first assign device qubits to device-independent (or program) qubits. 
This module contains three placement methods to perform such an assignment. 

For more on qubit placement (and routing in general) see the `qubit mapping and routing <https://docs.quantinuum.com/tket/user-guide/examples/circuit_compilation/mapping_example.html>`_ tutorial and the corresponding entry in the `user manual <https://docs.quantinuum.com/tket/user-guide/manual/manual_compiler.html#placement>`_.

.. currentmodule:: pytket._tket.placement

.. autoclass:: pytket.placement.GraphPlacement

   .. automethod:: __init__
   .. automethod:: to_dict

.. autoclass:: pytket.placement.LinePlacement

   .. automethod:: __init__
   .. automethod:: to_dict

.. autoclass:: pytket.placement.NoiseAwarePlacement

   .. automethod:: __init__
   .. automethod:: to_dict

.. autoclass:: pytket.placement.Placement

   .. automethod:: __init__
   .. automethod:: from_dict
   .. automethod:: get_placement_map
   .. automethod:: get_placement_maps
   .. automethod:: place
   .. automethod:: place_with_map
   .. automethod:: to_dict

.. automethod:: pytket.placement.place_fully_connected

.. automethod:: pytket.placement.place_with_map
