pytket.partition
==================================
.. currentmodule:: pytket._tket.partition

.. autoenum:: pytket.partition.GraphColourMethod

.. autoclass:: pytket.partition.MeasurementBitMap

   .. automethod:: __init__
   .. automethod:: from_dict
   .. automethod:: to_dict
   .. autoproperty:: bits
   .. autoproperty:: circ_index
   .. autoproperty:: invert

.. autoclass:: pytket.partition.MeasurementSetup

   .. automethod:: __init__
   .. automethod:: add_measurement_circuit
   .. automethod:: add_result_for_term
   .. automethod:: from_dict
   .. automethod:: to_dict
   .. automethod:: verify
   .. autoproperty:: measurement_circs
   .. autoproperty:: results

.. autoenum:: pytket.partition.PauliPartitionStrat

.. automethod:: pytket.partition.measurement_reduction

.. automethod:: pytket.partition.term_sequence
