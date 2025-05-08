pytket.pauli
==================================
.. currentmodule:: pytket._tket.pauli

.. automodule:: pytket.pauli
.. automodule:: pytket._tket.pauli

.. autoenum:: pytket.pauli.Pauli

.. autoclass:: pytket.pauli.PauliStabiliser

   .. automethod:: __init__
   .. autoproperty:: coeff
   .. autoproperty:: string

.. autoclass:: pytket.pauli.QubitPauliString

   .. automethod:: __init__
   .. automethod:: commutes_with
   .. automethod:: compress
   .. automethod:: dot_state
   .. automethod:: from_list
   .. autoproperty:: map
   .. automethod:: state_expectation
   .. automethod:: to_list
   .. automethod:: to_sparse_matrix

.. autoclass:: pytket.pauli.QubitPauliTensor

   .. automethod:: __init__
   .. autoproperty:: coeff
   .. automethod:: commutes_with
   .. automethod:: compress
   .. automethod:: dot_state
   .. automethod:: state_expectation
   .. autoproperty:: string
   .. automethod:: to_sparse_matrix
