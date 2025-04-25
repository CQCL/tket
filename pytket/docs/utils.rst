pytket.utils
==================================
.. automodule:: pytket.utils
    :members: expectation_from_shots, expectation_from_counts, get_pauli_expectation_value, get_operator_expectation_value, append_pauli_measurement, counts_from_shot_table, probs_from_counts, probs_from_state, permute_qubits_in_statevector, permute_rows_cols_in_unitary, permute_basis_indexing, gen_term_sequence_circuit, readout_counts, compare_statevectors, compare_unitaries, prepare_circuit

.. automodule:: pytket.utils.outcomearray
.. autoclass:: pytket.utils.outcomearray.OutcomeArray
    :special-members: __init__
    :members:

.. autofunction:: pytket.utils.outcomearray.readout_counts

.. automodule:: pytket.utils.operators
.. autoclass:: pytket.utils.operators.QubitPauliOperator
    :special-members: __init__
    :members:

.. automodule:: pytket.utils.graph
.. autoclass:: pytket.utils.graph.Graph
    :special-members: __init__
    :members:

pytket.utils.distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: pytket.utils.distribution.EmpiricalDistribution
    :special-members: __eq__, __getitem__, __add__
    :members:

.. autoclass:: pytket.utils.distribution.ProbabilityDistribution
    :special-members: _eq__, __getitem__, __add__, __mul__
    :members:

.. automodule:: pytket.utils.distribution
    :members: convex_combination

pytket.utils.expectations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.utils.expectations
    :members:

pytket.utils.measurements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.utils.measurements
    :members:

pytket.utils.prepare
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.utils.prepare
    :members:

pytket.utils.results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.utils.results
    :members:

pytket.utils.serialization.migration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket._tket.utils_serialization
.. automodule:: pytket.utils.serialization
.. automodule:: pytket.utils.serialization.migration
    :members: circuit_dict_from_pytket1_dict

pytket.utils.spam
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.utils.spam
    :members: compress_counts

.. autoclass:: pytket.utils.spam.SpamCorrecter
    :special-members: __init__
    :members:

pytket.utils.stats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.utils.stats
    :members: gate_counts

pytket.utils.symbolic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.utils.symbolic
    :members:
    :special-members: SymGateFunc

pytket.utils.term_sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.utils.term_sequence
    :members:
