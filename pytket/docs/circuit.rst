pytket.circuit
==================================
.. toctree::
   :caption: Classes:
   :maxdepth: 1

   circuit_class.rst
   optype.rst
   classical.rst
   display.md

.. currentmodule:: pytket._tket.circuit

.. automodule:: pytket.circuit
.. automodule:: pytket._tket.circuit

.. automethod:: pytket.circuit.fresh_symbol

.. autoclass:: pytket.circuit.Op

   .. automethod:: create
   .. automethod:: free_symbols
   .. automethod:: get_name
   .. automethod:: get_unitary
   .. automethod:: is_clifford
   .. automethod:: is_clifford_type
   .. automethod:: is_gate
   .. autoproperty:: dagger
   .. autoproperty:: n_qubits
   .. autoproperty:: params
   .. autoproperty:: transpose
   .. autoproperty:: type

.. autoclass:: pytket.circuit.Command

   .. automethod:: free_symbols
   .. autoproperty:: args
   .. autoproperty:: bits
   .. autoproperty:: op
   .. autoproperty:: opgroup
   .. autoproperty:: qubits

.. autoenum:: pytket.circuit.BasisOrder

.. autoenum:: pytket.circuit.CXConfigType

.. autoenum:: pytket.circuit.EdgeType

.. autoclass:: pytket.circuit.CircBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: symbol_substitution
   .. autoproperty:: circuit_name

.. autoclass:: pytket.circuit.Unitary1qBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_matrix

.. autoclass:: pytket.circuit.Unitary2qBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_matrix

.. autoclass:: pytket.circuit.Unitary3qBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_matrix

.. autoclass:: pytket.circuit.ExpBox

   .. automethod:: __init__
   .. automethod:: get_circuit

.. autoclass:: pytket.circuit.PauliExpBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_cx_config
   .. automethod:: get_paulis
   .. automethod:: get_phase

.. autoclass:: pytket.circuit.PauliExpPairBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_cx_config
   .. automethod:: get_paulis_pair
   .. automethod:: get_phase_pair

.. autoclass:: pytket.circuit.PauliExpCommutingSetBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_cx_config
   .. automethod:: get_paulis

.. autoclass:: pytket.circuit.TermSequenceBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_cx_config
   .. automethod:: get_depth_weight
   .. automethod:: get_graph_colouring_method
   .. automethod:: get_partition_strategy
   .. automethod:: get_paulis
   .. automethod:: get_synthesis_strategy

.. autoenum:: pytket.circuit.ToffoliBoxSynthStrat

.. autoclass:: pytket.circuit.ToffoliBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_permutation
   .. automethod:: get_rotation_axis
   .. automethod:: get_strat

.. autoclass:: pytket.circuit.QControlBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_control_state
   .. automethod:: get_control_state_bits
   .. automethod:: get_n_controls
   .. automethod:: get_op

.. autoclass:: pytket.circuit.CustomGateDef

   .. automethod:: define
   .. automethod:: from_dict
   .. automethod:: to_dict
   .. autoproperty:: args
   .. autoproperty:: arity
   .. autoproperty:: definition
   .. autoproperty:: name

.. autoclass:: pytket.circuit.CustomGate

   .. automethod:: get_circuit
   .. autoproperty:: gate
   .. autoproperty:: name
   .. autoproperty:: params

.. autoclass:: pytket.circuit.Conditional

   .. autoproperty:: op
   .. autoproperty:: value
   .. autoproperty:: width

.. autoclass:: pytket.circuit.ClExprOp

   .. automethod:: __init__
   .. autoproperty:: expr
   .. autoproperty:: type

.. autoclass:: pytket.circuit.WiredClExpr

   .. automethod:: __init__
   .. automethod:: from_dict
   .. automethod:: to_dict
   .. autoproperty:: bit_posn
   .. autoproperty:: expr
   .. autoproperty:: output_posn
   .. autoproperty:: reg_posn

.. autoclass:: pytket.circuit.ClExpr

   .. automethod:: __init__
   .. automethod:: as_qasm
   .. autoproperty:: args
   .. autoproperty:: op

.. autoenum:: pytket.circuit.ClOp

.. autoclass:: pytket.circuit.ClBitVar

   .. automethod:: __init__
   .. autoproperty:: index

.. autoclass:: pytket.circuit.ClRegVar

   .. automethod:: __init__
   .. autoproperty:: index

.. autoclass:: pytket.circuit.PhasePolyBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. autoproperty:: linear_transformation
   .. autoproperty:: n_qubits
   .. autoproperty:: phase_polynomial
   .. autoproperty:: phase_polynomial_as_list
   .. autoproperty:: qubit_indices

.. autoclass:: pytket.circuit.ProjectorAssertionBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_matrix

.. autoclass:: pytket.circuit.StabiliserAssertionBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_stabilisers

.. autoclass:: pytket.circuit.WASMOp

   .. automethod:: __init__
   .. autoproperty:: func_name
   .. autoproperty:: input_widths
   .. autoproperty:: n_i32
   .. autoproperty:: num_bits
   .. autoproperty:: num_w
   .. autoproperty:: output_widths
   .. autoproperty:: wasm_uid

.. autoclass:: pytket.circuit.MultiBitOp

   .. automethod:: __init__
   .. autoproperty:: basic_op
   .. autoproperty:: multiplier

.. autoclass:: pytket.circuit.SetBitsOp

   .. automethod:: __init__
   .. autoproperty:: values

.. autoclass:: pytket.circuit.ClassicalEvalOp

   .. automethod:: __init__

.. autoclass:: pytket.circuit.ClassicalOp

   .. automethod:: __init__
   .. autoproperty:: n_input_outputs
   .. autoproperty:: n_inputs
   .. autoproperty:: n_outputs

.. autoclass:: pytket.circuit.CopyBitsOp

   .. automethod:: __init__

.. autoclass:: pytket.circuit.RangePredicateOp

   .. automethod:: __init__
   .. autoproperty:: lower
   .. autoproperty:: upper

.. autoclass:: pytket.circuit.MultiplexorBox

   .. automethod:: __init__
   .. automethod:: get_bitstring_op_pair_list
   .. automethod:: get_circuit
   .. automethod:: get_op_map

.. autoclass:: pytket.circuit.MultiplexedRotationBox

   .. automethod:: __init__
   .. automethod:: get_bitstring_op_pair_list
   .. automethod:: get_circuit
   .. automethod:: get_op_map

.. autoclass:: pytket.circuit.MultiplexedU2Box

   .. automethod:: __init__
   .. automethod:: get_bitstring_op_pair_list
   .. automethod:: get_circuit
   .. automethod:: get_impl_diag
   .. automethod:: get_op_map

.. autoclass:: pytket.circuit.MultiplexedTensoredU2Box

   .. automethod:: __init__
   .. automethod:: get_bitstring_op_pair_list
   .. automethod:: get_circuit
   .. automethod:: get_op_map

.. autoclass:: pytket.circuit.StatePreparationBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_statevector
   .. automethod:: is_inverse
   .. automethod:: with_initial_reset

.. autoclass:: pytket.circuit.DiagonalBox

   .. automethod:: __init__
   .. automethod:: get_circuit
   .. automethod:: get_diagonal
   .. automethod:: is_upper_triangle

.. autoclass:: pytket.circuit.ConjugationBox

   .. automethod:: __init__
   .. automethod:: get_action
   .. automethod:: get_circuit
   .. automethod:: get_compute
   .. automethod:: get_uncompute

.. autoclass:: pytket.circuit.ResourceBounds

   .. automethod:: __init__
   .. automethod:: get_max
   .. automethod:: get_min

.. autoclass:: pytket.circuit.ResourceData

   .. automethod:: __init__
   .. automethod:: __repr__
   .. automethod:: get_gate_depth
   .. automethod:: get_op_type_count
   .. automethod:: get_op_type_depth
   .. automethod:: get_two_qubit_gate_depth

.. autoclass:: pytket.circuit.DummyBox

   .. automethod:: __init__
   .. automethod:: get_n_bits
   .. automethod:: get_n_qubits
   .. automethod:: get_resource_data
  
.. autoclass:: pytket.circuit.BarrierOp

   .. automethod:: __init__
   .. autoproperty:: data

.. autoclass:: pytket.circuit.MetaOp

   .. automethod:: __init__
   .. autoproperty:: data

.. automodule:: pytket.circuit.named_types
    :members:
