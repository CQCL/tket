pytket.circuit.Circuit
======================
:py:class:`Circuit` objects provide an abstraction of quantum circuits. They consist of a set of qubits/quantum wires and a collection of operations applied to them in a given order. These wires have open inputs and outputs, rather than assuming any fixed input state.

See the `pytket User Manual <https://docs.quantinuum.com/tket/user-guide/manual/manual_circuit.html>`_ for a step-by-step tutorial on constructing circuits.

See also the notebook tutorials on `circuit generation <https://docs.quantinuum.com/tket/user-guide/examples/circuit_construction/circuit_generation_example.html>`_ and `circuit analysis <https://docs.quantinuum.com/tket/user-guide/examples/circuit_construction/circuit_analysis_example.html>`_.


Many of the :py:class:`~pytket.circuit.Circuit` methods described below append a gate or box to
the end of the circuit. Where ``kwargs`` are indicated in these methods, the
following keyword arguments are supported:

- ``opgroup`` (:py:class:`str`): name of the associated operation group, if any
- ``condition`` (:py:class:`Bit`, :py:class:`BitLogicExp` or :py:class:`Predicate`): classical condition for applying operation
- ``condition_bits`` (list of :py:class:`Bit`): classical bits on which to condition operation
- ``condition_value`` (:py:class:`int`): required value of condition bits (little-endian), defaulting to all-1s if not specified

(Thus there are two ways to express classical conditions: either using a general
``condition``, or using the pair ``condition_bits`` and ``condition_value`` to
condition on a specified set of bit values.)

..
   Sphinx doesn't seem to offer much control over how the methods and properties are ordered in the docs
   We list the methods and properties manually (for now) to ensure that the most important methods (e.g. Circuit.add_gate) 
   are closer to the top of the page. Some less important methods like Circuit.YYPhase
   and Circuit.add_multiplexed_tensored_u2 are near the bottom.
   Since the Circuit class is so big we use the check_circuit_class_docs.py script to check that we haven't missed anything. 

.. currentmodule:: pytket.circuit.Circuit
.. autoclass:: pytket.circuit.Circuit

   .. automethod:: __init__
   
   .. automethod:: __iter__

   .. automethod:: __rshift__

   .. automethod:: add_gate

   .. automethod:: append

   .. automethod:: add_circuit

   .. automethod:: add_circuit_with_map

   .. automethod:: add_qubit

   .. automethod:: add_bit

   .. automethod:: add_phase

   .. automethod:: add_clexpr_from_logicexp

   .. automethod:: wasm_uid

   .. autoproperty:: name

   .. autoproperty:: n_qubits

   .. autoproperty:: n_bits

   .. autoproperty:: phase

   .. autoproperty:: qubits

   .. autoproperty:: bits

   .. autoproperty:: n_gates

   .. autoproperty:: is_symbolic

   .. autoproperty:: has_implicit_wireswaps

   .. automethod:: add_q_register

   .. automethod:: get_q_register

   .. automethod:: add_c_register

   .. automethod:: get_c_register

   .. autoproperty:: q_registers

   .. autoproperty:: c_registers

   .. automethod:: rename_units

   .. automethod:: add_blank_wires

   .. automethod:: remove_blank_wires

   .. automethod:: flatten_registers

   .. autoproperty:: qubit_readout

   .. autoproperty:: bit_readout

   .. autoproperty:: opgroups

   .. autoproperty:: is_simple

   .. autoproperty:: qubit_to_bit_map

   .. automethod:: commands_of_type

   .. automethod:: ops_of_type

   .. automethod:: n_gates_of_type

   .. automethod:: n_1qb_gates

   .. automethod:: n_2qb_gates

   .. automethod:: n_nqb_gates

   .. automethod:: free_symbols

   .. automethod:: symbol_substitution

   .. automethod:: substitute_named

   .. automethod:: depth

   .. automethod:: depth_2q

   .. automethod:: depth_by_type

   .. automethod:: get_commands

   .. automethod:: add_barrier

   .. automethod:: add_conditional_barrier

   .. automethod:: add_wasm

   .. automethod:: add_wasm_to_reg

   .. automethod:: from_dict

   .. automethod:: to_dict

   .. automethod:: get_statevector

   .. automethod:: get_unitary

   .. automethod:: get_unitary_times_other

   .. automethod:: dagger

   .. automethod:: transpose

   .. automethod:: copy

   .. automethod:: get_resources

   .. automethod:: add_c_and

   .. automethod:: add_c_not

   .. automethod:: add_c_or

   .. automethod:: add_c_xor
   
   .. automethod:: add_c_range_predicate

   .. automethod:: add_c_and_to_registers

   .. automethod:: add_c_or_to_registers

   .. automethod:: add_c_xor_to_registers

   .. automethod:: add_c_not_to_registers

   .. automethod:: add_c_copybits

   .. automethod:: add_c_copyreg

   .. automethod:: add_c_setreg

   .. automethod:: add_c_setbits

   .. automethod:: add_c_transform

   .. automethod:: add_c_modifier
      
   .. automethod:: add_c_predicate

   .. automethod:: add_clexpr

   .. automethod:: qubit_create

   .. automethod:: qubit_create_all

   .. automethod:: qubit_discard

   .. automethod:: qubit_discard_all

   .. automethod:: qubit_is_created
   
   .. automethod:: qubit_is_discarded

   .. autoproperty:: created_qubits

   .. autoproperty:: discarded_qubits

   .. autoproperty:: valid_connectivity
   
   .. automethod:: replace_SWAPs

   .. automethod:: replace_implicit_wire_swaps

   .. automethod:: implicit_qubit_permutation

   .. automethod:: to_latex_file


   Convenience methods for appending gates
   ---------------------------------------

   .. Note:: For adding gates to a circuit the :py:meth:`Circuit.add_gate` method is sufficient to append any :py:class:`OpType` to a :py:class:`Circuit`. 
      Some gates can only be added with :py:meth:`Circuit.add_gate`. For other more commonly used operations these can be added to a :py:class:`Circuit` directly using the convenience methods below.

   .. automethod:: H

   .. automethod:: X

   .. automethod:: Y

   .. automethod:: Z

   .. automethod:: S

   .. automethod:: Sdg

   .. automethod:: SX

   .. automethod:: SXdg

   .. automethod:: T

   .. automethod:: Tdg

   .. automethod:: V

   .. automethod:: Vdg

   .. automethod:: Rx

   .. automethod:: Ry

   .. automethod:: Rz

   .. automethod:: PhasedX

   .. automethod:: TK1

   .. automethod:: TK2

   .. automethod:: U1

   .. automethod:: U2
   
   .. automethod:: U3

   .. automethod:: CX

   .. automethod:: CY

   .. automethod:: CZ

   .. automethod:: CS

   .. automethod:: CSdg

   .. automethod:: CV

   .. automethod:: CVdg

   .. automethod:: CSX

   .. automethod:: CSXdg

   .. automethod:: CH

   .. automethod:: ECR

   .. automethod:: CRx

   .. automethod:: CRy

   .. automethod:: CRz

   .. automethod:: CU1

   .. automethod:: CU3

   .. automethod:: Measure

   .. automethod:: measure_all

   .. automethod:: measure_register

   .. automethod:: Reset

   .. automethod:: Phase
   
   .. automethod:: SWAP

   .. automethod:: CCX

   .. automethod:: CSWAP

   .. automethod:: ESWAP

   .. automethod:: ISWAP

   .. automethod:: ISWAPMax
   
   .. automethod:: PhasedISWAP

   .. automethod:: FSim

   .. automethod:: Sycamore
   
   .. automethod:: XXPhase

   .. automethod:: XXPhase3

   .. automethod:: YYPhase

   .. automethod:: ZZPhase

   .. automethod:: ZZMax

   .. automethod:: AAMS

   .. automethod:: GPI

   .. automethod:: GPI2

   Methods for appending circuit boxes
   -----------------------------------

   .. Note:: For adding boxes to a circuit the :py:meth:`Circuit.add_gate` method is sufficient to append any :py:class:`OpType` to a :py:class:`Circuit`.

   ::

      from pytket.circuit import Circuit, CircBox
     
      sub_circ = Circuit(2).CX(0, 1).Rz(0.25, 1).CX(0, 1)

      box = CircBox(sub_circ)

      bigger_circ = Circuit(3)

      # Equivalent to bigger_circ.add_circbox(box, [0, 1, 2])
      bigger_circ.add_gate(box, [0, 1, 2])

   .. automethod:: add_circbox

   .. automethod:: add_circbox_regwise

   .. automethod:: add_circbox_with_regmap

   .. automethod:: add_unitary1qbox
   
   .. automethod:: add_unitary2qbox
   
   .. automethod:: add_unitary3qbox

   .. automethod:: add_expbox

   .. automethod:: add_pauliexpbox

   .. automethod:: add_pauliexppairbox

   .. automethod:: add_pauliexpcommutingsetbox

   .. automethod:: add_termsequencebox

   .. automethod:: add_phasepolybox

   .. automethod:: add_toffolibox
   
   .. automethod:: add_dummybox

   .. automethod:: add_qcontrolbox
   
   .. automethod:: add_custom_gate

   .. automethod:: add_assertion

   .. automethod:: add_multiplexor

   .. automethod:: add_multiplexedrotation
   
   .. automethod:: add_multiplexedu2

   .. automethod:: add_multiplexed_tensored_u2
      
   .. automethod:: add_state_preparation_box

   .. automethod:: add_diagonal_box

   .. automethod:: add_conjugation_box
   

