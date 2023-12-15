pytket.circuit.Circuit
==================================
:py:class:`Circuit` objects provide an abstraction of quantum circuits. They consist of a set of qubits/quantum wires and a collection of operations applied to them in a given order. These wires have open inputs and outputs, rather than assuming any fixed input state.

See the `pytket User Manual <https://tket.quantinuum.com/user-manual/manual_circuit.html>`_ for a step-by-step tutorial on constructing circuits.

See also the notebook tutorials on `circuit generation <https://github.com/CQCL/pytket/blob/main/examples/circuit_generation_example.ipynb>`_ and `circuit analysis <https://github.com/CQCL/pytket/blob/main/examples/circuit_analysis_example.ipynb>`_.


Many of the :py:class:`Circuit` methods described below append a gate or box to
the end of the circuit. Where ``kwargs`` are indicated in these methods, the
following keyword arguments are supported:

- ``opgroup`` (:py:class:`str`): name of the associated operation group, if any
- ``condition`` (:py:class:`Bit`, :py:class:`BitLogicExp` or :py:class:`Predicate`): classical condition for applying operation
- ``condition_bits`` (list of :py:class:`Bit`): classical bits on which to condition operation
- ``condition_value`` (:py:class:`int`): required value of condition bits (little-endian), defaulting to all-1s if not specified

(Thus there are two ways to express classical conditions: either using a general
``condition``, or using the pair ``condition_bits`` and ``condition_value`` to
condition on a specified set of bit values.)

.. currentmodule:: pytket.circuit.Circuit
.. autoclass:: pytket.circuit.Circuit
   :special-members: __init__, __eq__, __iter__, __mul__, __repr__, __rshift__, __str__, __getstate__, __setstate__, __hash__

   .. automethod:: add_gate

   .. automethod:: append

   .. automethod:: add_circuit

   .. automethod:: add_qubit

   .. automethod:: add_bit

   .. automethod:: add_q_register

   .. automethod:: get_q_register

   .. automethod:: add_c_register

   .. automethod:: get_c_register

   .. automethod:: add_phase

   .. automethod:: rename_units

   .. automethod:: add_blank_wires

   .. automethod:: remove_blank_wires

   .. automethod:: flatten_registers

   .. autoproperty:: n_qubits

   .. autoproperty:: n_bits

   .. autoproperty:: phase

   .. autoproperty:: qubits

   .. autoproperty:: bits

   .. autoproperty:: q_registers

   .. autoproperty:: c_registers

   .. autoproperty:: n_gates

   .. autoproperty:: qubit_readout

   .. autoproperty:: opgroups

   .. autoproperty:: is_simple

   .. automethod:: commands_of_type

   .. automethod:: ops_of_type

   .. automethod:: n_gates_of_type

   .. automethod:: n_1qb_gates

   .. automethod:: n_2qb_gates

   .. automethod:: n_nqb_gates

   .. automethod:: free_symbols

   .. automethod:: depth

   .. automethod:: depth_2q

   .. automethod:: depth_by_type 

   .. automethod:: get_commands

   .. automethod:: add_barrier

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

   .. automethod:: add_classicalexpbox_bit

   .. automethod:: add_classicalexpbox_register

   .. automethod:: qubit_create

   .. automethod:: qubit_create_all

   .. automethod:: qubit_discard

   .. automethod:: qubit_discard_all

   .. automethod:: qubit_is_created
   
   .. automethod:: qubit_is_discarded
   
   .. automethod:: replace_SWAPs

   .. automethod:: replace_implicit_wire_swaps

   .. automethod:: implicit_qubit_permutation

   .. automethod:: to_latex_file


   Convenience methods for appending circuit operations

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

   .. automethod:: Measure

   .. automethod:: measure_all

   .. automethod:: measure_register

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

   Methods for adding circuit boxes
   --------------------------------

   .. Note:: For adding boxes to a circuit the :py:meth:`Circuit.add_gate` method is sufficient to append any :py:class:`OpType` to a :py:class:`Circuit`.

   .. jupyter-input::

      from pytket.circuit import Circuit, CircBox
     
      sub_circ = Circuit(2).CX(0, 1).Rz(0.25, 1).CX(0, 1)

      box = CircBox(sub_circ)

      bigger_circ = Circuit(3)

      # Equivalent to bigger_circ.add_circbox(box, [0, 1, 2])
      bigger_circ.add_gate(box, [0, 1, 2])

   .. automethod:: add_circbox

   .. automethod:: add_unitary1qbox
   
   .. automethod:: add_unitary2qbox
   
   .. automethod:: add_unitary3qbox

   .. automethod:: add_expbox

   .. automethod:: add_pauliexpbox

   .. automethod:: add_pauliexppairbox

   .. automethod:: add_pauliexpcommutingsetbox

   .. automethod:: add_toffolibox
   
   .. automethod:: add_dummybox

   .. automethod:: add_qcontrolbox
   
   .. automethod:: add_custom_gate

   .. automethod:: add_assertion

   .. automethod:: add_multiplexor

   .. automethod:: add_multiplexedrotation
   
   .. automethod:: add_multiplexedu2
      
   .. automethod:: add_state_preparation_box

   .. automethod:: add_diagonal_box

   .. automethod:: add_conjugation_box
   

