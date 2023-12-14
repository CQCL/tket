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

.. currentmodule:: pytket._tket.circuit.Circuit
.. autoclass:: pytket._tket.circuit.Circuit
   :special-members: __init__, __eq__, __iter__, __mul__, __repr__, __rshift__, __str__, __getstate__, __setstate__, __hash__

   .. automethod:: add_gate

   .. automethod:: add_circuit

   .. automethod:: add_bit

   .. automethod:: add_qubit

   .. automethod:: add_phase
   
   .. automethod:: add_blank_wires

   .. automethod:: commands_of_type

   .. automethod:: depth
   
   .. automethod:: depth_2q

   .. automethod:: depth_by_type 

   .. automethod:: add_barrier

   .. automethod:: flatten_registers

   .. automethod:: from_dict

   .. automethod:: to_dict

   .. automethod:: free_symbols

   .. automethod:: get_commands

   .. automethod:: get_statevector

   .. automethod:: get_unitary

   .. automethod:: get_unitary_times_other

   .. automethod:: add_c_and

   .. automethod:: add_c_and_to_registers

   .. authomethod:: add_c_copybits


   
