pytket.predicates
==================================

In pytket, predicates enforce properties of circuits. Each pytket :py:class:`Backend` has its own set of predicates which must be satisfied before a quantum circuit can be executed. There are predicates that enforce restrictions including gateset, number of qubits and classical control.

For more on predicates read the corresponding section of the `user manual <https://docs.quantinuum.com/tket/user-guide/manual/manual_compiler.html#compilation-predicates>`_. See also the `Compilation example <https://docs.quantinuum.com/tket/user-guide/examples/circuit_compilation/compilation_example.html>`_ notebook.

.. currentmodule:: pytket._tket.predicates

.. autoclass:: pytket.predicates.CliffordCircuitPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.CommutableMeasuresPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.CompilationUnit

   .. automethod:: __init__
   .. automethod:: check_all_predicates
   .. autoproperty:: circuit
   .. autoproperty:: final_map
   .. autoproperty:: initial_map

.. autoclass:: pytket.predicates.ConnectivityPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.DefaultRegisterPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.DirectednessPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.GateSetPredicate

   .. automethod:: __init__
   .. autoproperty:: gate_set

.. autoclass:: pytket.predicates.MaxNClRegPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.MaxNQubitsPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.MaxTwoQubitGatesPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.NoBarriersPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.NoClassicalBitsPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.NoClassicalControlPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.NoFastFeedforwardPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.NoMidMeasurePredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.NoSymbolsPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.NoWireSwapsPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.NormalisedTK2Predicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.PlacementPredicate

   .. automethod:: __init__

.. autoclass:: pytket.predicates.Predicate

   .. automethod:: __init__
   .. automethod:: from_dict
   .. automethod:: implies
   .. automethod:: to_dict
   .. automethod:: verify

.. autoclass:: pytket.predicates.UserDefinedPredicate

   .. automethod:: __init__
