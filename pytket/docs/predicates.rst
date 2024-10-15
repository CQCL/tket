pytket.predicates
==================================

In pytket, predicates enforce properties of circuits. Each pytket :py:class:`Backend` has its own set of predicates which must be satisfied before a quantum circuit can be executed. There are predicates that enforce restrictions including gateset, number of qubits and classical control.

For more on predicates read the corresponding section of the `user manual <https://docs.quantinuum.com/tket/user-guide/manual/manual_compiler.html#compilation-predicates>`_. See also the `Compilation example <https://docs.quantinuum.com/tket/user-guide/examples/circuit_compilation/compilation_example.html>`_ notebook.

.. automodule:: pytket._tket.predicates
    :members:
    :special-members: __init__
