pytket.passes
==================================

In pytket, compilation passes perform in-place transformations of circuits. From a user's point of view, passes are similar to `transforms <https://cqcl.github.io/tket/pytket/api/transform.html#>`_; however passes allow for additional predicate checking and compositionality. 

There are passes such as `FullPeepholeOptimise <https://cqcl.github.io/tket/pytket/api/passes.html#pytket.passes.FullPeepholeOptimise>`_ and  `KAKDecomposition <https://cqcl.github.io/tket/pytket/api/passes.html#pytket.passes.KAKDecomposition>`_ which are designed for general purpose circuit optimisation.

Also there are special purpose passes such as `OptimisePhaseGadgets <https://cqcl.github.io/tket/pytket/api/passes.html#pytket.passes.OptimisePhaseGadgets>`_ and `PauliSimp <https://cqcl.github.io/tket/pytket/api/passes.html#pytket.passes.PauliSimp>`_ which perform optimisation by targeting phase gadget and Pauli gadget structures within circuits. For more on these optimisation techniques see the `corresponding publication <https://arxiv.org/abs/1906.01734>`_.

Rebase passes can be used to convert a circuit to a desired gateset. See `RebaseCustom <https://cqcl.github.io/tket/pytket/api/passes.html#pytket.passes.RebaseCustom>`_ and `auto_rebase_pass <https://cqcl.github.io/tket/pytket/api/passes.html#module-pytket.passes.auto_rebase>`_.

For more on pytket passes see the `compilation <https://cqcl.github.io/pytket/manual/manual_compiler.html>`_ section of the user manual or the `notebook tutorials <https://tket.quantinuum.com/examples>`_


.. automodule:: pytket._tket.passes
    :members:
    :special-members: __init__

.. automodule:: pytket.passes
    :members: PassSelector

pytket.passes.script
~~~~~~~~~~~~~~~~~~~~

.. automodule:: pytket.passes.script
    :members: compilation_pass_from_script, compilation_pass_grammar

pytket.passes.auto_rebase
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.passes.auto_rebase
    :members: auto_rebase_pass, auto_squash_pass
    
