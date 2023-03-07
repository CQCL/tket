pytket.passes
==================================

In pytket, compilation passes perform in-place transformations of circuits. From a user's point of view, passes are similar to `transforms <https://cqcl.github.io/tket/pytket/api/transform.html#>`_; however passes allow for additional predicate checking and compositionality. 

There are passes such as `FullPeepholeOptimise <https://cqcl.github.io/tket/pytket/api/passes.html#pytket.passes.FullPeepholeOptimise>`_ and  `KAKDecomposition <https://cqcl.github.io/tket/pytket/api/passes.html#pytket.passes.KAKDecomposition>`_ which are designed for general purpose circuit optimisation.

Also there are special purpose passes such as `OptimisePhaseGadgets <https://cqcl.github.io/tket/pytket/api/passes.html#pytket.passes.OptimisePhaseGadgets>`_ and `PauliSimp <https://cqcl.github.io/tket/pytket/api/passes.html#pytket.passes.PauliSimp>`_ which perform optimisation by targeting phase gadget and Pauli gadget structures within circuits. For more on these optimisation techniques see the `corresponding publication <https://arxiv.org/abs/1906.01734>`_.

Rebase passes can be used to convert a circuit to a desired gateset. See `RebaseCustom <https://cqcl.github.io/tket/pytket/api/passes.html#pytket.passes.RebaseCustom>`_ and `auto_rebase_pass <https://cqcl.github.io/tket/pytket/api/passes.html#module-pytket.passes.auto_rebase>`_.

For more on pytket passes see the `compilation <https://cqcl.github.io/pytket/manual/manual_compiler.html>`_ section of the user manual or the following notebook tutorials:

1. `Compilation example <https://github.com/CQCL/pytket/blob/main/examples/compilation_example.ipynb>`_ - Introduces different pytket optimisation passes and shows how combine them along with predicate checking.
2. `Symbolics example <https://github.com/CQCL/pytket/blob/bca57d776a3e4fd497907b54b902fefdf39684d3/examples/symbolics_example.ipynb#L1>`_ - Covers `symbolic compilation <https://cqcl.github.io/pytket/manual/manual_compiler.html#compiling-symbolic-circuits>`_ in pytket ,a technique often used in variational experiments.
3. `Contextual optimisation <https://github.com/CQCL/pytket/blob/main/examples/contextual_optimization.ipynb>`_ - Covers contextual optimisation, an advanced tket feature. Also see the corresponding `manual section <https://cqcl.github.io/pytket/manual/manual_compiler.html#contextual-optimisations>`_. 
4. `UCC VQE example <https://github.com/CQCL/pytket/blob/main/examples/ucc_vqe.ipynb>`_ - Demonstrates various compilation passes in the context of a variational quantum eigensolver (VQE) experiment.

.. automodule:: pytket._tket.passes
    :members:
    :special-members: __init__

pytket.passes.script
~~~~~~~~~~~~~~~~~~~~

.. automodule:: pytket.passes.script
    :members: compilation_pass_from_script, compilation_pass_grammar

pytket.passes.auto_rebase
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.passes.auto_rebase
    :members: auto_rebase_pass, auto_squash_pass
    
