pytket.passes
==================================

In pytket, compilation passes perform in-place transformations of circuits. From a user's point of view, passes are similar to `transforms <https://docs.quantinuum.com/tket/api-docs/transform.html>`_; however passes allow for additional predicate checking and compositionality. 

There are passes such as `FullPeepholeOptimise <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket.passes.FullPeepholeOptimise>`_ and  `KAKDecomposition <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket.passes.KAKDecomposition>`_ which are designed for general purpose circuit optimisation.

Also there are special purpose passes such as `OptimisePhaseGadgets <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket.passes.OptimisePhaseGadgets>`_ and `PauliSimp <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket.passes.PauliSimp>`_ which perform optimisation by targeting phase gadget and Pauli gadget structures within circuits. For more on these optimisation techniques see the `corresponding publication <https://arxiv.org/abs/1906.01734>`_.

Rebase passes can be used to convert a circuit to a desired gateset. See `RebaseCustom <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket.passes.RebaseCustom>`_ and `AutoRebase <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket._tket.passes.AutoRebase>`_.

For more on pytket passes see the `compilation <https://docs.quantinuum.com/tket/user-guide/manual/manual_compiler.html>`_ section of the user manual or the `notebook tutorials <https://docs.quantinuum.com/tket/examples>`_


.. automodule:: pytket._tket.passes
    :members:
    :special-members: __init__

.. autofunction:: pytket.passes.scratch_reg_resize_pass

.. autoclass:: pytket.passes.PassSelector
    :special-members: __init__
    :members:

pytket.passes.script
~~~~~~~~~~~~~~~~~~~~

.. automodule:: pytket.passes.script
    :members: compilation_pass_from_script, compilation_pass_grammar
