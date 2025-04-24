pytket.passes
==================================

In pytket, compilation passes perform in-place transformations of circuits. From a user's point of view, passes are similar to `transforms <https://docs.quantinuum.com/tket/api-docs/transform.html>`_; however passes allow for additional predicate checking and compositionality. 

There are passes such as `FullPeepholeOptimise <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket.passes.FullPeepholeOptimise>`_ and  `KAKDecomposition <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket.passes.KAKDecomposition>`_ which are designed for general purpose circuit optimisation.

Also there are special purpose passes such as `OptimisePhaseGadgets <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket.passes.OptimisePhaseGadgets>`_ and `PauliSimp <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket.passes.PauliSimp>`_ which perform optimisation by targeting phase gadget and Pauli gadget structures within circuits. For more on these optimisation techniques see the `corresponding publication <https://arxiv.org/abs/1906.01734>`_.

Rebase passes can be used to convert a circuit to a desired gateset. See `RebaseCustom <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket.passes.RebaseCustom>`_ and `AutoRebase <https://docs.quantinuum.com/tket/api-docs/passes.html#pytket._tket.passes.AutoRebase>`_.

For more on pytket passes see the `compilation <https://docs.quantinuum.com/tket/user-guide/manual/manual_compiler.html>`_ section of the user manual or the `notebook tutorials <https://docs.quantinuum.com/tket/examples>`_


.. currentmodule:: pytket._tket.passes

.. autoclass:: pytket.passes.BasePass

   .. automethod:: __init__
   .. automethod:: apply
   .. automethod:: from_dict
   .. automethod:: get_gate_set
   .. automethod:: get_postconditions
   .. automethod:: get_preconditions
   .. automethod:: to_dict

.. autoenum:: pytket.passes.CNotSynthType

.. autoclass:: pytket.passes.RepeatPass

   .. automethod:: __init__
   .. automethod:: get_pass

.. autoclass:: pytket.passes.RepeatUntilSatisfiedPass

   .. automethod:: __init__
   .. automethod:: get_pass
   .. automethod:: get_predicate

.. autoclass:: pytket.passes.RepeatWithMetricPass

   .. automethod:: __init__
   .. automethod:: get_metric
   .. automethod:: get_pass

.. autoenum:: pytket.passes.SafetyMode

.. autoclass:: pytket.passes.SequencePass

   .. automethod:: __init__
   .. automethod:: get_sequence
   .. automethod:: to_dict

.. automethod:: pytket.passes.AASRouting
.. automethod:: pytket.passes.AutoRebase
.. automethod:: pytket.passes.AutoSquash
.. automethod:: pytket.passes.CXMappingPass
.. automethod:: pytket.passes.CliffordPushThroughMeasures
.. automethod:: pytket.passes.CliffordResynthesis
.. automethod:: pytket.passes.CliffordSimp
.. automethod:: pytket.passes.CnXPairwiseDecomposition
.. automethod:: pytket.passes.CommuteThroughMultis
.. automethod:: pytket.passes.ComposePhasePolyBoxes
.. automethod:: pytket.passes.ContextSimp
.. automethod:: pytket.passes.CustomPass
.. automethod:: pytket.passes.CustomPassMap
.. automethod:: pytket.passes.CustomRoutingPass
.. automethod:: pytket.passes.DecomposeArbitrarilyControlledGates
.. automethod:: pytket.passes.DecomposeBoxes
.. automethod:: pytket.passes.DecomposeClassicalExp
.. automethod:: pytket.passes.DecomposeMultiQubitsCX
.. automethod:: pytket.passes.DecomposeSingleQubitsTK1
.. automethod:: pytket.passes.DecomposeSwapsToCXs
.. automethod:: pytket.passes.DecomposeSwapsToCircuit
.. automethod:: pytket.passes.DecomposeTK2
.. automethod:: pytket.passes.DefaultMappingPass
.. automethod:: pytket.passes.DelayMeasures
.. automethod:: pytket.passes.EulerAngleReduction
.. automethod:: pytket.passes.FlattenRegisters
.. automethod:: pytket.passes.FlattenRelabelRegistersPass
.. automethod:: pytket.passes.FullMappingPass
.. automethod:: pytket.passes.FullPeepholeOptimise
.. automethod:: pytket.passes.GreedyPauliSimp
.. automethod:: pytket.passes.GuidedPauliSimp
.. automethod:: pytket.passes.KAKDecomposition
.. automethod:: pytket.passes.NaivePlacementPass
.. automethod:: pytket.passes.NormaliseTK2
.. automethod:: pytket.passes.OptimisePhaseGadgets
.. automethod:: pytket.passes.PauliExponentials
.. automethod:: pytket.passes.PauliSimp
.. automethod:: pytket.passes.PauliSquash
.. automethod:: pytket.passes.PeepholeOptimise2Q
.. automethod:: pytket.passes.PlacementPass
.. automethod:: pytket.passes.RebaseCustom
.. automethod:: pytket.passes.RebaseTket
.. automethod:: pytket.passes.RemoveBarriers
.. automethod:: pytket.passes.RemoveDiscarded
.. automethod:: pytket.passes.RemoveImplicitQubitPermutation
.. automethod:: pytket.passes.RemovePhaseOps
.. automethod:: pytket.passes.RemoveRedundancies
.. automethod:: pytket.passes.RenameQubitsPass
.. automethod:: pytket.passes.RoundAngles
.. automethod:: pytket.passes.RoutingPass
.. automethod:: pytket.passes.SimplifyInitial
.. automethod:: pytket.passes.SimplifyMeasured
.. automethod:: pytket.passes.SquashCustom
.. automethod:: pytket.passes.SquashRzPhasedX
.. automethod:: pytket.passes.SquashTK1
.. automethod:: pytket.passes.SynthesiseTK
.. automethod:: pytket.passes.SynthesiseTket
.. automethod:: pytket.passes.ThreeQubitSquash
.. automethod:: pytket.passes.ZXGraphlikeOptimisation
.. automethod:: pytket.passes.ZZPhaseToRz

.. autofunction:: pytket.passes.scratch_reg_resize_pass

.. autoclass:: pytket.passes.PassSelector

   .. automethod:: __init__
   .. automethod:: apply
   .. automethod:: get_scores

pytket.passes.script
~~~~~~~~~~~~~~~~~~~~

.. automodule:: pytket.passes.script
    :members: compilation_pass_from_script, compilation_pass_grammar
