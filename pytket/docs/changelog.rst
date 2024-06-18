Changelog
=========

1.29.2 (June 2024)
------------------

Feature:

* Revert keeping of blank classical wires when running
  ``FlattenRelabelRegistersPass``.

1.29.1 (June 2024)
------------------

Features:

* Improve depth of circuit produced by ``MultiplexedTensoredU2Box``.
* Revert support of classical transforms and predicates, and QASM registers,
  with up to 64 bits. (Revert maximum width to 32.)

1.29.0 (June 2024)
------------------

Features:

* Add ``OpType.CnRx`` and ``OpType.CnRz``.
* Add ``AutoRebase`` and ``AutoSquash`` passes.
  Deprecate ``auto_rebase_pass`` and ``auto_squash_pass``.
* Add new parameter to `remove_blank_wires` to allow to keep empty classical bits
* Support classical transforms and predicates, and QASM registers, with up to 64
  bits.

Fixes:

* Allow barriers when dagger or transpose a circuit.
* Keep blank classical wires when running `FlattenRelabelRegistersPass`
* Handle Clifford-angle ``NPhasedX`` gates in Clifford resynthesis.

1.28.0 (May 2024)
-----------------

Features:

* Update to pytket-circuit-renderer 0.8.
* Add two new status values for circuits on backends: "CANCELLING" and "RETRYING".
* Use `lark` package instead of deprecated `lark-parser`.
* Add ``GreedyPauliSimp`` optimisation pass.
* Add ``BitWiseOp.ZERO`` and ``BitWiseOp.ONE`` to allow construction of constant
  conditional expressions.
* Add target gateset ``(GPI, GPI2, AAMS)`` to ``auto_rebase_pass``.
* Add ``RebaseToIonQ`` transform.

Fixes:

* Escape underscores in qubit and bit names when converting to latex.

1.27.0 (April 2024)
-------------------

General:

* Remove deprecated ``SynthesiseHQS`` pass.

Features:

* Add ``circuit_name`` property to ``CircBox``.
* Enable pickling of ``Bit`` objects.
* New optimisation ``Transform.PushCliffordsThroughMeasures()`` and pass 
  ``CliffordPushThroughMeasures`` that optimises Clifford subcircuits 
  before end of circuit measurement gates.
* Add ``OpType.GPI``, ``OpType.GPI2`` and ``OpType.AAMS``.
* Allow construction of ``SequencePass`` without predicate checks, by means of
  new ``strict`` argument to the constructor (defaulting to ``True``).

Fixes:

* Correct handling of ``CustomGate`` when converting from pytket to QASM.
* Ensure that ECR, CS and CSdg operations have gate definitions in QASM
  conversion.
* Correct position of custom gate definitions needed for conditional operations
  in QASM conversion.
* Fix ``DelayMeasures()`` pass for circuits where bits are reused as measurement
  targets.
* When adding operations to a circuit, check for invalid wires before adding a
  vertex to the circuit.
* Make ``RemoveRedundancies`` pass remove ``OpType.Phase`` gates.
* Remove support for wasm functions with multiple return values.

Deprecations:

* Deprecate ``SynthesiseOQC`` pass.

1.26.0 (March 2024)
-------------------

Features:

* Allow ``CircBox`` containing non-default registers.
* Add new methods ``Circuit.add_circbox_regwise()`` and
  ``Circuit.add_circbox_with_regmap()`` for adding a ``CircBox`` to a circuit
  providing either an ordered sequence of registers or a mapping of registers
  from the box to the containing circuit.
* Add ``CliffordResynthesis`` pass to apply Clifford resynthesis (optionally
  with a user-defined resynthesis method) on all Clifford subcircuits.
* Add optional ``min_p`` argument to
  ``BackendResult.get_probability_distribution()`` and to the constructor of a
  ``ProbabilityDistribution``, defaulting to zero. (Previously probabilities
  below 1e-10 were by default treated as zero.)
* Add python binding for ``UnitaryRevTableau``.
* Add ``TermSequenceBox``, for circuit synthesis of a series of Pauli 
  Exponentials, where the ordering of terms can be changed.

Fixes:

* Add missing op types to methods for converting Clifford circuits to unitary
  tableaux.
* Require scipy >= 1.13 and quimb >= 1.8 for ZX module.

1.25.0 (February 2024)
----------------------

Features:

* Add ``WasmFileHandler.bytecode()`` method to retrieve the WASM as bytecode.

Fixes:

* Fix bug in ``PauliExponentials()`` pass affecting circuits containing
  ``PhasedX`` gates containing Clifford angles.

1.24.0 (January 2024)
---------------------

General:

* Python 3.12 support added; 3.9 dropped.

Features:

* Accept ``OpType.Phase`` in circuits passed to ``ZXGraphlikeOptimisation``.

Fixes:

* Handle a missing edge case in decomposition of single-qubit rotations.
* Add missing ``OpType.ConjugationBox``.

1.23.0 (January 2024)
---------------------

API changes:

* Make the ``architecture`` field in ``BackendInfo`` optional.

Deprecations:

* Deprecate ``SynthesiseHQS`` pass.
  
Fixes:

* Ensure that squashing long sequences of gates via unitary multiplication does
  not produce non-unitary results due to rounding errors.
* Fix `PauliFrameRandomisation.sample_circuits`.
* For `Circuit` with no 2-qubit gates, `NoiseAwarePlacement` now assigns `Qubit` to `Node` in `Architecture`
  with lowest reported error rates.
* Fix invalid registers returned by ``Circuit.q_registers`` and ``Circuit.c_registers``.
* Fix regression (introduced in 1.22.0) in compilation performance with certain
  sequences of passes.


1.22.0 (November 2023)
----------------------

Minor new features:

* Add optional parameter to QASM conversion methods to set the maximum allowed
  width of classical registers (default 32).
* New ``OpType.CS`` and ``OpType.CSdg``.
* New classes ``ResourceBounds``, ``ResourceData`` and ``DummyBox``, and method
  ``Circuit.get_resources()``, allowing reasoning about resource requirements
  on circuit templates.

Fixes:

* When converting QASM expressions to ``ClassicalExpBox``, preserve the ordering
  of the bits in the expression in the resulting ``cmd.args``
* Fix incorrect serialisation of ``PauliExpPairBox`` when the Pauli strings are of
  length 2.
* Fix incorrect controlled ``ConjugationBox`` handling.

General:

* Drop support for MacOS 11.

`Full changelog <https://github.com/CQCL/tket/compare/v1.21.0...v1.22.0>`_

1.21.0 (October 2023)
---------------------

Minor new features:

* Add optional ``strict_check`` parameter to ``RepeatPass`` to force stopping when
  the circuit is unchanged.
* Add optional parameters ``excluded_types`` and ``excluded_opgroups``
  to ``DecomposeBoxes``.
* More efficient decomposition for quantum controlled ``ConjugationBox``es.
* New ``PassSelector`` for automatically compiling with the best pass from a list
* ``PauliExpBox``, ``PauliExpPairBox``, and ``PauliExpCommutingSetBox`` are now
  decomposed into a single ``ConjugationBox``.
* Make ``SquashRzPhasedX`` pass always squash symbols.
* Add in-place symbol_substition method for ``CircBox``
* Add rendering support for 0-valued control-type gates.
* Typing improvements
* Make ``BitRegister`` and ``QubitRegister`` iterable

Fixes:

* Handle symbolic angles in ``ZZPhaseToRz`` pass.
* Bind ``sympy.exp()``.
* Ensure determinate command order for circuits containing Phase operations.

1.20.1 (September 2023)
-----------------------

Fixes:

* Fix ``Op.get_unitary()`` runtime error for non gate ``Op``s.
* Fix ``CliffordSimp`` slow runtime issue.
* Correct implementation of ``free_symbols()`` and ``symbol_substitution()`` for
  ``ConjugationBox``.
* Fix pytket-to-QASM conversion when individual bits of registers used in
  range predicates are later set.

1.20.0 (September 2023)
-----------------------

Fixes:

* Mixed up function index in wasm file check
* Fix handling of scratch bits in pytket-to-QASM conversion when the source bit
  for the scratch is overwritten before the scratch bit is used in a
  conditional.

Minor new features:

* ``Circuit.add_conditional_barrier``
* Add ``apply_clifford_basis_change_tensor`` method

API changes:

* barrier changed from MetaOp to be a BarrierOp


1.19.1 (September 2023)
-----------------------

Fixes:

* Fix `RebaseCustom()` rebasing of `TK2` gates.
* Correct implementation of `symbol_substitution()` for box types that cannot
  contain symbols.

1.19.0 (September 2023)
-----------------------

Major new features:

* Add ``ConjugationBox`` to express circuits that follow
  the compute-action-uncompute pattern.
* Added typing support for compiled modules

Minor new features:

* Implement equality checking for all boxes.
* Add ``Op.is_clifford`` to python binding.
* Single-qubit squashing ignores chains of symbolic gates if squashing them
  would increase the overall complexity of the expressions. This behaviour can
  be overridden using the ``always_squash_symbols`` parameter to
  ``SquashCustom``.
* Add ``control_state`` argument to ``QControlBox``.
* Add ``QubitPauliTensor`` (combining ``QubitPauliString`` with a complex
  coefficient) to python binding. This is incorporated into ``UnitaryTableau`` 
  row inspection for phase tracking.

Fixes:

* Allow ``BackendResult`` objects containing no results.

1.18.0 (August 2023)
--------------------

Minor new features:

* Add circuit method ``depth_2q``.
* Add ``allow_swaps`` parameter to ``auto_rebase_pass``.

Fixes:

* Fix slow ``Circuit.get_statevector()``.


1.17.1 (July 2023)
------------------

General:

* Fix issue with installing recent pytket versions on macos x86_64 in conda
  environments.

Minor new features:

* New constructor for ``ToffoliBox`` that allows switching between two decomposition strategies:
  ``ToffoliBoxSynthStrat.Matching`` and ``ToffoliBoxSynthStrat.Cycle``.
* Prefer ``ZZPhase`` to ``CX`` or ``ZZMax`` when using ``auto_rebase_pass()``.

1.17.0 (July 2023)
------------------

Minor new features:

* `Circuit.get_unitary()` and `Circuit.get_statevector()` now work for circuits
  containing boxes.
* New Box type `PauliExpPairBox`.
* New Box type `PauliExpCommutingSetBox`.
* New pass `PauliExponentials` that rewrites a circuit to a sequence of `PauliExpBox`,
  `PauliExpPairBox`, `PauliExpCommutingSetBox` and a Clifford circuit.

1.16.0 (June 2023)
------------------

Minor new features:

* Support ``allow_swaps`` parameter for ``PeepholeOptimise2Q``.
* Add missing add box methods that accept qubit indices as arguments.
* Add ``with_initial_reset`` parameter to ``StatePreparationBox`` to permit
  state preparation starting from unknown state.
* New method ``utils.stats.gate_counts`` to count gates of all types.

Fixes:

* Fix ``FlattenRegisters`` not updating ``ClassicalExpBox``.
* Fix missing default argument value to ``FlattenRelabelRegistersPass``.
* Fix ``auto_rebase_pass`` rebasing via TK2 even if CX is the only target 2q gate.
* Fix ``QControlBox`` not identifying SU(2) unitaries.

1.15.0 (May 2023)
-----------------

Major new features:

* Add new ``MultiplexedTensoredU2Box`` that synthesises multiplexed tensor product of U2 gates.

Minor new features:

* Add new ``MaxNClRegPredicate`` that checks that there are at most n classical
  registers in the circuit.
* Allow barriers in ``QControlBoxes``. Barriers are left in place.
* Add ``Circuit.TK1`` and ``Circuit.TK2`` methods that take ``Qubit`` arguments.
* Expose ``CircuitRenderer`` instance so users can set their own default options.
* QASM to circuit converters now recognise ``Rxxyyzz`` as ``OpType.TK2``. Circuit
  to QASM converters with the "hqslib1" header now map ``OpType.TK2`` to ``Rxxyyzz``.
* Add new transform ``round_angles`` and pass ``RoundAngles`` to remove angles
  below a threshold and/or round angles to a dyadic fraction of pi throughout a
  circuit.

Fixes:

* Fix bug in `get_operator_expectation_value()` computation when operator
  includes `Pauli.I` terms.
* Fix bug in routing code occurring in ``Circuits`` with qubit wires with no operations
  and some (other or same) qubits pre-labelled as "Node" from the ``Architecture`` being routed to.

1.14.0 (April 2023)
-------------------

Major new features:

* Support for ARM Linux platforms.
* Updated implementation of ``ToffoliBox`` utilising multiplexors
  for improved decomposition.
* Add new ``DiagonalBox`` that synthesises a diagonal unitary matrix
  into a sequence of multiplexed-Rz gates.

1.13.2 (March 2023)
-------------------

Minor new features:

* Update to networkx 3.
* Add "label" argument to ``SquareGrid``, ``RingArch`` and ``FullyConnected`` 
  ``Architecture`` classes to give custom name to constructed ``Node``.
* Add ``FlattenRelabelRegistersPass`` to remove empty quantum wires and relabel all
  qubits to a default register named after a passed label.

Fixes:

* Multiply symbolic parameters in auto-generated gate definitions by "/pi" in ``circuit_to_qasm_io``

1.13.1 (March 2023)
-------------------

Fixes:

* Throw error rather than abort when trying to add qubit or bit with existing name.

1.13.0 (March 2023)
-------------------

Major new features:

* New ``StatePreparationBox`` to prepare arbitrary quantum states.
* New WasmWire interface to keep all wasm operation in the initial order
* New ``ZXGraphlikeOptimisation`` compilation pass for optimising the circuit by
  simplifying in ZX calculus and extracting back out

Minor new features:

* New ``CommutableMeasuresPredicate`` predicate, added as precondition to the
  ``DelayMeasures`` pass.
* Added an ``allow_partial`` parameter to the ``DelayMeasures`` pass to delay
  the measurements as much as possible when they cannot be fully delayed to the
  end.
* Update to ``pytket-circuit-renderer`` 0.5.
* Support ``allow_swaps`` parameter for ``FullPeepholeOptimise`` even when
  targeting ``OpType.TK2``.

Fixes:

* ``DelayMeasures`` pass now correctly handles circuits with ``CircBox``es.
* ``get_op_map`` in multiplexor boxes return unhashable python dictionaries.


1.11.1 (January 2023)
---------------------

General:

* Support for MacOS >= 11.0 on both x86_64 and arm64.

1.11.0 (January 2023)
---------------------

Major new features:

* New boxes to implement multiplexor gates (i.e. uniformly controlled operations):
  ``MultiplexorBox``, ``MultiplexedRotationBox`` and ``MultiplexedU2Box``.

General:

* Python 3.11 support added; 3.8 dropped.

Minor new features:

* Circuit methods ``qubit_readout`` and ``qubit_to_bit_map`` now ignore barriers.
* New pass ``RemoveImplicitQubitPermutation``.
* ``PauliSimp`` pass accepts circuits containing implicit wire swaps.

Fixes:

* ``MultiGateReorderRoutingMethod`` raising unknown edge missing error.
* ``LexiRouteLabellingMethod`` hitting assertion during dynamic qubit allocation.
* ``PauliSimp`` pass preserves circuit name.

1.10.0 (December 2022)
----------------------

Minor new features:

* Add support for PhasedX gates in Pauli graph synthesis.

Fixes:

* Handle 0-qubit operations in connectivity check.
* Fix handling of Tdg, CY, ZZMax and Clifford-angle YYPhase gates in Pauli
  graph synthesis.
* Disallow conversion to QASM of operations conditioned on strict subregisters
  larger than one bit, or reordered registers.

1.9.1 (December 2022)
---------------------

Minor new features:

* New ``view_browser`` function for opening a browser with circuit render.

Fixes:

* Warn rather than abort when significant rounding errors are detected in
  TK2-to-CX rebase.
* Fix incorrect QASM output for ``OpType.CopyBits``.
* Fix incorrect QASM read in ``OpType.ZZPhase``.

1.9.0 (November 2022)
---------------------

Fixes:

* Rebase and synthesis passes now respect conditional phase, by adding
  conditional ``OpType.Phase`` operations to the rebased circuit. Any code that
  relies on the circuit having gates only in the specified gate set should be
  updated to handle ``OpType.Phase`` as well when conditional operations are
  present.
* A bug where the sequence of ``RoutingMethod`` used in ``DefaultMappingPass`` could 
  add a cycle to the ``Circuit`` DAG has been fixed.
* Fix support for ECR gate in QASM converters.

API changes:

* The default value of ``optimisation_level`` in ``Backend`` methods that have
  this parameter (such as ``get_compiled_circuit()``) has been changed from 1 to
  2.

Minor new features:

* Added shortcuts for adding ``U1``, ``U2``, ``U3``, ``TK1``, ``TK2``, ``CU1``, 
  ``CU3``, ``ISWAP``, ``PhasedISWAP``, ``ESWAP``, ``PhasedX``, ``FSim``, ``Sycamore``
  and ``ISWAPMax`` gates to a ``pytket`` ``Circuit``.
* New ``Circuit`` methods ``n_1qb_gates``, ``n_2qb_gates``, ``n_nqb_gates``.
* New ``EmpriricalDistribution`` and ``ProbabilityDistribution`` utility classes
  for manipulating distributions, and methods to extract them from
  ``BackendResult`` objects.

1.8.1 (November 2022)
---------------------

Fixes:

* Incorrect qasm filtering.
* Make graph placement work with multi-qubit barriers.

1.8.0 (November 2022)
---------------------

Minor new features:

* New ``OpType::Phase`` 0-qubit gate affecting global phase.
* New ``CnXPairwiseDecomposition`` pass.
* Allow ``QControlBox`` with implicit wire swaps to be decomposed.
* New ``Circuit`` methods ``replace_SWAPs`` and ``replace_implicit_wire_swaps``.

Fixes:

* Remove unused ``tk_SCRATCH_BIT`` registers from qasm output.
* Update the ``LogicExp`` in every ``ClassicalExpBox`` when calling ``Circuit.rename_units``.
* Fix the json schema for ``LinePlacement``
* Fix issue with ``QControlBox`` throwing error during decomposition
  if the controlled circuit contains identity gates.
* Fix issue with ``KAKDecomposition`` raising exception if the circuit contains ``ClassicalExpBox``.

1.7.3 (October 2022)
--------------------

Minor new features:

* New ``Circuit`` properties ``created_qubits`` and ``discarded_qubits``.
* Barrier operations inside QASM custom gates are now accepted.
* Added wasm functions will be checked if the signatures are supported

Fixes:

* Circuit equality check now takes into account qubit creations and qubit discards.
* Created qubits and discarded qubits are now shown in ``Circuit.__repr__`` and ``Circuit.to_dict``.
* Allow symbolic operations in initial simplification.
* Fix the json schema for compiler passes.
* Fix ``SquashRzPhasedX`` so it now preserves phase.

1.6.1 (September 2022)
----------------------

Minor new features:

* New ``OpType.CnY`` and ``OpType.CnZ``.
* Update ``DecomposeArbitrarilyControlledGates`` pass to decompose ``CnX``,
  ``CnY``, and ``CnZ`` gates.

Fixes:

* ``Circuit.get_unitary()`` and ``Circuit.get_statevector()`` now throw an error
  when the circuit contains measurements.
* Fix critical issue with compilation of circuits containing conditional gates.

1.6.0 (September 2022)
----------------------

* New ``ToffoliBox`` for constructing circuits that implement permutations of
  basis states.

1.5.2 (August 2022)
-------------------

Minor new features:

* Prefer `ZZPhase` in ``DecomposeTK2`` if it results in the same fidelity but
  fewer two-qubit gates.

* Add ``SquashRzPhasedX`` pass to squash single qubit gates into
  ``Rz`` and ``PhasedX`` gates while trying to commute ``Rz``s to the back. 

1.5.1 (August 2022)
-------------------

Minor new features:

* Improve ``FullPeepholeOptimise`` performance.

Fixes:

* Squash two-qubit circuits properly in ``FullPeepholeOptimise`` for parameter
  `target_2qb_gate=OpType.TK2`.
* Floating point inaccuracies in ``NormalisedTK2Predicate``.

1.5.0 (August 2022)
-------------------

Minor new features:

* Add support for TK2 gate in ``KAKDecomposition``.
* ``Transform.ThreeQubitSquash()`` can now use TK2 gates as an alternative to CX
  gates.
* ``Unitary3qBox.get_circuit()`` decomposes the circuit using (at most 15) TK2
  gates.
* New ``CustomPass()`` accepting a user-supplied circuit transformation
  function.
* ``measure_register`` now allows using an existing classical register
* Provide an additional ``RebaseCustom`` constructor that takes a
  TK2-replacement instead of a CX-replacement function.
* New ``int_dist_from_state`` function in ``pytket.utils.results`` to convert
  a statevector to the probability distribution over its indices.
* The precondition for ``CliffordSimp`` and ``KAKDecomposition`` has been relaxed
  to accept classical controlled operations. ``ThreeQubitSquash`` and ``FullPeepholeOptimise``
  now accept classical operations.
* Improve ``QControlBox`` decomposition.
* New ``allow_swaps`` flag in ``KAKDecomposition`` and ``DecomposeTK2`` to
  decompose two-qubit operations up to implicit wire swaps.
* Add support for TK2 gate in ``FullPeepholeOptimise``.

Fixes:

* ``FullPeepholeOptimise`` failure on conditional circuits.

1.4.3 (July 2022)
-----------------

Fixes:

* Further relax assertion in ``replace_TK2_2CX``.

1.4.2 (July 2022)
-----------------

Fixes:

* Relax assertion in replace_TK2_2CX to avoid crash due to rounding errors.

1.4.1 (July 2022)
-----------------

Minor new features:

* New ``NormalisedTK2Predicate`` predicate and ``NormaliseTK2`` pass.
* New ``ZZPhaseToRz`` pass.
* Circuit to QASM converters with the "hqslib1" header now fix ZZPhase angles
  to be between -1 and 1 half-turns.

Fixes:

* Ensure TK2 angles are normalised before decomposing TK2 gates in passes.

1.3.0 (June 2022)
-----------------

Minor new features:

* New ``circuit_to_zx`` function to convert ``Circuit`` to ``ZXDiagram``, and
  ``to_circuit`` to extract from a unitary diagram.
* New ``to_graphviz_str`` method for ``ZXDiagram`` to generate a source string
  that can be rendered by the ``graphviz`` package.
* New pass and transform `DecomposeTK2` to decompose TK2 gates using the
  approximate KAK decomposition.
* Pass and transform ``GlobalisePhasedX`` use fewer Rz rotations.
* Improved decomposition for CnX gates.

Fixes:

* Fix serialization of `BackendInfo` for `RingArch` and `FullyConnected`
  architectures.

1.2.2 (May 2022)
----------------

Minor new features:

* The ``GlobalisePhasedX`` transform and homonymous pass take a new optional
  ``squash`` parameter. ``squash=true`` (default) implements a new algorithm
  that significantly reduces the number of ``NPhasedX`` gates synthesised.
* New ``DecomposeNPhasedX`` transform and pass replaces all ``NPhasedX`` gates
  with single-qubit ``PhasedX`` gates.
* Extend range of Clifford operations recognized by
  ``CliffordCircuitPredicate``.
* New ``circuit_from_qasm_wasm`` function to parse QASM files containing
  external WASM calls.
* Faster QASM parsing, capable of parsing extended grammar.

1.2.1 (May 2022)
----------------

Minor new features:

* Added explicit constructors for various Python classes.
* New ``measure_register`` method for measuring registers.
* Added ``OpType.TK2``, a three-parameter two-qubit gate.
* New pass ``SynthesiseTK`` and transform ``OptimiseStandard`` to synthesize
  TK2 gates.
* Add ``Optype.WASM``, adding a classical wasm function call to the circuit
* Add optype for existing PhasePolyBox ``OpType.PhasePolyBox``

1.1.0 (April 2022)
------------------

Minor new features:

* new additional constructor for ``PhasePolyBox`` from a given ``Circuit``
* New compilation pass ``ComposePhasePolyBoxes`` for generating
  PhasePolyBoxes in a given circuit
* Add JSON serialization methods for ``Predicate``, ``MeasurementSetup`` and ``MeasurementBitMap``.
* Add ``NoBarriersPredicate``.

Fixes:

* Fix qubit order in ``QubitPauliOperator.to_sparse_matrix()``.
* Fix issue with "nan" values appearing after symbolic substitution following
  compilation of some symbolic circuits.
* ``PhasePolyBox`` constructor is not accepting invalid boxes anymore

1.0.1 (March 2022)
------------------

Fixes:

* Fix problem with unassigned ancilla qubits during mapping.

1.0.0 (March 2022)
------------------

API changes:

* ``Rebase<Target>`` and ``SquashHQS`` methods are removed. Specifically:

  * ``RebaseHQS``
  * ``RebaseProjectQ``
  * ``RebasePyZX``
  * ``RebaseQuil``
  * ``RebaseUMD``
  * ``RebaseUFR``
  * ``RebaseOQC``

* The deprecated ``QubitPauliString.to_dict`` method is removed. (Use the
  ``map`` property instead.)
* The deprecated ``Backend.compile_circuit`` method is removed. (Use
  ``get_compiled_circuit`` instead.)
* The ``routing`` module is removed.
* ``Placement``, ``LinePlacement``, ``GraphPlacement`` and ``NoiseAwarePlacement`` 
  are now imported from the ``placement`` module.
* ``Architecture``, ``SquareGrid``, ``RingArch`` and ``FullyConnected`` are now 
  imported from the ``architecture`` module.
* Methods for mapping logical to physical circuits are now available in the
  ``mapping`` module, with a new API and new functionality.
* The keyword parameter and property ``def`` is now called ``definition`` in 
  ``Circuit.add_custom_gate`` and ``CustomGateDef``.
* ``RebaseCustom`` takes one allowed gateset parameter rather than separate single qubit and multiqubit gatesets.
* The ``Backend.characterisation`` property is removed. (Use
  ``Backend.backend_info`` instead.)
* The ``QubitPauliOperator.from_OpenFermion`` and
  ``QubitPauliOperator.to_OpenFermion`` methods are removed.
* The ``pytket.program`` module is removed.
* The ``pytket.telemetry`` module is removed.

Major new features:

* New methods for mapping logical to physical circuits for some ``Architecture``.
  The new method will use a list of user-given methods, each of them suitable only 
  for a specific set of subcircuits. Users can add their own methods if they want to.
  All compiler passes in pytket are updated to use the new methods.
  The methods already given by pytket are ``LexiRouteRoutingMethod``,
  ``LexiLabellingMethod``, ``MultiGateReorderRoutingMethod``,
  ``AASRouteRoutingMethod``, ``BoxDecompositionRoutingMethod``, and ``AASLabellingMethod``.

Minor new features:

* Add ``delay_measures`` option to ``DefaultMappingPass``.
* New ``pytket.passes.auto_rebase_pass`` and ``pytket.passes.auto_squash_pass``
  which attempt to construct rebase and squash passess given a target gate set from known
  decompositions.
* Add ``get_c_register``, ``get_q_register``, ``c_registers`` and ``q_registers`` methods to ``Circuit``.
* New ``pytket.passes.NaivePlacementPass`` which completes a basic relabelling of all Circuit Qubit
  not labelled as some Architecture Node to any available Architecture Node
* Add ``opgroups`` property to ``Circuit``.
* ``Architecture`` has new ``valid_operation`` method which returns true if passed UnitIDs that respect 
  architecture constraints.
* ``CircuitStatus`` has several new optional properties such as time-stamps associated with status changes,
  queue position or detailed error information.

Fixes:

* ``ConnectivityPredicate.implies()`` checks for existence of isolated nodes as
  well as edges in second architecture.
  
0.19.2 (February 2022)
----------------------

Fixes:

* Fix issue with jinja2 by updating dependency.

0.19.1 (February 2022)
----------------------

Fixes:

* Fix regression in ``Circuit.symbol_substitution`` causing incorrect values to
  be substituted in some cases.

0.19.0 (February 2022)
----------------------

Major new features:

* New box types for Clifford tableaux.

Minor new features:

* Improve ``CnX`` gate decomposition for n=5,6,7.
* Add ``rebase_pass`` method to ``Backend``.
* Add ``is_clifford_type`` method to ``Op``.

General:

* Python 3.10 support added; 3.7 dropped.

0.18.0 (January 2022)
---------------------

Minor new features:

* Add ``NodeGraph`` as abstract base class for device connectivity graphs.
* Improved ``CnX`` gate decomposition.
* Squashing of adjacent ``PhasedX`` operations.
* Add pytket ``__version__`` attribute.

Fixes:

* Fix wire-swap handling in ``PhasePolyBox`` creation.

0.17.0 (November 2021)
----------------------

Major new features:

* New ``pytket.zx`` module for manipulating ZX diagrams.

Minor new features:

* New properties: :py:meth:``circuit.Op.dagger`` and :py:meth:``circuit.Op.transpose``.
* New methods: :py:meth:``routing.Placement.to_dict`` and :py:meth:``routing.Placement.from_dict``.
* New ``NPhasedX`` OpType.
* New ``GlobalPhasedXPredicate`` and ``GlobalisePhasedX`` (transform and pass).

Fixes:

* Fixed incorrect decomposition of ``QControlBox`` with more than one control
  acting on operation with global phase.

0.16.0 (October 2021)
---------------------

Minor new features:

* New :py:meth:``backends.Backend.run_circuit`` and
  :py:meth:``backends.Backend.run_circuits`` methods.
* New ``allow_swaps`` parameter to ``FullPeepholeOptimise`` pass controlling
  whether to allow introduction of implicit wire swaps (defaulting to ``True``
  to match existing behaviour).
* New ``Backend.available_devices`` method to retrieve available devices as a
  list of ``BackendInfo`` objects.

Fixes:

* Fixed bug in daggering of TK1 gates.

API changes:

* The deprecated ``get_shots``, ``get_counts`` and ``get_state`` methods of the
  ``Backend`` class are removed. Use ``run_circuits`` and the homonym methods of
  the :py:class:`backends.backendresult.BackendResult` class instead.

0.15.0 (September 2021)
-----------------------

Minor new features:

* Passes ``PauliSimp``, ``PauliSquash`` and ``GuidedPauliSimp`` can now
  decompose to three-qubit ``XXPhase3`` gates using the new
  ``CXConfigType.MultiQGate`` config type.
* New method ``compilation_pass_from_script`` to construct a compilation pass
  from a simple textual specification.
* New transform ``RebaseToTket`` and new pass ``SquashToTK1``.

API changes:

* The deprecated transform ``RebaseToQiskit`` and the deprecated passes
  ``DecomposeMultiQubitsIBM``, ``RebaseIBM``, ``SynthesiseIBM`` and
  ``USquashIBM`` are removed.
* The transform ``OptimisePostRouting`` transforms to TK1 instead of U gates.

0.14.0 (September 2021)
-----------------------

Major new features:

* New ``Circuit.add_assertion`` method for applying quantum assertions to circuits.
* Two new box types  ``StabiliserAssertionBox`` and ``ProjectorAssertionBox``.
* New ``BackendResult.get_debug_info`` method for summarising assertion results.
* New ``PauliStabiliser`` class.
* Native support for MacOS running on M1 (arm64) architecture (Python 3.8 and 3.9 only).
* New compilerpass for architecture aware synthesis of phase polynomials ``AASRouting``.

Minor new features:

* Update circuit display to include extra gate information and use ZX-style colours.
* `BackendInfo`, `Architecture` and `Node` are now JSON-serializable.
* `QubitPauliOperator` and `QubitPauliString` are now JSON-serializable.
* Equality checks on `Architecture` only consider node IDs and coupling.
* New pass `DecomposeMultiQubitsCX`, equivalent to `DecomposeMultiQubitsIBM` (which is deprecated).
* New pass `DecomposeSingleQubitsTK1`.
* New pass `SynthesiseTket`.
* New ``XXPhase3`` OpType.

API changes:

* The transforms `ReduceSingles`, `OptimisePauliGadgets` and `OptimisePhaseGadgets`, and the passes `CliffordSimp`, `PeepholeOptimise2Q`, `FullPeepholeOptimise` and `OptimisePhaseGadgets`, produce TK1 instead of U gates.
* The passes `O2Pass`, `O1Pass` and `DecomposeSingleQubitsIBM` are removed (use `FullPeepholeOptimise` and `SynthesiseTket` instead for the first two).
* `QubitPauliOperator.to_dict()` (deprecated) is replaced by the property `QubitPauliOperator.map`.

Deprecations:

* The passes`DecomposeMultiQubitsIBM` (equivalent to `DecomposeMultiQubitsCX`), `DecomposeSingleQubitsIBM`, `RebaseToQiskit`, `SynthesiseIBM`, `RebaseIBM` and `USquashIBM` are deprecated.


0.13.0 (July 2021)
------------------

Major new features:

* New circuit functions, e.g. ``get_unitary``, calculate numerical unitaries and statevectors from non-symbolic circuits.
* New serialization methods for compilation passes.

Minor new features:

* Additions to `BackendInfo`.
* More reliable handling of timeouts for placement.
* User-configurable placement timeout.

Fixes:

* Fixed occasional segfault in placement pass.
* Daggering or transposing circuits with CnX fixed to have valid operation arguments.

API changes:

* :py:meth:`Backend.compile_circuit` is deprecated,
  :py:meth:`Backend.get_compiled_circuit` and
  :py:meth:`Backend.get_compiled_circuits` (for a sequence of circuits) replace
  it, do not act in place, returning the compiled circuit(s). In place
  compilation can still be achieved with `backend.default_compilation_pass().apply(circ)`

0.12.0 (June 2021)
------------------

Major new features:

* New ``ThreeQubitSquash`` compilation pass to simplify long three-qubit subcircuits.
* Three-qubit squash included in ``FullPeepholeOptimise`` pass; new ``PeepholeOptimise2Q`` pass corresponds to former ``FullPeepholeOptimise``.

Minor new features:

* add_phase now returns the circuit
* Option for `process_circuits` to take a list of `n_shots`.
* `Device` class removed, replaced with :py:class:`BackendInfo`.
* ``QubitErrorContainer`` removed.
* ``RoutingMethod`` removed.

Bugfixes and improvements:

* Barriers no longer count towards circuit depth.
* Squashing of rotations with symbolic angles now performs more simplification, leading to much shorter expressions, and works around a bug in symengine that caused invalid simplification of some expressions.

0.11.0 (May 2021)
-----------------
Major new features:

* New ``pytket.utils.symbolic`` module to generate symbolic unitaries and statevectors from symbolic circuits.
* New box type ``Unitary3qBox`` implementing arbitrary 3-qubit unitaries.

Minor new features:

* New ``ECR`` OpType.
* New ``SynthesiseOQC`` pass.
* New ``RebaseOQC`` pass.
  
0.10.1 (May 2021)
-----------------

Minor new features:

* New ``PauliSquash`` pass combining ``PauliSimp`` with ``FullPeepholeOptimise``.
* New options for ``SimplifyInitial``.

0.10.0 (April 2021)
-------------------

Major new features:

* HTML rendering of Circuit in Jupyter notebooks, ``pytket.circuit.display.render_circuit_jupyter``.

Minor new features:

* EulerAngleReduction pass uses multi-qubit commutativity to reduce rotation triplets to pairs
* EulerAngleReduction takes additional strictness parameter
* RemoveBarriers pass added.

API changes:

* Remove architecture classes :py:class:`TriangularGrid`, :py:class:`HexagonalGrid` and :py:class:`CyclicButterfly`

Fixes:

* Several small bugfixes.

0.9.0 (March 2021)
------------------

Major new features:

* Contextual optimizations based on knowledge of state.

Minor new features:

* New box type ``PhasePolyBox``.
* Refactored PytketConfig. `pytket-qiskit`, `pytket-honeywell`, `pytket-aqt`, `pytket-ionq`, `pytket-qsharp` and `pytket-braket`
  now all have authentication or workspace parameters that can be set in config files.

Fixes:

* Several small bugfixes.

0.8.0 (March 2021)
------------------

API changes:

* All extension modules moved to `pytket.extensions` namespace.

Compatible extension versions:

* ``pytket-aqt``: 0.5.0
* ``pytket-braket``: 0.4.0
* ``pytket-cirq``: 0.8.0
* ``pytket-honeywell``: 0.7.0
* ``pytket-ionq``: 0.3.0
* ``pytket-projectq``: 0.7.0
* ``pytket-pyquil``: 0.8.0
* ``pytket-pyzx``: 0.7.0
* ``pytket-qiskit``: 0.8.0
* ``pytket-qsharp``: 0.9.0
* ``pytket-qulacs``: 0.5.0

0.7.2 (February 2021)
---------------------

Major new features:

* Support for Python 3.9, dropping 3.6.

Fixes:

* Fix memory corruption with symbolic circuits on Windows.

0.7.1 (February 2021)
--------------------------

Minor new features:

* Option to store encrypted Honeywell password (not recommended).
* Automatic retries for Honeywell result retrieval.

Fixes:

* Drop dependency on OpenFermion (conversions work with separate installation).
* Fix reset breaking ``AerBackend`` ``_process_model``.
* Fix ``IBMQEmulatorBackend`` not being initialised with noise model.


Compatible extension versions:

* ``pytket-aqt``: 0.4.0
* ``pytket-braket``: 0.3.0
* ``pytket-cirq``: 0.7.0
* ``pytket-honeywell``: 0.6.1
* ``pytket-ionq``: 0.2.0
* ``pytket-projectq``: 0.6.0
* ``pytket-pyquil``: 0.7.0
* ``pytket-pyzx``: 0.6.0
* ``pytket-qiskit``: 0.7.1
* ``pytket-qsharp``: 0.8.2
* ``pytket-qulacs``: 0.4.0


0.7.0 (February 2021)
--------------------------

Major new features:

* Subsitution of named operations with other operations, boxes or circuits.
* New ability to condition operations on compound (AND, OR, XOR) operations on ``Bit`` and ``BitRegister``,
  which can be compiled with ``DecomposeClassicalExp`` and executed with ``HoneywellBackend``.

Minor new features:

* Direct creation of operator from gate type and parameters (``Op.create``).
* New methods ``Circuit.ops_of_type`` and ``Circuit.commands_of_type``.
* ``KAKDecomposition`` now accepts the estimated CX gate fidelity as parameter
  and performs an approximate decomposition in that case.
* Significant optimisation of SPAM correction methods.
* New GraphColourMethod.Exhaustive added to gen_term_sequence_circuit
  for partitioning Pauli tensors.
* New OpTypes ``CRx`` and ``CRy``.
* New OpTypes ``SX``, ``SXdg``, ``CSX``, ``CSXdg``, ``CV`` and ``CVdg``.
* New ``BasePass.get_config()`` method, which returns the name and parameters
  for a pass.
* New ``SequencePass.get_sequence()`` method, which returns the sequence of passes.
* New ``get_pass()`` method for ``RepeatPass``, ``RepeatWithMetricPass``, ``RepeatUntilSatisfiedPass``.
* New ``get_predicate()`` method for ``RepeatUntilSatisfiedPass``.
* New ``get_metric()`` method for ``RepeatWithMetricPass``.
* New ``backend`` parameter to ``SpamCorrecter`` constructor.

New supported backends:

* Support for Azure Quantum backends in the ``pytket-qsharp`` extension.

New features in extensions:

* Conversion of ``Reset`` and custom gates in ``pytket-qiskit``.
* Support for mid-circuit measurements on IBMQ premium devices via ``pytket-qiskit``.

API changes:

* Removal of "minimise" method for SPAM correction

Compatible extension versions:

* ``pytket-aqt``: 0.4.0
* ``pytket-braket``: 0.3.0
* ``pytket-cirq``: 0.7.0
* ``pytket-honeywell``: 0.6.0
* ``pytket-ionq``: 0.2.0
* ``pytket-projectq``: 0.6.0
* ``pytket-pyquil``: 0.7.0
* ``pytket-pyzx``: 0.6.0
* ``pytket-qiskit``: 0.7.0
* ``pytket-qsharp``: 0.8.0
* ``pytket-qulacs``: 0.4.0

0.6.1 (October 2020)
--------------------

Minor New Features:

* New pass generator ``RenameQubitsPass``

New Supported Backends:

* Devices from IonQ (via separate ``pytket-ionq`` module)

0.6.0 (September 2020)
----------------------

Major New Features:

* Windows support
* Phase-aware circuits
* New box type for applying quantum controls to arbitrary quantum operations
* New ``tailoring`` module containing tools for noise tailoring
* Circuit transpose method
* Optimization levels for default backend compilation passes
* New serialization methods for circuits and results
* New online user manual

Minor New Features:

* New gate type ``OpType.PhasedISWAP``
* Expectations of non-Hermitian operators (when supported by backend)
* Greater control over graph-colouring algorithms
* Improved Clifford simplification
* Retrieval of gate set from ``GateSetPredicate``
* New ``Backend.cancel`` method
* New ``name`` attribute for circuits.
* Backends can be wrapped as Qiskit backends for use in Qiskit software.
* IBMQEmulatorBackend added to emulate IBMQBackend behaviour, with simulator execution.

New supported backends:

* Devices and simulators from Amazon Braket (via separate ``pytket-braket``
  module)
* Qulacs simulator (via separate ``pytket-qulacs`` module)

.. * IonQ devices (via separate ``pytket-ionq`` module)

API changes:

* Retrieval of shots, counts, state and unitary directly from ``ResultHandle``
  is no longer supported: either use ``Backend.get_shots(Circuit)`` or
  ``Backend.get_result(ResultHandle).get_shots()`` (etc).
* ``Backend.default_compilation_pass`` is no longer a property but a method.
* ``QubitMap`` is replaced by a Python dictionary.
* Bit ordering of `condition_value` for conditionals now follows QASM convention
  (opposite to before, now `[0, 1]` corresponds to value 2).

Bugfixes:

* Various small bug fixes

Known issues:

* There is an `issue <https://github.com/CQCL/pytket/issues/24>`_ with the use
  of symbolic circuits on Windows, causing memory access violations in some
  circumstances.

Compatible extension versions:

* ``pytket-aqt``: 0.3.0
* ``pytket-braket``: 0.2.0
* ``pytket-cirq``: 0.5.0
* ``pytket-honeywell``: 0.4.0
* ``pytket-projectq``: 0.5.0
* ``pytket-pyquil``: 0.6.0
* ``pytket-pyzx``: 0.5.0
* ``pytket-qiskit``: 0.6.0
* ``pytket-qsharp``: 0.6.0
* ``pytket-qulacs``: 0.3.0

.. * ``pytket-ionq``: 0.1.0

0.5.7 (August 2020)
-------------------
Number of bugs fixed including:


* ``OpType.Reset`` added to QASM conversion
* Bugfix for ``CnX`` with n=4, n=5
* Correct Node IDS for ``FullyConnected`` Architecture.


0.5.5 (June 2020)
-----------------
Major New Features:

* Redesigned algorithm for ``CliffordSimp``, improving speed and identifying more cases for optimisation

Minor New Features:

* New gates added: ``OpType.Sycamore`` and ``OpType.ISWAPMax``
* New class ``Graph`` for visualising circuit structure

Updates:

* First parameter of ``OpType.FSim`` gate corrected to have range :math:`[0, 2\pi)`
* New ``QubitPauliOperator`` and related classes replace use of OpenFermion's ``QubitOperator``
* Significant optimisation of ``pauli_tensor_matrix`` and ``operator_matrix``


0.5.4 (May 2020)
------------------
Minor New Features:

* Method to generate a circuit from a sequence of ``QubitOperator`` terms

Updates:

* Rename ``measurement`` module to ``partition``

Bugfixes:

* Fix invalid cancellation of certain controlled rotations


0.5.2 (April 2020)
------------------
Major New Features:

* Routing, gate decomposition, and basic optimisations can work around conditional gates and mid-circuit measurements
* New high-level optimisation routine for Trotterised Hamiltonians
* Measurement reduction via Pauli term diagonalisation
* Inspection of the status of circuit execution on asynchronous backends
* Error mitigation facilities via the SPAM method
* Introduction of the :py:class:`Program` class for specifying routines with classical control flow

Minor New Features:

* Improved error messages when circuits cannot be run on a backend
* Generalised :py:meth:`Circuit.depth_by_type` to allow sets of gate types
* A selection of optimisation passes are parameterised by pattern for decomposing into CXs
* New :py:class:`Architecture` subclass, :py:class:`FullyConnected`, added
* New gates added: `OpType.ESWAP` and `OpType.FSim`
* Additional utility methods for permuting qubits of statevectors
* Inspection of any implicit permutations within the :py:class:`Circuit` dag structure
* Inspection of free symbols in a circuit
* Inspection of detailed gate errors from a :py:class:`Device`
* Additional methods for parsing/producing QASM through strings and streams
* Ability to enable internal logs

Updates:

* Cleaner addition of conditions to gates via kwargs
* :py:class:`UnitID` objects are specialised into either :py:class:`Qubit` or :py:class:`Bit` objects, with more natural constructors
* Renamed many passes to give a uniform naming convention
* Getters on :py:class:`Architecture`, :py:class:`Device`, :py:class:`GateError`, and :py:class:`QubitErrorContainer` made into readonly properties
* Backend-specific runtime arguments (e.g. simulator seeds) are now passed in via kwargs
* Stability improvements and bug fixes
* Updated documentation and additional examples
* Stricter namespacing (most classes must be imported from submodules rather than top level)
* Python 3.8 support

Deprecations:

* Calling :py:meth:`get_counts`, :py:meth:`get_shots` or :py:meth:`get_state` on a :py:class:`Backend` object with a :py:class:`Circuit` argument is deprecated in favour of :py:class:`ResultHandle`.

New supported backends:

* AQT devices and simulators (via separate ``pytket_aqt`` module)
* Honeywell devices (via separate ``pytket_honeywell`` module)
* Q# simulators and resource estimator (via separate ``pytket_qsharp`` module)

0.4.1 (December 2019)
---------------------
New Features:

* New classes for placement of logical qubits from :py:class:`Circuit` to physical qubits from :py:class:`Device` or :py:class:`Architecture`
* Data from backends can be returned in either increasing lexicographical order of (qu)bit identifiers (the familiar ordering used in most textbooks) or decreasing order (popular with other quantum software platforms) using the :py:class:`BasisOrder` enum

Updates:

* Updated documentation and additional examples
* OptimiseCliffordsZX pass removed, FullPeepholeOptimise pass added
* New architectures added, including :py:class:`SquareGrid`, :py:class:`HexagonalGrid`, :py:class:`RingArch`, :py:class:`TriangularGrid` and :py:class:`CyclicButterfly`
* Device information from :py:class:`Device` can now be returned
* Stability improvements and bug fixes

0.4.0 (November 2019)
---------------------
New Features:

* Contractural compilation passes with guarantees on how they transform circuits that satisfy their preconditions. This provides a uniform interface for optimisations, routing, and other stages of compilation
* New "Box" gate types for encapsulating high-level structures (arbitrary subcircuits, parameterised composite gate definitions, unitaries, Pauli operators)
* Simpler and more flexible structure for registers and names of qubits/bits, allowing for non-contiguous and multi-dimensional indices (referring to individual units, linear registers, grids, etc.)
* Latex diagram output using Quantikz
* The :py:class:`Device` class to build on top of :py:class:`Architecture` with error and timing information
* Initial and final maps tracked throughout the entire compilation procedure using the :py:class:`CompilationUnit` wrapper
* Import circuits from Quipper source files
* Utility methods for processing data from Backends

Updates:

* All Backends refactored for more consistent interfaces, separation of data processing, and introducing batch circuit processing when possible
* Routing improved to use distributed CX (BRIDGE) gates in addition to SWAP insertion
* Cost function for noise-aware allocation of qubits improved to consider more sources of noise
* :py:class:`Architecture` objects can be specified with arbitrary node names, using the same :py:class:`UnitID` objects and qubits/bits
* Removed the :py:class:`PhysicalCircuit` class in preference of just using :py:class:`Circuit` objects
* Generalised and sped up the gate commutation pass
* Optimisation for redundant gate removal now removes diagonal gates before measurements
* Support for custom gate definitions in QASM input
* Support for a greater fragment of sympy expressions in gate parameters
* Stability improvements and bug fixes
* Updated documentation and additional examples

0.3.0 (August 2019)
-------------------
New Features:

* More options for circuit routing, including noise-aware allocation of qubits
* Basic support for generating circuits with classical conditions and multiple registers
* ForestBackend for running circuits on Rigett's QVM simulators and QCS
* AerUnitaryBackend for inspecting the full unitary of a circuit
* Chaining gate commands
* Primitive QASM<->Circuit (import and export)

Updates:

* Simplified conversions for pytket_qiskit, going straight to/from QuantumCircuit rather than DAGCircuit
* CSWAP gate added

0.2.3 (July 2019)
------------------
New Features:

* Decomposition `Transform` for controlled gates

Updates:

* Exposed additional gate types into Pytket
* Fixed bug in `add_circuit`
* Fixed routing bug
* Made `run` behaviour more sensible for backends

0.2.2 (June 2019)
------------------
Updates:

* Minor bug fixes, examples and documentation

0.2.1 (June 2019)
------------------
Updates:

* Extra support for appending Circuits from Matrices and Exponents
* More docs and examples
* Fixed bugs in backends

0.2.0 (June 2019)
------------------
New Features:

* Support for circuits and simulation using ProjectQ (0.4.2)
* Support for conversion to and from PyZX (https://github.com/Quantomatic/pyzx)
* Interface to many new optimisation passes, allowing for custom passes
* Circuit compilation using symbolic parameters
* New interface to routing
* Enabled noise modeling in the AerBackend module

Updates:

* Qiskit support updated for Qiskit 0.10.1 and Qiskit Chemistry 0.5
* Pytket Chemistry module has been removed, to be part of the separate Eumen package
* Bug fixes and performance improvements to routing

0.1.6 (April 2019)
------------------
Updates:

* Routing can return SWAP gates rather than decomposing to CNOTs
* Decomposition and routing bug fixes

0.1.5 (April 2019)
------------------
New Features:

* Enabled conversions from 4x4 unitary matrices to 2 qubit circuit

0.1.4 (April 2019)
------------------
Updates:

* Bug fix patch for routing and performance improvements

0.1.3 (March 2019)
------------------
Updates:

* Qiskit support updated for Terra 0.7.3, Aqua 0.4.1, and Chemistry 0.4.2
* Bug fixes in routing

0.1.2 (February 2019)
---------------------
New Features:

* Support for circuits from Rigetti pyQuil (2.3)
* New interface for constructing and analysing circuits in pytket directly
* Named classical registers for measurements

Updates:

* Documentation and tutorial improvements
* Bug fixes in routing and optimisations
* Minor API changes for notational consistency

0.1.0 (December 2018)
---------------------
New Features:

* Support for circuits and architectures from IBM Qiskit (0.7)
* ``pytket.qiskit.TketPass`` allows pytket to be plugged in to the Qiskit compilation stack to take advantage of tket's routing and optimisations
* New Chemistry package featuring an implementation of the Quantum Subspace Expansion to work within or alongside Qiskit Aqua (0.4)
* Optimisation passes introduced for powerful circuit rewriting before routing, and safe rewriting after routing

Updates:

* Cirq functionality supports Cirq 0.4
* Refactoring into modules

0.0.1 (July 2018)
-----------------
New Features:

* Support for circuits and architectures from Google Cirq (0.3)
* Routing and placement procedures available for manipulating circuits to satisfy device specifications
