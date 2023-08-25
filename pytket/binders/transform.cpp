// Copyright 2019-2023 Cambridge Quantum Computing
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "tket/Transformations/Transform.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <functional>

#include "tket/Circuit/Circuit.hpp"
#include "tket/Transformations/BasicOptimisation.hpp"
#include "tket/Transformations/Combinator.hpp"
#include "tket/Transformations/ContextualReduction.hpp"
#include "tket/Transformations/Decomposition.hpp"
#include "tket/Transformations/OptimisationPass.hpp"
#include "tket/Transformations/PauliOptimisation.hpp"
#include "tket/Transformations/Rebase.hpp"
#include "tket/Transformations/ThreeQubitSquash.hpp"
#include "typecast.hpp"

namespace py = pybind11;

namespace tket {

// given keyword arguments for DecomposeTK2, return a TwoQbFidelities struct
Transforms::TwoQbFidelities get_fidelities(const py::kwargs &kwargs) {
  Transforms::TwoQbFidelities fid;
  for (const auto &kwarg : kwargs) {
    const std::string kwargstr = py::cast<std::string>(kwarg.first);
    using Func = std::function<double(double)>;
    if (kwargstr == "CX_fidelity") {
      fid.CX_fidelity = py::cast<double>(kwarg.second);
    } else if (kwargstr == "ZZMax_fidelity") {
      fid.ZZMax_fidelity = py::cast<double>(kwarg.second);
    } else if (kwargstr == "ZZPhase_fidelity") {
      fid.ZZPhase_fidelity = py::cast<Func>(kwarg.second);
    } else {
      throw py::type_error(
          "got an unexpected keyword argument '" + kwargstr + "'");
    }
  }
  return fid;
}

PYBIND11_MODULE(transform, m) {
  // py::module::import("pytket._tket.circuit");
  // py::module::import("pytket._tket.architecture");
  py::enum_<Transforms::PauliSynthStrat>(
      m, "PauliSynthStrat",
      "Enum for available strategies to synthesise Pauli gadgets")
      .value(
          "Individual", Transforms::PauliSynthStrat::Individual,
          "Synthesise gadgets individually")
      .value(
          "Pairwise", Transforms::PauliSynthStrat::Pairwise,
          "Synthesise gadgets using an efficient pairwise strategy "
          "from Cowtan et al (https://arxiv.org/abs/1906.01734)")
      .value(
          "Sets", Transforms::PauliSynthStrat::Sets,
          "Synthesise gadgets in commuting sets");

  py::class_<Transform>(
      m, "Transform", "An in-place transformation of a :py:class:`Circuit`.")
      .def(py::init<const Transform::SimpleTransformation &>())
      .def(
          "apply",
          [](const Transform &tr, Circuit &circ) { return tr.apply(circ); },
          "Performs the transformation on the circuit in "
          "place.\n\n:param circuit: The circuit to be "
          "transformed\n"
          ":return: True if any changes were made, else False",
          py::arg("circuit"))

      /* COMBINATORS */
      .def(
          py::self >> py::self,
          "Composes two Transforms together in sequence.\n\n>>> a >> "
          "b\n\nis equivalent to\n\n>>> sequence([a,b])")
      .def_static(
          "sequence", &Transforms::sequence,
          "Composes a list of Transforms together in sequence. The "
          ":py:meth:`apply` method will return ``True`` if ANY of "
          "the individual Transforms returned ``True``."
          "\n\n:param sequence: The list of Transforms to be "
          "composed\n:return: the combined Transform",
          py::arg("sequence"))
      .def_static(
          "repeat", &Transforms::repeat,
          "Applies a given Transform repeatedly to a circuit until "
          "no further changes are made (i.e. it no longer returns "
          "``True``). "
          ":py:meth:`apply` will return ``True`` if at least one run "
          "returned ``True``."
          "\n\n:param transform: The Transform to be applied "
          "repeatedly\n:return: a new Transform representing the "
          "iteration",
          py::arg("transform"))
      .def_static(
          "while_repeat", &Transforms::repeat_while,
          "Repeatedly applies the `condition` Transform until it "
          "returns ``False``, running `body` in between each "
          "`condition` application. "
          "Intuitively, this corresponds to \"WHILE `condition` DO "
          "`body`\"."
          "\n\n:param condition: The Transform to be applied "
          "repeatedly as the condition of a loop"
          "\n:param body: The Transform to be applied after each "
          "successful test of the condition"
          "\n:return: a new Transform representing the iteration",
          py::arg("condition"), py::arg("body"))

      /* REBASE TRANSFORMS */
      .def_static(
          "RebaseToTket", &Transforms::rebase_tket,
          "Rebase from any gate set into TK1, CX.")
      .def_static(
          "RebaseToRzRx", &Transforms::decompose_ZX,
          "Rebase single qubit gates into Rz, Rx.")
      .def_static(
          "RebaseToCliffordSingles", &Transforms::decompose_cliffords_std,
          "Replace all single-qubit unitary gates outside the set {Z, X, S, V} "
          "that are recognized as Clifford operations with an equivalent "
          "sequence of gates from that set.")
      .def_static(
          "RebaseToCirq", &Transforms::rebase_cirq,
          "Rebase from any gate set into PhasedX, Rz, CZ.")
      .def_static(
          "RebaseToQuil", &Transforms::rebase_quil,
          "Rebase from any gate set into Rx, Rz, CZ.")
      .def_static(
          "RebaseToPyZX", &Transforms::rebase_pyzx,
          "Rebase from any gate set into the gate set supported by "
          "PyZX (Rx, Rz, X, Z, S, T, H, CX, CZ, SWAP).")
      .def_static(
          "RebaseToProjectQ", &Transforms::rebase_projectq,
          "Rebase from any gate set into the gate set supported by "
          "ProjectQ (Rx, Ry, Rz, X, Y, Z, S, T, V, H, CX, CZ, CRz, "
          "SWAP).")
      .def_static(
          "DecomposeCCX", &Transforms::decomp_CCX,
          "Decomposes all 3-qubit Toffoli (CCX) gates into "
          "Clifford+T gates.")
      .def_static(
          "DecomposeControlledRys", &Transforms::decomp_controlled_Rys,
          "Decomposes all arbitrarily-quantum-controlled Rys into CX "
          "and Ry gates.")
      .def_static(
          "DecomposeSWAP", &Transforms::decompose_SWAP,
          "Decomposes all SWAP gates to provided replacement "
          "circuit.\n\n:param circuit: A circuit that is logically "
          "equivalent to a SWAP operation",
          py::arg("circuit"))
      .def_static(
          "DecomposeSWAPtoCX", &Transforms::decompose_SWAP_to_CX,
          "Decomposes all SWAP gates into triples of CX gates. "
          "If the SWAP is adjacent to a CX, it will prefer to insert "
          "in the direction that allows for gate cancellation. "
          "If an :py:class:`Architecture` is provided, this will "
          "prefer to insert the CXs such that fewer need redirecting."
          "\n\n:param arc: Device architecture used to specify a "
          "preference for CX direction",
          py::arg("arc"))
      .def_static(
          "DecomposeBRIDGE", &Transforms::decompose_BRIDGE_to_CX,
          "Decomposes all BRIDGE gates into CX gates.")
      .def_static(
          "DecomposeCXDirected", &Transforms::decompose_CX_directed,
          "Decompose CX gates to H+CX to match the direction of the "
          "CXs to edges of the :py:class:`Architecture` `arc`. "
          "Assumes the circuit already satisfies the connectivity of "
          "`arc`.\n\n:param arc: The architecture for which CXs "
          "should be redirected",
          py::arg("arc"))
      .def_static(
          "DecomposeBoxes", &Transforms::decomp_boxes,
          "Decomposes all Boxed operations into elementary gates.")
      .def_static(
          "DecomposeTK2",
          [](bool allow_swaps, const py::kwargs &kwargs) {
            return Transforms::decompose_TK2(
                get_fidelities(kwargs), allow_swaps);
          },
          "Decompose each TK2 gate into two-qubit gates."
          "\n\nWe currently support CX, ZZMax and ZZPhase."
          "\n\nIf one or more gate fidelities are provided, the two-qubit gate "
          "type achieving the highest fidelity will be chosen for the "
          "decomposition, as measured using squared trace fidelity. "
          "If no fidelities are provided, the TK2 gates will be decomposed "
          "exactly using CX gates. For equal fidelities, ZZPhase will be "
          "prefered over ZZMax and CX if the decomposition results in fewer "
          "two-qubit gates.\n\n"
          "All TK2 gate parameters must be normalised, i.e. they must satisfy "
          "`NormalisedTK2Predicate`."
          "\n\n"
          "Gate fidelities are passed as keyword arguments to perform "
          "noise-aware decompositions. "
          "We currently support `CX_fidelity`, `ZZMax_fidelity` and "
          "`ZZPhase_fidelity`. If provided, the `CX` and `ZZMax` fidelities "
          "must be given by a single floating point fidelity. The `ZZPhase` "
          "fidelity is given as a lambda float -> float, mapping a ZZPhase "
          "angle parameter to its fidelity. These parameters will be used "
          "to return the optimal decomposition of each TK2 gate, taking "
          "noise into consideration.\n\n"
          "Using the `allow_swaps=True` (default) option, qubits will be "
          "swapped when convenient to reduce the two-qubit gate count of the "
          "decomposed TK2.\n\n"
          "If the TK2 angles are symbolic values, the decomposition will "
          "be exact (i.e. not noise-aware). It is not possible in general "
          "to obtain optimal decompositions for arbitrary symbolic parameters, "
          "so consider substituting for concrete values if possible."
          "\n\n:param allow_swaps: Whether to allow implicit wire swaps.",
          py::arg("allow_swaps") = true)
      .def_static(
          "NormaliseTK2", &Transforms::normalise_TK2,
          "Normalises all TK2 gates.\n\n"
          "TK2 gates have three angles in the interval [0, 4], but these can "
          "always be normalised to be within the so-called Weyl chamber by "
          "adding single-qubit gates.\n\n"
          "More precisely, the three angles a, b, c of TK2(a, b, c) are "
          "normalised exactly when the two following conditions are met:\n"
          " - numerical values must be in the Weyl chamber, "
          "ie `1/2 >= a >= b >= |c|`,\n"
          " - symbolic values must come before any numerical value in the "
          "array.\n\n"
          "After this transform, all TK2 angles will be normalised and the "
          "circuit will satisfy `NormalisedTK2Predicate`.")

      /* OPTIMISATION TRANSFORMS */
      .def_static(
          "OptimiseStandard", &Transforms::synthesise_tk,
          "Fast optimisation pass, performing basic simplifications. "
          "Works on any circuit, giving the result in TK1 and TK2 gates. "
          "Preserves connectivity of circuit.")
      .def_static(
          "OptimisePostRouting", &Transforms::synthesise_tket,
          "Fast optimisation pass, performing basic simplifications. "
          "Works on any circuit, giving the result in TK1 and CX gates. "
          "If all multi-qubit gates are CXs, then this preserves "
          "their placement and orientation, so it is safe to perform "
          "after routing.")
      .def_static(
          "OptimisePhaseGadgets", &Transforms::optimise_via_PhaseGadget,
          "An optimisation pass that starts by identifying "
          "subcircuits corresponding to phase gadgets (see Cowtan, "
          "Duncan, Dilkes, Simmons, & Sivarajah "
          "https://arxiv.org/abs/1906.01734) and resynthesises them "
          "in a balanced-tree form, followed by applying "
          "OptimisePostRouting. "
          "Results use TK1 and CX gates. This will not preserve "
          "CX placement or orientation.",
          py::arg("cx_config") = CXConfigType::Snake)
      .def_static(
          "OptimiseCliffords", &Transforms::clifford_simp,
          "An optimisation pass that performs a number of rewrite "
          "rules for simplifying Clifford gate sequences, similar to "
          "Duncan & Fagan (https://arxiv.org/abs/1901.10114). "
          "Given a circuit with CXs and any single-qubit gates, "
          "produces a circuit with Z, X, S, V, U1, U2, U3, CX gates. "
          "This will not preserve CX placement or orientation and "
          "may introduce implicit wire swaps."
          "\n\n:param allow_swaps: dictates whether the rewriting "
          "will disregard CX placement or orientation and introduce "
          "wire swaps.",
          py::arg("allow_swaps") = true)
      .def_static(
          "OptimisePauliGadgets", &Transforms::pairwise_pauli_gadgets,
          "An optimisation pass that identifies the Pauli gadgets "
          "corresponding to any non-Clifford rotations and "
          "synthesises them pairwise (see Cowtan, Duncan, Dilkes, "
          "Simmons, & Sivarajah https://arxiv.org/abs/1906.01734). "
          "Results use TK1, CX gates.",
          py::arg("cx_config") = CXConfigType::Snake)
      .def_static(
          "RemoveRedundancies", &Transforms::remove_redundancies,
          "Applies a collection of simple optimisations, such as "
          "removing gate-inverse pairs, merging similar rotation "
          "gates, and removing identity gates. "
          "Preserves the gate set and any placement/orientation of "
          "multi-qubit gates.")
      .def_static(
          "ReduceSingles", &Transforms::squash_1qb_to_tk1,
          "Reduces each sequence of single-qubit rotations into a single TK1.")
      .def_static(
          "CommuteThroughMultis", &Transforms::commute_through_multis,
          "Applies a collection of commutation rules to move single "
          "qubit operations past multiqubit operations they commute "
          "with, towards the front of the circuit.")
      .def_static(
          "KAKDecomposition",
          py::overload_cast<OpType, double, bool>(
              &Transforms::two_qubit_squash),
          "Squash sequences of two-qubit operations into minimal form.\n\n"
          "Squash together sequences of single- and two-qubit gates "
          "into minimal form. Can decompose to TK2 or CX gates.\n\n"
          "Two-qubit operations can always be expressed in a minimal form of "
          "maximum three CXs, or as a single TK2 gate (a result also known "
          "as the KAK or Cartan decomposition).\n\n"
          "It is in general recommended to squash to TK2 gates, and to then use"
          " the `DecomposeTK2` pass for noise-aware decompositions to other "
          "gatesets. For backward compatibility, decompositions to CX are also "
          "supported. In this case, `cx_fidelity` can be provided to perform "
          "approximate decompositions to CX gates.\n\n"
          "When decomposing to TK2 gates, any sequence of two or more two-qubit"
          " gates on the same set of qubits is replaced by a single TK2 gate. "
          "When decomposing to CX, the substitution is only performed if it "
          "results in a reduction of the number of CX gates, or if at least "
          "one of the two-qubit passes is not a CX.\n\n"
          "Using the `allow_swaps=True` (default) option, qubits will be "
          "swapped when convenient to further reduce the two-qubit gate count. "
          "(only applicable when decomposing to CX gates).\n\n"
          ":param target_2qb_gate: OpType to decompose to. Either TK2 or CX.\n"
          ":param cx_fidelity: Estimated CX gate fidelity, used when "
          "target_2qb_gate=CX.\n"
          ":param allow_swaps: Whether to allow implicit wire swaps.",
          py::arg("target_2qb_gate") = OpType::CX, py::arg("cx_fidelity") = 1.,
          py::arg("allow_swaps") = true)
      .def_static(
          "KAKDecomposition",
          [](double cx_fidelity) {
            return Transforms::two_qubit_squash(OpType::CX, cx_fidelity);
          },
          py::arg("cx_fidelity"))
      .def_static(
          "ThreeQubitSquash", &Transforms::three_qubit_squash,
          "Squash three-qubit subcircuits into subcircuits having fewer "
          "2-qubit gates of the target type, when possible. The supported "
          "target types are CX (default) and TK2.",
          py::arg("target_2qb_gate") = OpType::CX)
      .def_static(
          "CommuteSQThroughSWAP",
          [](const avg_node_errors_t &avg_node_errors) {
            return Transforms::commute_SQ_gates_through_SWAPS(avg_node_errors);
          },
          "Commutes single qubit gates through SWAP gates, leaving "
          "them on the physical qubit with best fidelity for given "
          "gate type. "
          "Assumes the circuit is already mapped onto the "
          "architecture."
          "\n\n:param avg_node_errors: a dict mapping Nodes to average "
          "single-qubit gate errors",
          py::arg("avg_node_errors"))
      .def_static(
          "CommuteSQThroughSWAP",
          [](const op_node_errors_t &op_node_errors) {
            return Transforms::commute_SQ_gates_through_SWAPS(op_node_errors);
          },
          "Commutes single qubit gates through SWAP gates, leaving "
          "them on the physical qubit with best fidelity for given "
          "gate type. "
          "Assumes the circuit is already mapped onto the "
          "architecture."
          "\n\n:param avg_node_errors: a dict of dicts, mapping Nodes to dicts "
          "of OpType to single-qubit gate error maps",
          py::arg("op_node_errors"))
      .def_static(
          "DecomposeNPhasedX", &Transforms::decompose_NPhasedX,
          "Decompose NPhasedX gates into single-qubit PhasedX gates.")
      .def_static(
          "GlobalisePhasedX", &Transforms::globalise_PhasedX,
          "Turns all PhasedX and NPhasedX gates into global gates\n\n"
          "Replaces any PhasedX gates with global NPhasedX gates. "
          "By default, this transform will squash all single-qubit gates "
          "to PhasedX and Rz gates before proceeding further. "
          "Existing non-global NPhasedX will not be preserved. "
          "This is the recommended setting for best "
          "performance. If squashing is disabled, each non-global PhasedX gate "
          "will be replaced with two global NPhasedX, but any other gates will "
          "be left untouched."
          "\n\n:param squash: Whether to squash the circuit in pre-processing "
          "(default: true)."
          "\n\nIf squash=true (default), the `GlobalisePhasedX` transform's "
          "`apply` method "
          "will always return true. "
          "For squash=false, `apply()` will return true if the circuit was "
          "changed and false otherwise.\n\n"
          "It is not recommended to use this transformation with symbolic "
          "expressions, as in certain cases a blow-up in symbolic expression "
          "sizes may occur.",
          py::arg("squash") = true)
      .def_static(
          "SynthesisePauliGraph", &Transforms::synthesise_pauli_graph,
          "Synthesises Pauli Graphs.",
          py::arg("synth_strat") = Transforms::PauliSynthStrat::Sets,
          py::arg("cx_config") = CXConfigType::Snake)
      .def_static(
          "UCCSynthesis", &Transforms::special_UCC_synthesis,
          "Synthesises UCC circuits in the form that Term Sequencing "
          "provides them.",
          py::arg("synth_strat") = Transforms::PauliSynthStrat::Sets,
          py::arg("cx_config") = CXConfigType::Snake)
      .def_static(
          "ZZPhaseToRz", &Transforms::ZZPhase_to_Rz,
          "Fixes all ZZPhase gate angles to [-1, 1) half turns.")
      .def_static(
          "CnXPairwiseDecomposition", &Transforms::cnx_pairwise_decomposition,
          "Decompose CnX gates to 2-qubit gates and single qubit gates. "
          "For every two CnX gates, reorder their control qubits to improve "
          "the chance of gate cancellation.")
      .def_static(
          "round_angles", &Transforms::round_angles,
          "Rounds angles to the nearest :math:`\\pi / 2^n`."
          "\n\n:param n: precision parameter, must be >= 0 and < 32",
          "\n\n:param only_zeros: if True, only round angles less than "
          ":math:`\\pi / 2^{n+1}` to zero, leave other angles alone (default "
          "False)",
          py::arg("n"), py::arg("only_zeros") = false);
  m.def(
      "separate_classical", &Transforms::separate_classical,
      "Separate the input circuit into a 'main' circuit and a classical "
      "'post-processing' circuit, which are equivalent to the original "
      "when composed."
      "\n\n:param circ: circuit to be separated",
      py::arg("circ"));
}

}  // namespace tket
