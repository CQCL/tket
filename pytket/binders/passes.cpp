// Copyright 2019-2022 Cambridge Quantum Computing
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

#include <pybind11/functional.h>

#include "ArchAwareSynth/SteinerForest.hpp"
#include "Mapping/LexiLabelling.hpp"
#include "Mapping/LexiRoute.hpp"
#include "Mapping/RoutingMethod.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Predicates/PassLibrary.hpp"
#include "Transformations/ContextualReduction.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/PauliOptimisation.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/Json.hpp"
#include "binder_json.hpp"
#include "typecast.hpp"

namespace py = pybind11;
using json = nlohmann::json;

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

static PassPtr gen_cx_mapping_pass_kwargs(
    const Architecture &arc, const PlacementPtr &placer, py::kwargs kwargs) {
  std::vector<RoutingMethodPtr> config = {
      std::make_shared<LexiLabellingMethod>(),
      std::make_shared<LexiRouteRoutingMethod>()};
  if (kwargs.contains("config")) {
    config = py::cast<std::vector<RoutingMethodPtr>>(kwargs["config"]);
  }
  bool directed_cx = false;
  if (kwargs.contains("directed_cx")) {
    directed_cx = py::cast<bool>(kwargs["directed_cx"]);
  }
  bool delay_measures = true;
  if (kwargs.contains("delay_measures")) {
    delay_measures = py::cast<bool>(kwargs["delay_measures"]);
  }
  return gen_cx_mapping_pass(arc, placer, config, directed_cx, delay_measures);
}

static PassPtr gen_default_routing_pass(const Architecture &arc) {
  return gen_routing_pass(
      arc, {std::make_shared<LexiLabellingMethod>(),
            std::make_shared<LexiRouteRoutingMethod>()});
}

static PassPtr gen_default_aas_routing_pass(
    const Architecture &arc, py::kwargs kwargs) {
  unsigned lookahead = 1;
  aas::CNotSynthType cnotsynthtype = aas::CNotSynthType::Rec;

  if (kwargs.contains("lookahead"))
    lookahead = py::cast<unsigned>(kwargs["lookahead"]);

  if (kwargs.contains("cnotsynthtype"))
    cnotsynthtype = py::cast<aas::CNotSynthType>(kwargs["cnotsynthtype"]);

  if (lookahead == 0) {
    throw std::invalid_argument(
        "[AAS]: invalid input, the lookahead must be > 0");
  }

  return gen_full_mapping_pass_phase_poly(arc, lookahead, cnotsynthtype);
}

static PassPtr gen_composephasepolybox_pass(py::kwargs kwargs) {
  unsigned min_size = 0;

  if (kwargs.contains("min_size")) {
    min_size = py::cast<unsigned>(kwargs["min_size"]);
  }

  return RebaseUFR() >> ComposePhasePolyBoxes(min_size);
}

const PassPtr &DecomposeClassicalExp() {
  // a special box decomposer for Circuits containing
  // ClassicalExpBox<py::object>
  static const PassPtr pp([]() {
    Transform t = Transform([](Circuit &circ) {
      py::module decomposer =
          py::module::import("pytket.circuit.decompose_classical");
      const py::tuple result = decomposer.attr("_decompose_expressions")(circ);
      const bool success = result[1].cast<bool>();
      if (success) {
        circ = result[0].cast<Circuit>();
      }
      return success;
    });
    PredicatePtrMap s_ps;
    /**
     * Preserves Max2QubitGatesPredicate since any box with >2 qubits is
     * already invalid.
     * Preserves ConnectivityPredicate (and DirectednessPredicate) since the
     * verification looks inside CircBoxes and any other boxes with >2
     * qubits are already invalid.
     * Most others are preserved since the predicates look within CircBoxes.
     *
     * Invalidates GateSetPredicate because it adds Classical OpTypes
     */
    PredicateClassGuarantees g_postcons = {
        {typeid(GateSetPredicate), Guarantee::Clear},
    };
    PostConditions postcon{s_ps, g_postcons, Guarantee::Preserve};
    json j;
    j["name"] = "DecomposeClassicalExp";
    return std::make_shared<StandardPass>(s_ps, t, postcon, j);
  }());
  return pp;
}

PYBIND11_MODULE(passes, m) {
  py::enum_<SafetyMode>(m, "SafetyMode")
      .value(
          "Audit", SafetyMode::Audit,
          "Checks which predicates a circuit satisfies after the "
          "application of each base pass")
      .value(
          "Default", SafetyMode::Default,
          "Only check that a circuit satisfies the preconditions of "
          "the overall pass at the start and the postconditions at "
          "the end")
      // .value("Off", SafetyMode::Off) // not currently supported
      .export_values();

  py::enum_<aas::CNotSynthType>(m, "CNotSynthType")
      .value(
          "SWAP", aas::CNotSynthType::SWAP,
          "swap-based algorithm for CNOT synthesis")
      .value(
          "HamPath", aas::CNotSynthType::HamPath,
          "Hamilton-path-based method for CNOT synthesis; this method will "
          "fail if there is no Hamilton path in the given architecture")
      .value(
          "Rec", aas::CNotSynthType::Rec,
          "recursive Steiner--Gauss method for CNOT synthesis")
      .export_values();

  /* Compiler passes */

  class PyBasePass : public BasePass {
   public:
    using BasePass::BasePass;

    /* Trampolines (need one for each virtual function */
    virtual bool apply(
        CompilationUnit &c_unit, SafetyMode safe_mode = SafetyMode::Default,
        const PassCallback &before_apply = trivial_callback,
        const PassCallback &after_apply = trivial_callback) const override {
      PYBIND11_OVERLOAD_PURE(
          bool,     /* Return type */
          BasePass, /* Parent class */
          apply,    /* Name of function in C++ (must match Python name) */
          c_unit, before_apply, after_apply, safe_mode /* Argument(s) */
      );
    }
    virtual std::string to_string() const override {
      PYBIND11_OVERLOAD_PURE(
          std::string, /* Return type */
          BasePass,    /* Parent class */
          to_string    /* Name of function in C++ (must match Python name) */
      );
    }
    virtual json get_config() const override {
      PYBIND11_OVERLOAD_PURE(
          json,      /* Return type */
          BasePass,  /* Parent class */
          get_config /* Name of function in C++ (must match Python name) */
      );
    }
  };

  py::class_<BasePass, PassPtr, PyBasePass>(
      m, "BasePass", "Base class for passes.")
      .def(
          "apply",
          [](const BasePass &pass, CompilationUnit &cu,
             SafetyMode safety_mode) { return pass.apply(cu, safety_mode); },
          "Apply to a :py:class:`CompilationUnit`.\n\n"
          ":return: True if pass modified the circuit, else False",
          py::arg("compilation_unit"),
          py::arg("safety_mode") = SafetyMode::Default)
      .def(
          "apply",
          [](const BasePass &pass, Circuit &circ) {
            CompilationUnit cu(circ);
            bool applied = pass.apply(cu);
            circ = cu.get_circ_ref();
            return applied;
          },
          "Apply to a :py:class:`Circuit` in-place.\n\n"
          ":return: True if pass modified the circuit, else False",
          py::arg("circuit"))
      .def(
          "apply",
          [](const BasePass &pass, Circuit &circ,
             const PassCallback &before_apply,
             const PassCallback &after_apply) {
            CompilationUnit cu(circ);
            bool applied =
                pass.apply(cu, SafetyMode::Default, before_apply, after_apply);
            circ = cu.get_circ_ref();
            return applied;
          },
          "Apply to a :py:class:`Circuit` in-place and invoke "
          "callbacks "
          "for all nested passes.\n\n"
          "\n:param before_apply: Invoked before a pass is applied. "
          "The CompilationUnit and a summary of the pass "
          "configuration are passed into the callback."
          "\n:param after_apply: Invoked after a pass is applied. "
          "The CompilationUnit and a summary of the pass "
          "configuration are passed into the callback."
          "\n:return: True if pass modified the circuit, else False",
          py::arg("circuit"), py::arg("before_apply"), py::arg("after_apply"))
      .def("__str__", [](const BasePass &) { return "<tket::BasePass>"; })
      .def("__repr__", &BasePass::to_string)
      .def(
          "to_dict", &BasePass::get_config,
          ":return: A JSON serializable dictionary representation of the Pass.")
      .def_static(
          "from_dict", [](const json &j) { return j.get<PassPtr>(); },
          "Construct a new Pass instance from a JSON serializable dictionary "
          "representation.")
      .def(py::pickle(
          [](py::object self) {  // __getstate__
            return py::make_tuple(self.attr("to_dict")());
          },
          [](const py::tuple &t) {  // __setstate__
            const json j = t[0].cast<json>();
            return j.get<PassPtr>();
          }));
  py::class_<SequencePass, std::shared_ptr<SequencePass>, BasePass>(
      m, "SequencePass", "A sequence of compilation passes.")
      .def(
          py::init<const std::vector<PassPtr> &>(),
          "Construct from a list of compilation passes arranged in "
          "order of application.",
          py::arg("pass_list"))
      .def("__str__", [](const BasePass &) { return "<tket::SequencePass>"; })
      .def(
          "get_sequence", &SequencePass::get_sequence,
          ":return: The underlying sequence of passes.");
  py::class_<RepeatPass, std::shared_ptr<RepeatPass>, BasePass>(
      m, "RepeatPass", "Repeat a pass until it has no effect.")
      .def(
          py::init<const PassPtr &>(), "Construct from a compilation pass.",
          py::arg("compilation_pass"))
      .def("__str__", [](const RepeatPass &) { return "<tket::BasePass>"; })
      .def(
          "get_pass", &RepeatPass::get_pass,
          ":return: The underlying compilation pass.");
  py::class_<
      RepeatWithMetricPass, std::shared_ptr<RepeatWithMetricPass>, BasePass>(
      m, "RepeatWithMetricPass",
      "Repeat a compilation pass until the given metric stops "
      "decreasing.")
      .def(
          py::init<const PassPtr &, const Transform::Metric &>(),
          "Construct from a compilation pass and a metric function.",
          py::arg("compilation_pass"), py::arg("metric"))
      .def(
          "__str__",
          [](const BasePass &) { return "<tket::RepeatWithMetricPass>"; })
      .def(
          "get_pass", &RepeatWithMetricPass::get_pass,
          ":return: The underlying compilation pass.")
      .def(
          "get_metric", &RepeatWithMetricPass::get_metric,
          ":return: The underlying metric.");
  py::class_<
      RepeatUntilSatisfiedPass, std::shared_ptr<RepeatUntilSatisfiedPass>,
      BasePass>(
      m, "RepeatUntilSatisfiedPass",
      "Repeat a compilation pass until a predicate on the circuit is "
      "satisfied.")
      .def(
          py::init<const PassPtr &, const PredicatePtr &>(),
          "Construct from a compilation pass and a predicate.",
          py::arg("compilation_pass"), py::arg("predicate"))
      .def(
          py::init<
              const PassPtr &, const std::function<bool(const Circuit &)> &>(),
          "Construct from a compilation pass and a user-defined "
          "function from :py:class:`Circuit` to `bool`.",
          py::arg("compilation_pass"), py::arg("check_function"))
      .def(
          "__str__",
          [](const BasePass &) { return "<tket::RepeatUntilSatisfiedPass>"; })
      .def(
          "get_pass", &RepeatUntilSatisfiedPass::get_pass,
          ":return: The underlying compilation pass.")
      .def(
          "get_predicate", &RepeatUntilSatisfiedPass::get_predicate,
          ":return: The underlying predicate.");

  /* Pass library */
  m.def(
      "KAKDecomposition", &KAKDecomposition,
      "Construct an optimisation pass that performs a Cartan/KAK "
      "Decomposition "
      "for 2 qubit gate sequences, well explained in Robert Tucci "
      "(https://arxiv.org/abs/quant-ph/0507171). "
      "Given a circuit with CXs, SWAPs and any single-qubit gates, "
      "produces a circuit with the same gates. "
      "The optional parameter `cx_fidelity` must be between 0 and 1 and "
      "should be the estimated CX gate fidelity of the noisy quantum "
      "hardware to be compiled for. It is used to trade off decomposition "
      "accuracy for smaller CX gate count. "
      "If `cx_fidelity` < 1, the resulting decomposition will optimise "
      "the expected circuit fidelity. In this case, the output circuit "
      "will not be logically equivalent."
      "This will not preserve CX orientation."
      "\n\n:param cx_fidelity: The estimated gate fidelity"
      "\n:return: a KAK Decomposition pass using the given CX gate "
      "fidelity",
      py::arg("cx_fidelity") = 1.);
  // TODO: Add NormalisedTK2 predicate support.
  m.def(
      "DecomposeTK2",
      [](const py::kwargs &kwargs) {
        return DecomposeTK2(get_fidelities(kwargs));
      },
      "Decompose each TK2 gate into two-qubit gates."
      "\n\nGate fidelities can be passed as keyword arguments to perform "
      "noise-aware decompositions. If the fidelities of several gate types "
      "are provided, the best will be chosen.\n\n"
      "We currently support `CX_fidelity`, `ZZMax_fidelity` and "
      "`ZZPhase_fidelity`. If provided, the `CX` and `ZZMax` fidelities "
      "must be given by a single floating point fidelity. The `ZZPhase` "
      "fidelity is given as a lambda float -> float, mapping a ZZPhase "
      "angle parameter to its fidelity. These parameters will be used "
      "to return the optimal decomposition of each TK2 gate, taking "
      "noise into consideration.\n\n"
      "If no fidelities are provided, the TK2 gates will be decomposed "
      "exactly using CX gates.\n\n"
      "All TK2 gate parameters must be normalised to the Weyl chamber."
      "\n\nIf the TK2 angles are symbolic values, the decomposition will "
      "be exact (i.e. not noise-aware). It is not possible in general "
      "to obtain optimal decompositions for arbitrary symbolic parameters, "
      "so consider substituting for concrete values if possible.");
  m.def(
      "ThreeQubitSquash", &ThreeQubitSquash,
      "Squash three-qubit subcircuits into subcircuits having fewer CX gates, "
      "when possible, and apply Clifford simplification."
      "\n\n:param allow_swaps: whether to allow implicit wire swaps",
      py::arg("allow_swaps") = true);
  m.def(
      "CommuteThroughMultis", &CommuteThroughMultis,
      "Moves single-qubit operations past multi-qubit operations that they "
      "commute with, towards the front of the circuit.");
  m.def(
      "DecomposeArbitrarilyControlledGates",
      &DecomposeArbitrarilyControlledGates,
      "Decomposes CnX and CnRy gates into Ry, CX, H, T and Tdg gates.");
  m.def(
      "DecomposeBoxes", &DecomposeBoxes,
      "Replaces all boxes by their decomposition into circuits.");
  m.def(
      "DecomposeClassicalExp", &DecomposeClassicalExp,
      "Replaces each :py:class:`ClassicalExpBox` by a sequence of "
      "classical gates.");
  m.def(
      "DecomposeMultiQubitsCX", &DecomposeMultiQubitsCX,
      "Converts all multi-qubit gates into CX and single-qubit gates.");
  m.def(
      "GlobalisePhasedX", &GlobalisePhasedX,
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
      "\n\nIf squash=true (default), the `GlobalisePhasedX().apply` method "
      "will always return true. "
      "For squash=false, `apply()` will return true if the circuit was "
      "changed and false otherwise.\n\n"
      "It is not recommended to use this pass with symbolic expressions, as"
      " in certain cases a blow-up in symbolic expression sizes may occur.",
      py::arg("squash") = true);
  m.def(
      "DecomposeSingleQubitsTK1", &DecomposeSingleQubitsTK1,
      "Converts all single-qubit gates into TK1 gates.");
  m.def(
      "PeepholeOptimise2Q", &PeepholeOptimise2Q,
      "Performs peephole optimisation including resynthesis of 2-qubit "
      "gate sequences, and converts to a circuit containing only CX and TK1 "
      "gates.");
  m.def(
      "FullPeepholeOptimise", &FullPeepholeOptimise,
      "Performs peephole optimisation including resynthesis of 2- and 3-qubit "
      "gate sequences, and converts to a circuit containing only CX and TK1 "
      "gates."
      "\n\n:param allow_swaps: whether to allow implicit wire swaps",
      py::arg("allow_swaps") = true);
  m.def("RebaseTket", &RebaseTket, "Converts all gates to CX and TK1.");
  m.def(
      "RemoveRedundancies", &RemoveRedundancies,
      "Removes gate-inverse pairs, merges rotations, removes identity "
      "rotations, and removes redundant gates before measurement. Does not "
      "add any new gate types.\n\n"
      "When merging rotations with the same op group name, the merged"
      "operation keeps the same name.");
  m.def(
      "SynthesiseHQS", &SynthesiseHQS,
      "Optimises and converts a circuit consisting of CX and single-qubit "
      "gates into one containing only ZZMax, PhasedX and Rz.");
  m.def(
      "SynthesiseTK", &SynthesiseTK,
      "Optimises and converts all gates to TK2 and TK1 gates.");
  m.def(
      "SynthesiseTket", &SynthesiseTket,
      "Optimises and converts all gates to CX and TK1 gates.");
  m.def(
      "SynthesiseOQC", &SynthesiseOQC,
      "Optimises and converts all gates to ECR, Rz and SX.");
  m.def(
      "SynthesiseUMD", &SynthesiseUMD,
      "Optimises and converts all gates to XXPhase, PhasedX and Rz.");
  m.def(
      "SquashTK1", &SquashTK1,
      "Squash sequences of single-qubit gates to TK1 gates.");
  m.def(
      "FlattenRegisters", &FlattenRegisters,
      "Merges all quantum and classical registers into their "
      "respective "
      "default registers with contiguous indexing.");
  m.def(
      "SquashCustom", &gen_squash_pass,
      "Squash sequences of single qubit gates from the target gate set "
      "into an optimal form given by `tk1_replacement`."
      "\n\n:param singleqs: The types of single qubit gates in the target "
      "gate set. This pass will only affect sequences of gates that are "
      "already in this set."
      "\n:param tk1_replacement: A function which, given the parameters of "
      "an Rz(a)Rx(b)Rz(c) triple, returns an equivalent circuit in the "
      "desired basis.",
      py::arg("singleqs"), py::arg("tk1_replacement"));
  m.def(
      "DelayMeasures", &DelayMeasures,
      "Commutes Measure operations to the end of the circuit. Throws an "
      "exception when this is not possible because of gates following the "
      "measure which are dependent on either the resulting quantum state "
      "or classical values.");
  m.def(
      "RemoveDiscarded", &RemoveDiscarded,
      "A pass to remove all operations that have no ``OpType.Output`` or "
      "``OpType.ClOutput`` in their causal future (in other words, all "
      "operations whose causal future is discarded).");
  m.def(
      "SimplifyMeasured", &SimplifyMeasured,
      "A pass to replace all 'classical maps' followed by measure "
      "operations whose quantum output is discarded with classical "
      "operations following the measure. (A 'classical map' is a quantum "
      "operation that acts as a permutation of the computational basis "
      "states followed by a diagonal operation.)");
  m.def(
      "RemoveBarriers", &RemoveBarriers,
      "A pass to remove all barrier instructions from the circuit.");

  /* Pass generators */

  m.def(
      "RebaseCustom", &gen_rebase_pass,
      "Construct a custom rebase pass. This pass:\n(1) decomposes "
      "multi-qubit gates not in the set of gate types `gateset` to CX "
      "gates;\n(2) if CX is not in `gateset`, replaces CX gates with "
      "`cx_replacement`;\n(3) converts any single-qubit gates not in the "
      "gate type set to the form "
      ":math:`\\mathrm{Rz}(a)\\mathrm{Rx}(b)\\mathrm{Rz}(c)` (in "
      "matrix-multiplication order, i.e. reverse order in the "
      "circuit);\n(4) applies the `tk1_replacement` function to each of "
      "these triples :math:`(a,b,c)` to generate replacement circuits."
      "\n\n:param gateset: The allowed operations in the rebased circuit."
      "\n:param cx_replacement: The equivalent circuit to replace a CX "
      "gate using two qubit gates from the desired basis (can use any single "
      "qubit OpTypes)."
      "\n:param tk1_replacement: A function which, given the parameters of "
      "an Rz(a)Rx(b)Rz(c) triple, returns an equivalent circuit in the "
      "desired basis."
      "\n:return: a pass that rebases to the given gate set",
      py::arg("gateset"), py::arg("cx_replacement"),
      py::arg("tk1_replacement"));

  m.def(
      "EulerAngleReduction", &gen_euler_pass,
      "Uses Euler angle decompositions to squash all chains of P and Q "
      "rotations, "
      "where P,Q ∈ {Rx,Ry,Rz}. "
      "By default (`strict=False`), this pass will try to decompose the chains "
      "into pairs of -P-Q- or -Q-P- rotations, commuting any third rotation "
      "past multi-qubit gates. "
      "If `strict=True`, all chains will be decomposed to P-Q-P triples "
      "and no further optimisation is performed."
      "\n\n:param q: The type of the Q rotation (Q ∈ {Rx,Ry,Rz})."
      "\n:param p: The type of the P rotation (P ∈ {Rx,Ry,Rz}, P ≠ Q)."
      "\n:param strict: Optionally performs strict P-Q-P Euler decomposition"
      "\n:return: a pass that squashes chains of P and Q rotations",
      py::arg("q"), py::arg("p"), py::arg("strict") = false);

  m.def(
      "CustomRoutingPass", &gen_routing_pass,
      "Construct a pass to route to the connectivity graph of an "
      ":py:class:`Architecture`. Edge direction is ignored."
      "\n:return: a pass that routes to the given device architecture",
      py::arg("arc"), py::arg("config"));

  m.def(
      "RoutingPass", &gen_default_routing_pass,
      "Construct a pass to route to the connectivity graph of an "
      ":py:class:`Architecture`. Edge direction is ignored. "
      "Uses :py:class:`LexiLabellingMethod` and "
      ":py:class:`LexiRouteRoutingMethod`."
      "\n:return: a pass that routes to the given device architecture",
      py::arg("arc"));

  m.def(
      "PlacementPass", &gen_placement_pass,
      ":param placer: The Placement used for relabelling."
      "\n:return: a pass to relabel :py:class:`Circuit` Qubits to "
      ":py:class:`Architecture` Nodes",
      py::arg("placer"));

  m.def(
      "NaivePlacementPass", &gen_naive_placement_pass,
      ":param architecture: The Architecture used for relabelling."
      "\n:return: a pass to relabel :py:class:`Circuit` Qubits to "
      ":py:class:`Architecture` Nodes",
      py::arg("arc"));

  m.def(
      "RenameQubitsPass", &gen_rename_qubits_pass, "Rename some or all qubits.",
      "\n\n:param qubit_map: map from old to new qubit names",
      py::arg("qubit_map"));

  m.def(
      "FullMappingPass", &gen_full_mapping_pass,
      "Construct a pass to relabel :py:class:`Circuit` Qubits to "
      ":py:class:`Architecture` Nodes, and then route to the connectivity "
      "graph "
      "of an :py:class:`Architecture`. Edge direction is ignored."
      "\n\n:param arc: The architecture to use for connectivity information. "
      "\n:param placer: The Placement used for relabelling."
      "\n:param config: Parameters for routing, a list of RoutingMethod, each "
      "method is checked and run if applicable in turn."
      "\n:return: a pass to perform the remapping",
      py::arg("arc"), py::arg("placer"), py::arg("config"));

  m.def(
      "DefaultMappingPass", &gen_default_mapping_pass,
      "Construct a pass to relabel :py:class:`Circuit` Qubits to "
      ":py:class:`Architecture` Nodes, and then route to the connectivity "
      "graph "
      "of the given :py:class: 'Architecture'. Edge direction is ignored. "
      "Placement used "
      "is GraphPlacement."
      "\n\n:param arc: The Architecture used for connectivity information."
      "\n:param delay_measures: Whether to commute measurements to the end "
      "of the circuit, defaulting to true."
      "\n:return: a pass to perform the remapping",
      py::arg("arc"), py::arg("delay_measures") = true);

  m.def(
      "AASRouting", &gen_default_aas_routing_pass,
      "Construct a pass to relabel :py:class:`Circuit` Qubits to "
      ":py:class:`Device` Nodes, and then use architecture-aware "
      "synthesis to route the circuit. In the steps of the pass the circuit "
      "will be converted to CX, Rz, H gateset. The limited connectivity of the "
      ":py:class:`Architecture` is used for the routing. "
      "The direction of the edges is ignored. The placement used "
      "is GraphPlacement. This pass can take a few parameters for the "
      "routing, described below:"
      "\n\n- (unsigned) lookahead=1: parameter for the recursive iteration"
      "\n- (CNotSynthType) cnotsynthtype=CNotSynthType.Rec: CNOT synthesis type"
      "\n\nNB: The circuit needs to have at most as "
      "many qubits as the architecture has nodes. The resulting circuit will "
      "always have the same number of qubits as the architecture has nodes, "
      "even if the input circuit had fewer."
      "\n\n:param arc: target architecture"
      "\n:param \\**kwargs: parameters for routing (described above)"
      "\n:return: a pass to perform the remapping",
      py::arg("arc"));

  m.def(
      "ComposePhasePolyBoxes", &gen_composephasepolybox_pass,
      "Pass to convert a given :py:class:`Circuit` to the CX, Rz, H gateset "
      "and compose "
      "phase polynomial boxes from the groups of the CX+Rz gates."
      "\n\n- (unsigned) min_size=0: minimal number of CX gates in each phase "
      "polynominal box: groups with a smaller number of CX gates are not "
      "affected by this transformation\n"
      "\n:param \\**kwargs: parameters for composition (described above)"
      "\n:return: a pass to perform the composition");

  m.def(
      "CXMappingPass", &gen_cx_mapping_pass_kwargs,
      "Construct a pass to convert all gates to CX, relabel "
      ":py:class:`Circuit` Qubits to :py:class'Architecture' Nodes, route to "
      "the "
      "connectivty graph of a :py:class:`Architecture` and decompose "
      "additional "
      "routing gates (SWAP and BRIDGE) to CX gates."
      "\n\n:param arc: The Architecture used for connectivity information."
      "\n:param placer: The placement used for relabelling."
      "\n:param \\**kwargs: Parameters for routing: "
      "(bool)directed_cx=false, (bool)delay_measures=true"
      "\n:return: a pass to perform the remapping",
      py::arg("arc"), py::arg("placer"));

  m.def(
      "CliffordSimp", &gen_clifford_simp_pass,
      "An optimisation pass that performs a number of rewrite rules for "
      "simplifying Clifford gate sequences, similar to Duncan & Fagan "
      "(https://arxiv.org/abs/1901.10114). "
      "Given a circuit with CXs and any single-qubit gates, produces a "
      "circuit with TK1, CX gates."
      "\n\n:param allow_swaps: dictates whether the rewriting will "
      "disregard CX placement or orientation and introduce wire swaps."
      "\n:return: a pass to perform the rewriting",
      py::arg("allow_swaps") = true);

  m.def(
      "DecomposeSwapsToCXs", &gen_decompose_routing_gates_to_cxs_pass,
      "Construct a pass to decompose SWAP and BRIDGE gates to CX gates, "
      "constraining connectivity to an :py:class:`Architecture`, optionally "
      "taking the directedness of the connectivity graph into account."
      "\n\n:param arc: The architecture to use for connectivity information."
      "\n:param respect_direction: Optionally takes the directedness of "
      "the connectivity graph into account."
      "\n:return: a pass to perform the decomposition",
      py::arg("arc"), py::arg("respect_direction") = false);

  m.def(
      "DecomposeSwapsToCircuit", &gen_user_defined_swap_decomp_pass,
      ":param replacement_circuit: An equivalent circuit to replace a "
      "SWAP gate with in the desired basis."
      "\n:return: a pass to replace all SWAP gates with the given circuit",
      py::arg("replacement_circuit"));
  m.def(
      "OptimisePhaseGadgets", &gen_optimise_phase_gadgets,
      "Construct a pass that synthesises phase gadgets and converts to a "
      "circuit containing only CX and TK1 gates."
      "\n\n:param cx_config: A configuration of CXs to convert phase "
      "gadgets into."
      "\n:return: a pass to perform the synthesis",
      py::arg("cx_config") = CXConfigType::Snake);
  m.def(
      "PauliSimp", &gen_synthesise_pauli_graph,
      "Construct a pass that converts a circuit into a graph of Pauli "
      "gadgets to account for commutation and phase folding, and "
      "resynthesises them as either individual gagdets, pairwise "
      "constructions, or by diagonalising sets of commuting gadgets.\n\n"
      "This pass will not preserve the global phase of the circuit."
      "\n\n:param strat: A synthesis strategy for the Pauli graph."
      "\n:param cx_config: A configuration of CXs to convert Pauli gadgets "
      "into."
      "\n:return: a pass to perform the simplification",
      py::arg("strat") = Transforms::PauliSynthStrat::Sets,
      py::arg("cx_config") = CXConfigType::Snake);
  m.def(
      "GuidedPauliSimp", &gen_special_UCC_synthesis,
      "Applies the ``PauliSimp`` optimisation pass to any region of the "
      "circuit contained within a :py:class:`CircBox`. This can be useful "
      "to focus the synthesis to target specific sets of commuting "
      "operations, rather than the default greedy approach."
      "\n\n:param strat: A synthesis strategy for the Pauli graph."
      "\n:param cx_config: A configuration of CXs to convert Pauli gadgets "
      "into."
      "\n:return: a pass to perform the simplification",
      py::arg("strat") = Transforms::PauliSynthStrat::Sets,
      py::arg("cx_config") = CXConfigType::Snake);
  m.def(
      "PauliSquash", &PauliSquash,
      "Applies :py:meth:`PauliSimp` followed by "
      ":py:meth:`FullPeepholeOptimise`."
      "\n\n:param strat: a synthesis strategy for the Pauli graph"
      "\n:param cx_config: a configuration of CXs to convert Pauli gadgets into"
      "\n:return: a pass to perform the simplification",
      py::arg("strat") = Transforms::PauliSynthStrat::Sets,
      py::arg("cx_config") = CXConfigType::Snake);
  m.def(
      "SimplifyInitial",
      [](bool allow_classical, bool create_all_qubits, bool remove_redundancies,
         std::shared_ptr<const Circuit> xcirc) -> PassPtr {
        PassPtr simpinit = gen_simplify_initial(
            allow_classical ? Transforms::AllowClassical::Yes
                            : Transforms::AllowClassical::No,
            create_all_qubits ? Transforms::CreateAllQubits::Yes
                              : Transforms::CreateAllQubits::No,
            xcirc);
        if (remove_redundancies) {
          std::vector<PassPtr> seq = {simpinit, RemoveRedundancies()};
          return std::make_shared<SequencePass>(seq);
        } else {
          return simpinit;
        }
      },
      "Simplify the circuit using knowledge of qubit state."
      "\n\n:param allow_classical: allow replacement of measurements on "
      "known state with classical set-bit operations"
      "\n:param create_all_qubits: automatically annotate all qubits as "
      "initialized to the zero state"
      "\n:param remove_redundancies: apply a :py:meth:`RemoveRedundancies` "
      "pass after the initial simplification"
      "\n:param xcirc: 1-qubit circuit implementing an X gate in the "
      "transformed circuit (if omitted, an X gate is used)"
      "\n:return: a pass to perform the simplification",
      py::arg("allow_classical") = true, py::arg("create_all_qubits") = false,
      py::arg("remove_redundancies") = true, py::arg("xcirc") = nullptr);
  m.def(
      "ContextSimp",
      [](bool allow_classical, std::shared_ptr<const Circuit> xcirc) {
        return gen_contextual_pass(
            allow_classical ? Transforms::AllowClassical::Yes
                            : Transforms::AllowClassical::No,
            xcirc);
      },
      "Applies simplifications enabled by knowledge of qubit state and "
      "discarded qubits."
      "\n\n:param allow_classical: allow replacement of measurements on "
      "known state with classical set-bit operations"
      "\n:param xcirc: 1-qubit circuit implementing an X gate in the "
      "transformed circuit (if omitted, an X gate is used)"
      "\n:return: a pass to perform the simplification",
      py::arg("allow_classical") = true, py::arg("xcirc") = nullptr);
}

}  // namespace tket
