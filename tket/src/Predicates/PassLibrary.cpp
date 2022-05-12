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

#include "PassLibrary.hpp"

#include <memory>

#include "Circuit/CircPool.hpp"
#include "PassGenerators.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Transformations/BasicOptimisation.hpp"
#include "Transformations/ControlledGates.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/MeasurePass.hpp"
#include "Transformations/OptimisationPass.hpp"
#include "Transformations/Rebase.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/Json.hpp"

namespace tket {

template <typename T>
static PassPtr gate_translation_pass(
    const T &transform, OpTypeSet after_set, bool respect_connectivity,
    const std::string &name) {
  after_set.insert(OpType::Measure);
  after_set.insert(OpType::Collapse);
  after_set.insert(OpType::Reset);
  PredicatePtrMap precons;
  std::type_index ti = typeid(ConnectivityPredicate);
  PredicatePtr out_gateset = std::make_shared<GateSetPredicate>(after_set);
  PredicatePtr max2qb = std::make_shared<MaxTwoQubitGatesPredicate>();
  PredicatePtrMap postcon_spec = {
      CompilationUnit::make_type_pair(out_gateset),
      CompilationUnit::make_type_pair(max2qb)};
  PredicateClassGuarantees g_postcons;
  if (!respect_connectivity)
    g_postcons.insert({ti, Guarantee::Clear});  // synthesis passes do not in
                                                // general preserve connectivity
  PostConditions postcon{postcon_spec, g_postcons, Guarantee::Preserve};
  // record pass config
  nlohmann::json j;
  j["name"] = name;
  PassPtr ptr = std::make_shared<StandardPass>(precons, transform, postcon, j);
  return ptr;
}

const PassPtr &SynthesiseTK() {
  static const PassPtr pp(gate_translation_pass(
      Transforms::synthesise_tk(), {OpType::TK1, OpType::TK2}, true,
      "SynthesiseTK"));
  return pp;
}
const PassPtr &SynthesiseTket() {
  static const PassPtr pp(gate_translation_pass(
      Transforms::synthesise_tket(), {OpType::TK1, OpType::CX}, true,
      "SynthesiseTket"));
  return pp;
}
const PassPtr &SynthesiseHQS() {
  static const PassPtr pp(gate_translation_pass(
      Transforms::synthesise_HQS(),
      {OpType::ZZMax, OpType::PhasedX, OpType::Rz}, false, "SynthesiseHQS"));
  return pp;
}
const PassPtr &SynthesiseOQC() {
  static const PassPtr pp(gate_translation_pass(
      Transforms::synthesise_OQC(), {OpType::Rz, OpType::SX, OpType::ECR}, true,
      "SynthesiseOQC"));
  return pp;
}
const PassPtr &SynthesiseUMD() {
  static const PassPtr pp(gate_translation_pass(
      Transforms::synthesise_UMD(),
      {OpType::XXPhase, OpType::PhasedX, OpType::Rz}, true, "SynthesiseUMD"));
  return pp;
}
const PassPtr &RebaseTket() {
  static const PassPtr pp(gate_translation_pass(
      Transforms::rebase_tket(), {OpType::CX, OpType::TK1}, true,
      "RebaseTket"));
  return pp;
}

const PassPtr &RebaseUFR() {
  static const PassPtr pp(gate_translation_pass(
      Transforms::rebase_UFR(), {OpType::CX, OpType::Rz, OpType::H}, true,
      "RebaseUFR"));
  return pp;
}

const PassPtr &PeepholeOptimise2Q() {
  OpTypeSet after_set = {
      OpType::TK1, OpType::CX, OpType::Measure, OpType::Collapse,
      OpType::Reset};
  PredicatePtrMap precons = {};
  std::type_index ti = typeid(ConnectivityPredicate);
  PredicatePtr out_gateset = std::make_shared<GateSetPredicate>(after_set);
  PredicatePtr max2qb = std::make_shared<MaxTwoQubitGatesPredicate>();
  PredicatePtrMap postcon_spec = {
      CompilationUnit::make_type_pair(out_gateset),
      CompilationUnit::make_type_pair(max2qb)};
  PredicateClassGuarantees g_postcons;
  g_postcons.insert({ti, Guarantee::Clear});
  PostConditions postcon{postcon_spec, g_postcons, Guarantee::Preserve};
  // record pass config
  nlohmann::json j;
  j["name"] = "PeepholeOptimise2Q";
  static const PassPtr pp(std::make_shared<StandardPass>(
      precons, Transforms::peephole_optimise_2q(), postcon, j));
  return pp;
}

const PassPtr &RemoveRedundancies() {
  static const PassPtr pp([]() {
    Transform t = Transforms::remove_redundancies();
    PostConditions postcon = {{}, {}, Guarantee::Preserve};
    PredicatePtrMap precons;
    // record pass config
    nlohmann::json j;
    j["name"] = "RemoveRedundancies";
    return std::make_shared<StandardPass>(precons, t, postcon, j);
  }());
  return pp;
}

const PassPtr &CommuteThroughMultis() {
  static const PassPtr pp([]() {
    Transform t = Transforms::commute_through_multis();
    PostConditions postcon = {{}, {}, Guarantee::Preserve};
    PredicatePtrMap precons;
    // record pass config
    nlohmann::json j;
    j["name"] = "CommuteThroughMultis";
    return std::make_shared<StandardPass>(precons, t, postcon, j);
  }());
  return pp;
}

const PassPtr &DecomposeArbitrarilyControlledGates() {
  static const PassPtr pp([]() {
    Transform t = Transforms::decomp_arbitrary_controlled_gates();
    PredicateClassGuarantees g_postcons = {
        {typeid(GateSetPredicate), Guarantee::Clear}};
    PostConditions postcon = {{}, g_postcons, Guarantee::Preserve};
    PredicatePtrMap precons;
    // record pass config
    nlohmann::json j;
    j["name"] = "DecomposeArbitrarilyControlledGates";
    return std::make_shared<StandardPass>(precons, t, postcon, j);
  }());
  return pp;
}

const PassPtr &DecomposeMultiQubitsCX() {
  static const PassPtr pp([]() {
    Transform t = Transforms::decompose_multi_qubits_CX();
    /* Spits out CX + any single qb gates */
    OpTypeSet ots = {OpType::CX};
    op_signature_t singleq = {EdgeType::Quantum};
    ots.insert(all_projective_types().begin(), all_projective_types().end());
    for (const std::pair<const OpType, OpTypeInfo> &ott : optypeinfo()) {
      if (!ott.second.signature || *ott.second.signature == singleq)
        ots.insert(ott.first);
    }
    PredicatePtrMap precons;
    PredicatePtr outp_gates = std::make_shared<GateSetPredicate>(ots);
    PredicatePtr twoqbpred = std::make_shared<MaxTwoQubitGatesPredicate>();
    PredicatePtrMap spec_postcons = {
        CompilationUnit::make_type_pair(outp_gates),
        CompilationUnit::make_type_pair(twoqbpred)};
    PredicateClassGuarantees g_postcons;
    PostConditions postcon{spec_postcons, g_postcons, Guarantee::Preserve};
    // record pass config
    nlohmann::json j;
    j["name"] = "DecomposeMultiQubitsCX";
    return std::make_shared<StandardPass>(precons, t, postcon, j);
  }());
  return pp;
}

const PassPtr &DecomposeSingleQubitsTK1() {
  static const PassPtr pp([]() {
    Transform t = Transforms::decompose_single_qubits_TK1();
    /* Spits out TK1 + any multi qb gates */
    OpTypeSet ots = {OpType::TK1};
    op_signature_t singleq = {EdgeType::Quantum};
    ots.insert(all_projective_types().begin(), all_projective_types().end());
    for (const std::pair<const OpType, OpTypeInfo> &ott : optypeinfo()) {
      if (!ott.second.signature || *ott.second.signature != singleq)
        ots.insert(ott.first);
    }
    PredicatePtrMap precons;
    PredicatePtr outp_gates = std::make_shared<GateSetPredicate>(ots);
    PredicatePtrMap spec_postcons = {
        CompilationUnit::make_type_pair(outp_gates)};
    PredicateClassGuarantees g_postcons;
    PostConditions postcon{spec_postcons, g_postcons, Guarantee::Preserve};
    // record pass config
    nlohmann::json j;
    j["name"] = "DecomposeSingleQubitsTK1";
    return std::make_shared<StandardPass>(precons, t, postcon, j);
  }());
  return pp;
}

PassPtr ComposePhasePolyBoxes(const unsigned min_size) {
  /**
   * converts a circuit containing all possible gates to a circuit
   * containing only phase poly boxes + H gates (and measure + reset + collapse
   * + barrier)
   * this pass will replace all wire swaps in the given circuit and they will be
   * included in the last or an additional phase poly boxes
   */

  Transform t =
      (Transforms::rebase_UFR() >>
       Transforms::compose_phase_poly_boxes(min_size));

  PredicatePtr noclas = std::make_shared<NoClassicalControlPredicate>();

  PredicatePtrMap precons{CompilationUnit::make_type_pair(noclas)};

  PredicatePtr no_wire_swap = std::make_shared<NoWireSwapsPredicate>();

  PredicatePtrMap s_postcons{
      CompilationUnit::make_type_pair(noclas),
      CompilationUnit::make_type_pair(no_wire_swap)};
  PostConditions postcon{s_postcons, {}, Guarantee::Preserve};

  nlohmann::json j;
  j["name"] = "ComposePhasePolyBoxes";
  j["min_size"] = min_size;

  return std::make_shared<StandardPass>(precons, t, postcon, j);
}

const PassPtr &DecomposeBoxes() {
  static const PassPtr pp([]() {
    Transform t = Transforms::decomp_boxes();
    PredicatePtrMap s_ps;
    /**
     * Preserves Max2QubitGatesPredicate since any box with >2 qubits is
     * already invalid.
     * Preserves ConnectivityPredicate (and DirectednessPredicate) since the
     * verification looks inside CircBoxes and any other boxes with >2
     * qubits are already invalid.
     * Most others are preserved since the predicates look within CircBoxes.
     *
     * Invalidates GateSetPredicate because it doesn't look inside boxes or
     * account for the gate set of their decomposition.
     */
    PredicateClassGuarantees g_postcons = {
        {typeid(GateSetPredicate), Guarantee::Clear},
    };
    PostConditions postcon{s_ps, g_postcons, Guarantee::Preserve};
    // record pass config
    nlohmann::json j;
    j["name"] = "DecomposeBoxes";
    return std::make_shared<StandardPass>(s_ps, t, postcon, j);
  }());
  return pp;
}

const PassPtr &SquashTK1() {
  static const PassPtr pp([]() {
    Transform t = Transforms::squash_1qb_to_tk1();
    PredicatePtrMap s_ps;
    PredicateClassGuarantees g_postcons{
        {typeid(GateSetPredicate), Guarantee::Clear}};
    PostConditions postcon{s_ps, g_postcons, Guarantee::Preserve};
    // record pass config
    nlohmann::json j;
    j["name"] = "SquashTK1";
    return std::make_shared<StandardPass>(s_ps, t, postcon, j);
  }());
  return pp;
}

const PassPtr &DecomposeBridges() {
  static const PassPtr pp([]() {
    Transform t = Transforms::decompose_BRIDGE_to_CX();
    PredicatePtrMap s_ps;
    PredicateClassGuarantees g_postcons{
        {typeid(GateSetPredicate), Guarantee::Clear},
        {typeid(DirectednessPredicate), Guarantee::Clear}};
    PostConditions postcon{s_ps, g_postcons, Guarantee::Preserve};
    nlohmann::json j;
    j["name"] = "DecomposeBridges";
    return std::make_shared<StandardPass>(s_ps, t, postcon, j);
  }());
  return pp;
}

const PassPtr &FlattenRegisters() {
  static const PassPtr pp([]() {
    Transform t =
        Transform([](Circuit &circ, std::shared_ptr<unit_bimaps_t> maps) {
          if (circ.is_simple()) return false;
          unit_map_t qmap = circ.flatten_registers();
          update_maps(maps, qmap, qmap);
          return true;
        });
    PredicatePtrMap s_ps;
    PredicatePtr simple = std::make_shared<DefaultRegisterPredicate>();
    PredicatePtrMap spec_postcons{CompilationUnit::make_type_pair(simple)};
    PredicateClassGuarantees g_postcons{
        {typeid(ConnectivityPredicate), Guarantee::Clear},
        {typeid(DirectednessPredicate), Guarantee::Clear}};
    PostConditions postcon{spec_postcons, g_postcons, Guarantee::Preserve};
    // record pass config
    nlohmann::json j;
    j["name"] = "FlattenRegisters";
    return std::make_shared<StandardPass>(s_ps, t, postcon, j);
  }());
  return pp;
}

const PassPtr &RemoveBarriers() {
  static const PassPtr pp([]() {
    Transform t = Transform([](Circuit &circ) {
      VertexList barriers;
      BGL_FORALL_VERTICES(v, circ.dag, DAG) {
        if (circ.get_OpType_from_Vertex(v) == OpType::Barrier) {
          barriers.push_back(v);
        }
      }
      circ.remove_vertices(
          barriers, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);
      return !barriers.empty();
    });
    const PredicatePtrMap no_precons;
    const PredicatePtr no_barriers = std::make_shared<NoBarriersPredicate>();
    PredicatePtrMap no_barriers_con{
        CompilationUnit::make_type_pair(no_barriers)};
    PredicateClassGuarantees preserve_all;
    PostConditions postcons{no_barriers_con, preserve_all, Guarantee::Preserve};
    nlohmann::json j;
    j["name"] = "RemoveBarriers";
    return std::make_shared<StandardPass>(no_precons, t, postcons, j);
  }());
  return pp;
}

const PassPtr &DelayMeasures() {
  static const PassPtr pp([]() {
    Transform t = Transforms::delay_measures();
    PredicatePtr midmeaspred = std::make_shared<NoMidMeasurePredicate>();
    PredicatePtrMap spec_postcons = {
        CompilationUnit::make_type_pair(midmeaspred)};
    PostConditions postcon = {spec_postcons, {}, Guarantee::Preserve};
    PredicatePtrMap precons;
    // record pass config
    nlohmann::json j;
    j["name"] = "DelayMeasures";
    return std::make_shared<StandardPass>(precons, t, postcon, j);
  }());
  return pp;
}

const PassPtr &RemoveDiscarded() {
  static const PassPtr pp([]() {
    Transform t = Transforms::remove_discarded_ops();
    PredicatePtrMap no_precons;
    PostConditions postcon = {{}, {}, Guarantee::Preserve};
    nlohmann::json j;
    j["name"] = "RemoveDiscarded";
    return std::make_shared<StandardPass>(no_precons, t, postcon, j);
  }());
  return pp;
}

const PassPtr &SimplifyMeasured() {
  static const PassPtr pp([]() {
    Transform t = Transforms::simplify_measured();
    PredicatePtrMap no_precons;

    // GateSetPredicate not preserved because classical gates may be
    // introduced.
    PredicateClassGuarantees g_postcons = {
        {typeid(GateSetPredicate), Guarantee::Clear}};

    PostConditions postcon = {{}, g_postcons, Guarantee::Preserve};
    nlohmann::json j;
    j["name"] = "SimplifyMeasured";
    return std::make_shared<StandardPass>(no_precons, t, postcon, j);
  }());
  return pp;
}

}  // namespace tket
