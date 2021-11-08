// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "PassGenerators.hpp"

#include "ArchAwareSynth/SteinerForest.hpp"
#include "Circuit/CircPool.hpp"
#include "Circuit/Circuit.hpp"
#include "Converters/PhasePoly.hpp"
#include "Mapping/MappingManager.hpp"
#include "Predicates/CompilationUnit.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Predicates/PassLibrary.hpp"
#include "Predicates/Predicates.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/Json.hpp"

namespace tket {

PassPtr gen_rebase_pass(
    const OpTypeSet& multiqs, const Circuit& cx_replacement,
    const OpTypeSet& singleqs,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement) {
  Transform t = Transform::rebase_factory(
      multiqs, cx_replacement, singleqs, tk1_replacement);

  PredicatePtrMap precons;
  OpTypeSet all_types(singleqs);
  all_types.insert(multiqs.begin(), multiqs.end());
  all_types.insert(OpType::Measure);
  all_types.insert(OpType::Collapse);
  all_types.insert(OpType::Reset);
  PredicatePtr postcon1 = std::make_shared<GateSetPredicate>(all_types);
  PredicatePtr postcon2 = std::make_shared<MaxTwoQubitGatesPredicate>();
  std::pair<const std::type_index, PredicatePtr> pair2 =
      CompilationUnit::make_type_pair(postcon1);
  PredicatePtrMap s_postcons{pair2, CompilationUnit::make_type_pair(postcon2)};
  PredicateClassGuarantees g_postcons{{pair2.first, Guarantee::Clear}};
  PostConditions pc = {s_postcons, g_postcons, Guarantee::Preserve};
  // record pass config
  nlohmann::json j;
  j["name"] = "RebaseCustom";
  j["basis_multiqs"] = multiqs;
  j["basis_cx_replacement"] = cx_replacement;
  j["basis_singleqs"] = singleqs;
  j["basis_tk1_replacement"] =
      "SERIALIZATION OF FUNCTIONS IS NOT YET SUPPORTED";
  return std::make_shared<StandardPass>(precons, t, pc, j);
}

PassPtr gen_squash_pass(
    const OpTypeSet& singleqs,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement) {
  Transform t = Transform::squash_factory(singleqs, tk1_replacement);
  PostConditions postcon = {{}, {}, Guarantee::Preserve};
  PredicatePtrMap precons;
  // record pass config
  nlohmann::json j;
  j["name"] = "SquashCustom";
  j["basis_singleqs"] = singleqs;
  j["basis_tk1_replacement"] =
      "SERIALIZATION OF FUNCTIONS IS NOT YET SUPPORTED";
  return std::make_shared<StandardPass>(precons, t, postcon, j);
}

// converting chains of p, q rotations to minimal triplets of p,q-rotations (p,
// q in {Rx,Ry,Rz})
PassPtr gen_euler_pass(const OpType& q, const OpType& p, bool strict) {
  PredicatePtr ccontrol_pred = std::make_shared<NoClassicalControlPredicate>();
  PredicatePtrMap precons{CompilationUnit::make_type_pair(ccontrol_pred)};

  Transform t = Transform::squash_1qb_to_pqp(q, p, strict);
  PostConditions pc{{}, {}, Guarantee::Preserve};
  // record pass config
  nlohmann::json j;
  j["name"] = "EulerAngleReduction";
  j["euler_q"] = q;
  j["euler_p"] = p;
  j["euler_strict"] = strict;
  return std::make_shared<StandardPass>(precons, t, pc, j);
}

PassPtr gen_clifford_simp_pass(bool allow_swaps) {
  // Expects: CX and any single-qubit gates,
  // but does not break if it encounters others
  Transform t = Transform::clifford_simp(allow_swaps);
  PredicatePtr ccontrol_pred = std::make_shared<NoClassicalControlPredicate>();
  PredicatePtrMap precons = {CompilationUnit::make_type_pair(ccontrol_pred)};
  PredicateClassGuarantees g_postcons;
  if (allow_swaps) {
    g_postcons = {
        {typeid(ConnectivityPredicate), Guarantee::Clear},
        {typeid(NoWireSwapsPredicate), Guarantee::Clear},
        {typeid(DirectednessPredicate), Guarantee::Clear}};
  }
  OpTypeSet ots2 = {OpType::CX, OpType::tk1};
  PredicatePtr outp_gates = std::make_shared<GateSetPredicate>(ots2);
  PredicatePtrMap spec_postcons = {CompilationUnit::make_type_pair(outp_gates)};
  PostConditions postcon{spec_postcons, g_postcons, Guarantee::Preserve};

  // record pass config
  nlohmann::json j;
  j["name"] = "CliffordSimp";
  j["allow_swaps"] = allow_swaps;

  return std::make_shared<StandardPass>(precons, t, postcon, j);
}

PassPtr gen_rename_qubits_pass(const std::map<Qubit, Qubit>& qm) {
  Transform t = Transform([=](Circuit& circ) { return circ.rename_units(qm); });
  PredicatePtrMap precons = {};
  PostConditions postcons = {
      {},
      {{typeid(DefaultRegisterPredicate), Guarantee::Clear}},
      Guarantee::Preserve};
  // record pass config
  nlohmann::json j;
  j["name"] = "RenameQubitsPass";
  j["qubit_map"] = qm;

  return std::make_shared<StandardPass>(precons, t, postcons, j);
}

PassPtr gen_placement_pass(const PlacementPtr& placement_ptr) {
  Transform::Transformation trans = [=](Circuit& circ) {
    return placement_ptr->place(circ);
  };
  Transform t = Transform(trans);
  PredicatePtr twoqbpred = std::make_shared<MaxTwoQubitGatesPredicate>();
  PredicatePtr n_qubit_pred = std::make_shared<MaxNQubitsPredicate>(
      placement_ptr->get_architecture_ref().n_uids());

  PredicatePtrMap precons{
      CompilationUnit::make_type_pair(twoqbpred),
      CompilationUnit::make_type_pair(n_qubit_pred)};
  PredicatePtr placement_pred = std::make_shared<PlacementPredicate>(
      placement_ptr->get_architecture_ref());
  PredicatePtrMap s_postcons{CompilationUnit::make_type_pair(placement_pred)};
  PostConditions pc{s_postcons, {}, Guarantee::Preserve};
  // record pass config
  nlohmann::json j;
  j["name"] = "PlacementPass";
  j["placement"] = placement_ptr;
  return std::make_shared<StandardPass>(precons, t, pc, j);
}

PassPtr gen_full_mapping_pass(
    const Architecture& arc, const PlacementPtr& placement_ptr,
    const std::vector<RoutingMethodPtr>& config) {
  return gen_placement_pass(placement_ptr) >> gen_routing_pass(arc, config);
}

PassPtr gen_default_mapping_pass(const Architecture& arc) {
  PlacementPtr pp = std::make_shared<GraphPlacement>(arc);
  RoutingMethodPtr rmw = std::make_shared<LexiRouteRoutingMethod>(100);
  return gen_full_mapping_pass(arc, pp, {rmw});
}

PassPtr gen_cx_mapping_pass(
    const Architecture& arc, const PlacementPtr& placement_ptr,
    const std::vector<RoutingMethodPtr>& config, bool directed_cx,
    bool delay_measures) {
  PassPtr rebase_pass = gen_rebase_pass(
      {OpType::CX}, CircPool::CX(), all_single_qubit_types(),
      Transform::tk1_to_tk1);

  PassPtr return_pass =
      rebase_pass >> gen_full_mapping_pass(arc, placement_ptr, config);
  if (delay_measures) return_pass = return_pass >> DelayMeasures();
  return_pass = return_pass >> rebase_pass >>
                gen_decompose_routing_gates_to_cxs_pass(arc, directed_cx);
  return return_pass;
}

PassPtr gen_routing_pass(
    const Architecture& arc, const std::vector<RoutingMethodPtr>& config) {
  Transform::Transformation trans = [=](Circuit& circ) {
    MappingManager mm(std::make_shared<Architecture>(arc));
    return mm.route_circuit(circ, config);
  };
  Transform t = Transform(trans);

  PredicatePtr twoqbpred = std::make_shared<MaxTwoQubitGatesPredicate>();
  PredicatePtr placedpred = std::make_shared<PlacementPredicate>(arc);
  PredicatePtr n_qubit_pred =
      std::make_shared<MaxNQubitsPredicate>(arc.n_uids());
  PredicatePtrMap precons{
      CompilationUnit::make_type_pair(placedpred),
      CompilationUnit::make_type_pair(twoqbpred),
      CompilationUnit::make_type_pair(n_qubit_pred)};

  PredicatePtr postcon1 = std::make_shared<ConnectivityPredicate>(arc);
  std::pair<const std::type_index, PredicatePtr> pair1 =
      CompilationUnit::make_type_pair(postcon1);
  PredicatePtr postcon2 = std::make_shared<NoWireSwapsPredicate>();
  PredicatePtrMap s_postcons{pair1, CompilationUnit::make_type_pair(postcon2)};
  std::type_index gateset_ti = typeid(GateSetPredicate);
  // clears all Routing predicates apart from the ones which are specified
  // clears all GateSet predicates (inserts SWAPs and BRIDGEs)
  PredicateClassGuarantees g_postcons{
      {pair1.first, Guarantee::Clear},
      {gateset_ti, Guarantee::Clear},
      {typeid(MaxTwoQubitGatesPredicate), Guarantee::Clear}};
  PostConditions pc{s_postcons, g_postcons, Guarantee::Preserve};

  // record pass config
  nlohmann::json j;
  j["name"] = "RoutingPass";
  nlohmann::json config_array;
  for (const auto& con : config) {
    config_array.push_back(*con);
  }
  j["routing_config"] = config_array;
  j["architecture"] = arc;

  return std::make_shared<StandardPass>(precons, t, pc, j);
}

PassPtr gen_placement_pass_phase_poly(const Architecture& arc) {
  Transform::Transformation trans = [=](Circuit& circ) {
    if (arc.n_uids() < circ.n_qubits()) {
      throw CircuitInvalidity(
          "Circuit has more qubits than the architecture has nodes.");
    }

    qubit_vector_t q_vec = circ.all_qubits();
    std::map<Qubit, Node> qubit_to_nodes;
    unsigned counter = 0;

    for (Node x : arc.get_all_uids_set()) {
      if (counter < circ.n_qubits()) {
        qubit_to_nodes.insert({q_vec[counter], x});
        ++counter;
      }
    }

    circ.rename_units(qubit_to_nodes);

    return true;
  };
  Transform t = Transform(trans);

  PredicatePtrMap precons{};
  PredicatePtr placement_pred = std::make_shared<PlacementPredicate>(arc);
  PredicatePtr n_qubit_pred =
      std::make_shared<MaxNQubitsPredicate>(arc.n_uids());
  PredicatePtrMap s_postcons{
      CompilationUnit::make_type_pair(placement_pred),
      CompilationUnit::make_type_pair(n_qubit_pred)};
  PostConditions pc{s_postcons, {}, Guarantee::Preserve};
  // record pass config
  nlohmann::json j;
  j["name"] = "PlacementPass";
  PlacementPtr pp = std::make_shared<GraphPlacement>(arc);
  j["params"]["placement"] = pp;
  return std::make_shared<StandardPass>(precons, t, pc, j);
}

PassPtr aas_routing_pass(
    const Architecture& arc, const unsigned lookahead,
    const aas::CNotSynthType cnotsynthtype) {
  Transform::Transformation trans = [=](Circuit& circ) {
    // check input:
    if (lookahead == 0) {
      throw std::logic_error("lookahead must be > 0");
    }
    if (arc.n_uids() < circ.n_qubits()) {
      throw CircuitInvalidity(
          "Circuit has more qubits than the architecture has nodes.");
    }

    qubit_vector_t all_qu = circ.all_qubits();

    Circuit input_circ = circ;

    VertexList bin;

    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (!circ.detect_boundary_Op(v)) {
        bin.push_back(v);
      }
    }

    circ.remove_vertices(
        bin, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);

    for (const Command& com : input_circ) {
      OpType ot = com.get_op_ptr()->get_type();
      unit_vector_t qbs = com.get_args();
      switch (ot) {
        case OpType::PhasePolyBox: {
          Op_ptr op = com.get_op_ptr();
          const PhasePolyBox& b = static_cast<const PhasePolyBox&>(*op);

          Circuit b_circ(*b.to_circuit());

          PhasePolyBox ppb(b_circ);

          Circuit result =
              aas::phase_poly_synthesis(arc, ppb, lookahead, cnotsynthtype);

          for (const Command& res_com : result) {
            OpType ot = res_com.get_op_ptr()->get_type();
            unit_vector_t res_qbs = res_com.get_args();

            switch (ot) {
              case OpType::CX: {
                circ.add_op<Node>(ot, {Node(res_qbs[0]), Node(res_qbs[1])});
                break;
              }
              case OpType::Rz: {
                auto angle = res_com.get_op_ptr()->get_params().at(0);
                circ.add_op<Node>(ot, angle, {Node(res_qbs[0])});
                break;
              }
              default: {
                TKET_ASSERT(!"Invalid gate type in phase poly box");
              }
            }
          }
          break;
        }
        case OpType::H: {
          circ.add_op<Node>(ot, {Node(qbs[0])});
          break;
        }
        case OpType::Collapse: {
          circ.add_op<Node>(ot, {Node(qbs[0])});
          break;
        }
        case OpType::Reset: {
          circ.add_op<Node>(ot, {Node(qbs[0])});
          break;
        }
        case OpType::Measure: {
          circ.add_measure(Node(qbs[0]), Bit(qbs[1]));
          break;
        }
        default: {
          throw NotImplemented(
              "Invalid gate in input of aas routing: " +
              com.get_op_ptr()->get_name());
        }
      }
    }

    return true;
  };
  Transform t = Transform(trans);

  PredicatePtr placedpred = std::make_shared<PlacementPredicate>(arc);
  PredicatePtr n_qubit_pred =
      std::make_shared<MaxNQubitsPredicate>(arc.n_uids());
  PredicatePtrMap precons{
      CompilationUnit::make_type_pair(placedpred),
      CompilationUnit::make_type_pair(n_qubit_pred)};

  PredicatePtr postcon1 = std::make_shared<ConnectivityPredicate>(arc);
  std::pair<const std::type_index, PredicatePtr> pair1 =
      CompilationUnit::make_type_pair(postcon1);
  PredicatePtr postcon2 = std::make_shared<NoWireSwapsPredicate>();
  PredicatePtrMap s_postcons{pair1, CompilationUnit::make_type_pair(postcon2)};
  std::type_index gateset_ti = typeid(GateSetPredicate);
  PredicateClassGuarantees g_postcons{
      {pair1.first, Guarantee::Clear}, {gateset_ti, Guarantee::Clear}};
  PostConditions pc{s_postcons, g_postcons, Guarantee::Preserve};

  nlohmann::json j;
  j["name"] = "AASRoutingPass";
  j["architecture"] = arc;

  return std::make_shared<StandardPass>(precons, t, pc, j);
}

PassPtr gen_full_mapping_pass_phase_poly(
    const Architecture& arc, const unsigned lookahead,
    const aas::CNotSynthType cnotsynthtype) {
  return RebaseUFR() >> ComposePhasePolyBoxes() >>
         gen_placement_pass_phase_poly(arc) >>
         aas_routing_pass(arc, lookahead, cnotsynthtype);
}

PassPtr gen_directed_cx_routing_pass(
    const Architecture& arc, const std::vector<RoutingMethodPtr>& config) {
  OpTypeSet multis = {OpType::CX, OpType::BRIDGE, OpType::SWAP};
  return gen_routing_pass(arc, config) >>
         gen_rebase_pass(
             multis, CircPool::CX(), all_single_qubit_types(),
             Transform::tk1_to_tk1) >>
         gen_decompose_routing_gates_to_cxs_pass(arc, true);
}

PassPtr gen_decompose_routing_gates_to_cxs_pass(
    const Architecture& arc, bool directed) {
  PredicateClassGuarantees g_postcons{
      {typeid(GateSetPredicate), Guarantee::Clear}};
  PredicatePtrMap precons;
  PredicatePtrMap s_postcons;
  Transform t = Transform::decompose_SWAP_to_CX(arc) >>
                Transform::decompose_BRIDGE_to_CX() >>
                Transform::remove_redundancies();
  if (directed) {
    OpTypeSet out_optypes{all_single_qubit_types()};
    out_optypes.insert(OpType::CX);
    OpTypeSet in_optypes = out_optypes;
    in_optypes.insert(OpType::SWAP);
    in_optypes.insert(OpType::BRIDGE);
    PredicatePtr twoqbpred = std::make_shared<MaxTwoQubitGatesPredicate>();
    PredicatePtr connected = std::make_shared<ConnectivityPredicate>(arc);
    PredicatePtr wireswaps = std::make_shared<NoWireSwapsPredicate>();
    PredicatePtr directedpred = std::make_shared<DirectednessPredicate>(arc);
    PredicatePtr ingates = std::make_shared<GateSetPredicate>(in_optypes);
    PredicatePtr outgates = std::make_shared<GateSetPredicate>(out_optypes);
    precons = {
        CompilationUnit::make_type_pair(connected),
        CompilationUnit::make_type_pair(wireswaps),
        CompilationUnit::make_type_pair(ingates)};
    s_postcons = {
        CompilationUnit::make_type_pair(directedpred),
        CompilationUnit::make_type_pair(outgates),
        CompilationUnit::make_type_pair(twoqbpred)};
    t = t >> Transform::decompose_CX_directed(arc) >>
        Transform::remove_redundancies();
  }
  PostConditions pc{s_postcons, g_postcons, Guarantee::Preserve};
  // record pass config
  nlohmann::json j;
  j["name"] = "DecomposeSwapsToCXs";
  j["directed"] = directed;
  j["architecture"] = arc;
  return std::make_shared<StandardPass>(precons, t, pc, j);
}

PassPtr gen_user_defined_swap_decomp_pass(const Circuit& replacement_circ) {
  Transform t = Transform::decompose_SWAP(replacement_circ);
  PredicateClassGuarantees g_postcons{
      {typeid(GateSetPredicate), Guarantee::Clear}};
  PostConditions pc{{}, g_postcons, Guarantee::Preserve};
  PredicatePtrMap precons;

  // record pass config
  nlohmann::json j;
  j["name"] = "DecomposeSwapsToCircuit";
  j["swap_replacement"] = replacement_circ;
  return std::make_shared<StandardPass>(precons, t, pc, j);
}

PassPtr KAKDecomposition(double cx_fidelity) {
  Transform t = Transform::two_qubit_squash(cx_fidelity);
  PredicatePtr ccontrol_pred = std::make_shared<NoClassicalControlPredicate>();
  OpTypeSet ots{all_single_qubit_types()};
  ots.insert(OpType::SWAP);
  ots.insert(OpType::CX);
  PredicatePtr gate_pred = std::make_shared<GateSetPredicate>(ots);
  PredicatePtrMap precons{
      CompilationUnit::make_type_pair(ccontrol_pred),
      CompilationUnit::make_type_pair(gate_pred)};
  PredicateClassGuarantees g_postcons = {
      {typeid(DirectednessPredicate), Guarantee::Clear},
      {typeid(CliffordCircuitPredicate), Guarantee::Clear}};
  PostConditions postcon = {{}, g_postcons, Guarantee::Preserve};

  // record pass config
  nlohmann::json j;
  j["name"] = "KAKDecomposition";
  j["fidelity"] = cx_fidelity;

  return std::make_shared<StandardPass>(precons, t, postcon, j);
}

PassPtr ThreeQubitSquash(bool allow_swaps) {
  Transform t = Transform::two_qubit_squash() >>
                Transform::three_qubit_squash() >>
                Transform::clifford_simp(allow_swaps);
  OpTypeSet ots{all_single_qubit_types()};
  ots.insert(OpType::CX);
  PredicatePtr gate_pred = std::make_shared<GateSetPredicate>(ots);
  PredicatePtrMap precons{CompilationUnit::make_type_pair(gate_pred)};
  PredicateClassGuarantees g_postcons = {
      {typeid(DirectednessPredicate), Guarantee::Clear},
      {typeid(CliffordCircuitPredicate), Guarantee::Clear}};
  PostConditions postcon = {{}, g_postcons, Guarantee::Preserve};
  nlohmann::json j;
  j["name"] = "ThreeQubitSquash";
  j["allow_swaps"] = allow_swaps;
  return std::make_shared<StandardPass>(precons, t, postcon, j);
}

PassPtr FullPeepholeOptimise(bool allow_swaps) {
  OpTypeSet after_set = {
      OpType::tk1, OpType::CX, OpType::Measure, OpType::Collapse,
      OpType::Reset};
  PredicatePtrMap precons = {};
  PredicatePtr out_gateset = std::make_shared<GateSetPredicate>(after_set);
  PredicatePtr max2qb = std::make_shared<MaxTwoQubitGatesPredicate>();
  PredicatePtrMap postcon_spec = {
      CompilationUnit::make_type_pair(out_gateset),
      CompilationUnit::make_type_pair(max2qb)};
  PredicateClassGuarantees g_postcons = {
      {typeid(ConnectivityPredicate), Guarantee::Clear}};
  PostConditions postcon = {postcon_spec, g_postcons, Guarantee::Preserve};
  nlohmann::json j;
  j["name"] = "FullPeepholeOptimise";
  j["allow_swaps"] = allow_swaps;
  return std::make_shared<StandardPass>(
      precons, Transform::full_peephole_optimise(allow_swaps), postcon, j);
}

PassPtr gen_optimise_phase_gadgets(CXConfigType cx_config) {
  Transform t = Transform::optimise_via_PhaseGadget(cx_config);
  PredicatePtr ccontrol_pred = std::make_shared<NoClassicalControlPredicate>();
  PredicatePtrMap precons{CompilationUnit::make_type_pair(ccontrol_pred)};
  OpTypeSet after_set{
      OpType::Measure, OpType::Collapse, OpType::Reset, OpType::tk1,
      OpType::CX};
  std::type_index ti = typeid(ConnectivityPredicate);
  PredicatePtr out_gateset = std::make_shared<GateSetPredicate>(after_set);
  PredicatePtr max2qb = std::make_shared<MaxTwoQubitGatesPredicate>();
  PredicatePtrMap postcon_spec = {
      CompilationUnit::make_type_pair(out_gateset),
      CompilationUnit::make_type_pair(max2qb)};
  PredicateClassGuarantees g_postcons{{ti, Guarantee::Clear}};
  PostConditions postcon{postcon_spec, g_postcons, Guarantee::Preserve};

  // record pass config
  nlohmann::json j;
  j["name"] = "OptimisePhaseGadgets";
  j["cx_config"] = cx_config;

  return std::make_shared<StandardPass>(precons, t, postcon, j);
}

PassPtr gen_pairwise_pauli_gadgets(CXConfigType cx_config) {
  Transform t = Transform::pairwise_pauli_gadgets(cx_config);
  PredicatePtr ccontrol_pred = std::make_shared<NoClassicalControlPredicate>();
  PredicatePtr simple = std::make_shared<DefaultRegisterPredicate>();
  PredicatePtrMap precons = {
      CompilationUnit::make_type_pair(simple),
      CompilationUnit::make_type_pair(ccontrol_pred)};
  PredicateClassGuarantees g_postcons = {
      {typeid(ConnectivityPredicate), Guarantee::Clear},
      {typeid(NoWireSwapsPredicate), Guarantee::Clear}};
  OpTypeSet ots2 = {OpType::CX, OpType::Z,  OpType::X,  OpType::S,
                    OpType::V,  OpType::U1, OpType::U2, OpType::U3};
  ots2.insert(all_projective_types().begin(), all_projective_types().end());
  PredicatePtr outp_gates = std::make_shared<GateSetPredicate>(ots2);
  PredicatePtrMap spec_postcons = {CompilationUnit::make_type_pair(outp_gates)};
  PostConditions postcon{spec_postcons, g_postcons, Guarantee::Preserve};
  nlohmann::json j;
  j["name"] = "OptimisePairwiseGadgets";
  j["cx_config"] = cx_config;
  return std::make_shared<StandardPass>(precons, t, postcon, j);
}

PassPtr gen_synthesise_pauli_graph(
    PauliSynthStrat strat, CXConfigType cx_config) {
  Transform t = Transform::synthesise_pauli_graph(strat, cx_config);
  PredicatePtr ccontrol_pred = std::make_shared<NoClassicalControlPredicate>();
  PredicatePtr mid_pred = std::make_shared<NoMidMeasurePredicate>();
  PredicatePtr wire_pred = std::make_shared<NoWireSwapsPredicate>();
  OpTypeSet ins = {OpType::Z,       OpType::X,           OpType::Y,
                   OpType::S,       OpType::Sdg,         OpType::V,
                   OpType::Vdg,     OpType::H,           OpType::CX,
                   OpType::CY,      OpType::CZ,          OpType::SWAP,
                   OpType::Rz,      OpType::Rx,          OpType::Ry,
                   OpType::T,       OpType::Tdg,         OpType::ZZMax,
                   OpType::ZZPhase, OpType::PhaseGadget, OpType::XXPhase,
                   OpType::YYPhase, OpType::PauliExpBox, OpType::Measure};
  PredicatePtr in_gates = std::make_shared<GateSetPredicate>(ins);
  PredicatePtrMap precons{
      CompilationUnit::make_type_pair(ccontrol_pred),
      CompilationUnit::make_type_pair(mid_pred),
      CompilationUnit::make_type_pair(wire_pred),
      CompilationUnit::make_type_pair(in_gates)};
  PredicateClassGuarantees g_postcons = {
      {typeid(ConnectivityPredicate), Guarantee::Clear},
      {typeid(NoWireSwapsPredicate), Guarantee::Clear}};
  PostConditions postcon{{}, g_postcons, Guarantee::Preserve};

  // record pass config
  nlohmann::json j;
  j["name"] = "PauliSimp";
  j["cx_config"] = cx_config;
  j["pauli_synth_strat"] = strat;

  return std::make_shared<StandardPass>(precons, t, postcon, j);
}

PassPtr gen_special_UCC_synthesis(
    PauliSynthStrat strat, CXConfigType cx_config) {
  Transform t = Transform::special_UCC_synthesis(strat, cx_config);
  PredicatePtr ccontrol_pred = std::make_shared<NoClassicalControlPredicate>();
  PredicatePtrMap precons{CompilationUnit::make_type_pair(ccontrol_pred)};
  PredicateClassGuarantees g_postcons = {
      {typeid(ConnectivityPredicate), Guarantee::Clear},
      {typeid(NoWireSwapsPredicate), Guarantee::Clear}};
  PostConditions postcon{{}, g_postcons, Guarantee::Preserve};

  // record pass config
  nlohmann::json j;
  j["name"] = "GuidedPauliSimp";
  j["cx_config"] = cx_config;
  j["pauli_synth_strat"] = strat;

  return std::make_shared<StandardPass>(precons, t, postcon, j);
}

PassPtr gen_simplify_initial(
    Transform::AllowClassical allow_classical,
    Transform::CreateAllQubits create_all_qubits,
    std::shared_ptr<const Circuit> xcirc) {
  Transform t =
      Transform::simplify_initial(allow_classical, create_all_qubits, xcirc);
  PredicatePtrMap no_precons;
  PredicateClassGuarantees g_postcons;

  // GateSetPredicate not preserved because X gates (or their specified
  // equivalents) may be introduced, and if allow_classical then classical
  // gates may also be introduced.
  g_postcons[typeid(GateSetPredicate)] = Guarantee::Clear;

  PostConditions postcon = {{}, g_postcons, Guarantee::Preserve};
  nlohmann::json j;
  j["name"] = "SimplifyInitial";
  j["allow_classical"] = (allow_classical == Transform::AllowClassical::Yes);
  j["create_all_qubits"] =
      (create_all_qubits == Transform::CreateAllQubits::Yes);
  if (xcirc) j["x_circuit"] = *xcirc;
  return std::make_shared<StandardPass>(no_precons, t, postcon, j);
}

PassPtr gen_contextual_pass(
    Transform::AllowClassical allow_classical,
    std::shared_ptr<const Circuit> xcirc) {
  std::vector<PassPtr> seq = {
      RemoveDiscarded(), SimplifyMeasured(),
      gen_simplify_initial(
          allow_classical, Transform::CreateAllQubits::No, xcirc),
      RemoveRedundancies()};
  return std::make_shared<SequencePass>(seq);
}

PassPtr PauliSquash(PauliSynthStrat strat, CXConfigType cx_config) {
  std::vector<PassPtr> seq = {
      gen_synthesise_pauli_graph(strat, cx_config), FullPeepholeOptimise()};
  return std::make_shared<SequencePass>(seq);
}

}  // namespace tket
