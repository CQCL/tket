// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tket/Transformations/GreedyPauliOptimisation.hpp"

#include <algorithm>

#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/Converters/Converters.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/PauliGraph/PauliGraph.hpp"
#include "tket/Transformations/CliffordOptimisation.hpp"
#include "tket/Transformations/GreedyPauliOptimisationLookupTables.hpp"
#include "tket/Transformations/Transform.hpp"

namespace tket {

namespace Transforms {

namespace GreedyPauliSimp {

static void apply_tqe_to_circ(const TQE& tqe, Circuit& circ) {
  auto [gate_type, a, b] = tqe;
  switch (gate_type) {
    case TQEType::XX:
      circ.add_op<unsigned>(OpType::H, {a});
      circ.add_op<unsigned>(OpType::CX, {a, b});
      circ.add_op<unsigned>(OpType::H, {a});
      break;
    case TQEType::XY:
      circ.add_op<unsigned>(OpType::H, {a});
      circ.add_op<unsigned>(OpType::CY, {a, b});
      circ.add_op<unsigned>(OpType::H, {a});
      break;
    case TQEType::XZ:
      circ.add_op<unsigned>(OpType::CX, {b, a});
      break;
    case TQEType::YX:
      circ.add_op<unsigned>(OpType::H, {b});
      circ.add_op<unsigned>(OpType::CY, {b, a});
      circ.add_op<unsigned>(OpType::H, {b});
      break;
    case TQEType::YY:
      circ.add_op<unsigned>(OpType::V, {a});
      circ.add_op<unsigned>(OpType::CY, {a, b});
      circ.add_op<unsigned>(OpType::Vdg, {a});
      break;
    case TQEType::YZ:
      circ.add_op<unsigned>(OpType::CY, {b, a});
      break;
    case TQEType::ZX:
      circ.add_op<unsigned>(OpType::CX, {a, b});
      break;
    case TQEType::ZY:
      circ.add_op<unsigned>(OpType::CY, {a, b});
      break;
    case TQEType::ZZ:
      circ.add_op<unsigned>(OpType::CZ, {a, b});
      break;
  }
}

// return the sum of the cost increases on remaining tableau nodes
static double default_tableau_tqe_cost(
    const std::vector<PauliNode_ptr>& rows,
    const std::vector<unsigned>& remaining_indices, const TQE& tqe) {
  double cost = 0;
  for (const unsigned& index : remaining_indices) {
    cost += rows[index]->tqe_cost_increase(tqe);
  }
  return cost;
}

// return the weighted sum of the cost increases on remaining nodes
// we discount the weight after each set
static double default_pauliexp_tqe_cost(
    const double discount_rate,
    const std::vector<std::vector<PauliNode_ptr>>& rotation_sets,
    const std::vector<PauliNode_ptr>& rows, const TQE& tqe) {
  double discount = 1 / (1 + discount_rate);
  double weight = 1;
  double exp_cost = 0;
  double tab_cost = 0;
  for (const std::vector<PauliNode_ptr>& rotation_set : rotation_sets) {
    for (const PauliNode_ptr& node : rotation_set) {
      exp_cost += weight * node->tqe_cost_increase(tqe);
    }
    weight *= discount;
  }
  for (const PauliNode_ptr& node : rows) {
    tab_cost += weight * node->tqe_cost_increase(tqe);
  }
  return exp_cost + tab_cost;
}

// given a map from TQE to a vector of costs, select the one with the minimum
// weighted sum of minmax-normalised costs
static TQE minmax_selection(
    const std::map<TQE, std::vector<double>>& tqe_candidates_cost,
    const std::vector<double>& weights) {
  TKET_ASSERT(tqe_candidates_cost.size() > 0);
  size_t n_costs = tqe_candidates_cost.begin()->second.size();
  TKET_ASSERT(n_costs == weights.size());
  // for each cost type, store its min and max
  std::vector<double> mins = tqe_candidates_cost.begin()->second;
  std::vector<double> maxs = tqe_candidates_cost.begin()->second;
  for (const auto& pair : tqe_candidates_cost) {
    TKET_ASSERT(pair.second.size() == n_costs);
    for (unsigned cost_index = 0; cost_index < n_costs; cost_index++) {
      if (pair.second[cost_index] < mins[cost_index]) {
        mins[cost_index] = pair.second[cost_index];
      }
      if (pair.second[cost_index] > maxs[cost_index]) {
        maxs[cost_index] = pair.second[cost_index];
      }
    }
  }
  // valid_indices stores the indices of the costs where min!=max
  std::vector<unsigned> valid_indices;
  for (unsigned cost_index = 0; cost_index < n_costs; cost_index++) {
    if (mins[cost_index] != maxs[cost_index]) {
      valid_indices.push_back(cost_index);
    }
  }
  // if all have the same cost, return the first one
  if (valid_indices.size() == 0) {
    TQE min_tqe = tqe_candidates_cost.begin()->first;
    return min_tqe;
  }
  // if only one cost variable, no need to normalise
  if (valid_indices.size() == 1) {
    auto it = tqe_candidates_cost.begin();
    double min_cost = it->second[valid_indices[0]];
    TQE min_tqe = it->first;
    for (; it != tqe_candidates_cost.end(); it++) {
      if (it->second[valid_indices[0]] < min_cost) {
        min_tqe = it->first;
        min_cost = it->second[valid_indices[0]];
      }
    }
    return min_tqe;
  }
  // find the tqe with the minimum normalised cost
  auto it = tqe_candidates_cost.begin();
  double min_cost = 0;
  TQE min_tqe = it->first;
  // initialise min_cost
  for (const auto& cost_index : valid_indices) {
    min_cost += weights[cost_index] *
                (it->second[cost_index] - mins[cost_index]) /
                (maxs[cost_index] - mins[cost_index]);
  }
  it++;
  // iterate all tqes
  for (; it != tqe_candidates_cost.end(); it++) {
    double cost = 0;
    for (const auto& cost_index : valid_indices) {
      cost += weights[cost_index] *
              (it->second[cost_index] - mins[cost_index]) /
              (maxs[cost_index] - mins[cost_index]);
    }
    if (cost < min_cost) {
      min_cost = cost;
      min_tqe = it->first;
    }
  }
  return min_tqe;
}

static TQE select_pauliexp_tqe(
    const std::map<TQE, std::vector<double>>& tqe_candidates_cost,
    double depth_weight) {
  return minmax_selection(tqe_candidates_cost, {1, depth_weight});
}

static TQE select_tableau_tqe(
    const std::map<TQE, std::vector<double>>& tqe_candidates_cost,
    double depth_weight) {
  return minmax_selection(tqe_candidates_cost, {1, depth_weight});
}

// simple struct that tracks the depth on each qubit
struct DepthTracker {
  std::vector<unsigned> qubit_depth;
  unsigned max_depth;
  DepthTracker(unsigned n) : qubit_depth(n, 0), max_depth(0) {};

  unsigned gate_depth(unsigned a, unsigned b) const {
    return std::max(qubit_depth[a], qubit_depth[b]) + 1;
  };
  void add_1q_gate(unsigned a) {
    qubit_depth[a]++;
    if (qubit_depth[a] > max_depth) {
      max_depth = qubit_depth[a];
    }
  };
  void add_2q_gate(unsigned a, unsigned b) {
    unsigned new_gate_depth = gate_depth(a, b);
    qubit_depth[a] = new_gate_depth;
    qubit_depth[b] = new_gate_depth;
    if (new_gate_depth > max_depth) {
      max_depth = new_gate_depth;
    }
  };
};

/**
 * @brief Synthesise a vector of PauliPropagation
 */
static void tableau_row_nodes_synthesis(
    std::vector<PauliNode_ptr>& rows, Circuit& circ, double depth_weight,
    DepthTracker& depth_tracker) {
  // only consider nodes with a non-zero cost
  std::vector<unsigned> remaining_indices;
  for (unsigned i = 0; i < rows.size(); i++) {
    if (rows[i]->tqe_cost() > 0) {
      remaining_indices.push_back(i);
    }
  }
  while (remaining_indices.size() != 0) {
    // get nodes with min cost
    std::vector<unsigned> min_nodes_indices = {remaining_indices[0]};
    unsigned min_cost = rows[remaining_indices[0]]->tqe_cost();
    for (unsigned i = 1; i < remaining_indices.size(); i++) {
      unsigned node_cost = rows[remaining_indices[i]]->tqe_cost();
      if (node_cost == min_cost) {
        min_nodes_indices.push_back(remaining_indices[i]);
      } else if (node_cost < min_cost) {
        min_nodes_indices = {remaining_indices[i]};
        min_cost = node_cost;
      }
    }
    // for each node with min cost, find the list of tqe gates that can reduce
    // its cost
    std::set<TQE> tqe_candidates;
    TKET_ASSERT(min_nodes_indices.size() > 0);
    for (const unsigned& index : min_nodes_indices) {
      std::vector<TQE> node_reducing_tqes = rows[index]->reduction_tqes();
      tqe_candidates.insert(
          node_reducing_tqes.begin(), node_reducing_tqes.end());
    }
    // for each tqe we compute a vector of cost factors which will
    // be combined to make the final decision.
    // we currently only consider tqe_cost and gate_depth.
    std::map<TQE, std::vector<double>> tqe_candidates_cost;
    for (const TQE& tqe : tqe_candidates) {
      tqe_candidates_cost.insert(
          {tqe,
           {default_tableau_tqe_cost(rows, remaining_indices, tqe),
            static_cast<double>(depth_tracker.gate_depth(
                std::get<1>(tqe), std::get<2>(tqe)))}});
    }
    TKET_ASSERT(tqe_candidates_cost.size() > 0);
    // select the best one
    TQE selected_tqe = select_tableau_tqe(tqe_candidates_cost, depth_weight);
    // apply TQE
    apply_tqe_to_circ(selected_tqe, circ);
    // update depth tracker
    depth_tracker.add_2q_gate(
        std::get<1>(selected_tqe), std::get<2>(selected_tqe));
    // remove finished nodes
    for (unsigned i = remaining_indices.size(); i-- > 0;) {
      unsigned node_index = remaining_indices[i];
      rows[node_index]->update(selected_tqe);
      if (rows[node_index]->tqe_cost() == 0) {
        remaining_indices.erase(remaining_indices.begin() + i);
      }
    }
  }
  // apply local Cliffords
  for (PauliNode_ptr& node_ptr : rows) {
    PauliPropagation& node = dynamic_cast<PauliPropagation&>(*node_ptr);
    auto [q_index, supp_z, supp_x] = node.first_support();
    // transform supp_z,supp_x to Z,X
    std::vector<OpType> optype_list = AA_TO_ZX.at({supp_z, supp_x});
    for (auto it = optype_list.rbegin(); it != optype_list.rend(); ++it) {
      circ.add_op<unsigned>(SQ_CLIFF_DAGGER.at(*it), {q_index});
      node.update(*it, q_index);
    }
    // remove signs
    if (!node.z_sign()) {
      circ.add_op<unsigned>(OpType::X, {q_index});
      node.update(OpType::X, q_index);
    }
    if (!node.x_sign()) {
      circ.add_op<unsigned>(OpType::Z, {q_index});
      node.update(OpType::Z, q_index);
    }
    if (q_index != node.qubit_index()) {
      circ.add_op<unsigned>(OpType::SWAP, {q_index, node.qubit_index()});
      for (PauliNode_ptr& node_ptr2 : rows) {
        node_ptr2->swap(q_index, node.qubit_index());
      }
    }
  }
}

/**
 * @brief Given a vector of sets of PauliRotation, implement any node in the
 * first set where the tqe_cost is zero. Remove implemented nodes and the first
 * set if empty.
 *
 * @param rotation_sets
 * @param circ
 * @return true if the first set is now empty and removed
 * @return false
 */
static void consume_available_rotations(
    std::vector<std::vector<PauliNode_ptr>>& rotation_sets, Circuit& circ,
    DepthTracker& depth_tracker) {
  if (rotation_sets.empty()) {
    return;
  }
  while (true) {
    std::vector<PauliNode_ptr>& first_set = rotation_sets[0];
    for (unsigned i = first_set.size(); i-- > 0;) {
      TKET_ASSERT(first_set[i]->get_type() == PauliNodeType::Rotation);
      PauliRotation& node = dynamic_cast<PauliRotation&>(*first_set[i]);
      if (node.tqe_cost() > 0) continue;
      auto [q_index, supp] = node.first_support();
      depth_tracker.add_1q_gate(q_index);
      OpType rot_type;
      switch (supp) {
        case Pauli::Y: {
          rot_type = OpType::Ry;
          break;
        }
        case Pauli::Z: {
          rot_type = OpType::Rz;
          break;
        }
        case Pauli::X: {
          rot_type = OpType::Rx;
          break;
        }
        default:
          // support can't be Pauli::I
          TKET_ASSERT(false);
      }
      if (node.sign()) {
        circ.add_op<unsigned>(rot_type, node.theta(), {q_index});
      } else {
        circ.add_op<unsigned>(rot_type, -node.theta(), {q_index});
      }
      first_set.erase(first_set.begin() + i);
    }
    if (first_set.empty()) {
      rotation_sets.erase(rotation_sets.begin());
      if (rotation_sets.empty()) {
        return;
      }
    } else {
      return;
    }
  }
}

/**
 * @brief Synthesise a vector of unordered rotation sets
 */
static void pauli_exps_synthesis(
    std::vector<std::vector<PauliNode_ptr>>& rotation_sets,
    std::vector<PauliNode_ptr>& rows, Circuit& circ, double discount_rate,
    double depth_weight, DepthTracker& depth_tracker) {
  while (true) {
    consume_available_rotations(rotation_sets, circ, depth_tracker);
    if (rotation_sets.empty()) break;
    std::vector<PauliNode_ptr>& first_set = rotation_sets[0];
    // get nodes with min cost
    std::vector<unsigned> min_nodes_indices = {0};
    unsigned min_cost = first_set[0]->tqe_cost();
    for (unsigned i = 1; i < first_set.size(); i++) {
      unsigned node_cost = first_set[i]->tqe_cost();
      if (node_cost == min_cost) {
        min_nodes_indices.push_back(i);
      } else if (node_cost < min_cost) {
        min_nodes_indices = {i};
        min_cost = node_cost;
      }
    }
    std::set<TQE> tqe_candidates;
    for (const unsigned& index : min_nodes_indices) {
      std::vector<TQE> node_reducing_tqes = first_set[index]->reduction_tqes();
      tqe_candidates.insert(
          node_reducing_tqes.begin(), node_reducing_tqes.end());
    }
    // for each tqe we compute costs which might subject to normalisation
    std::map<TQE, std::vector<double>> tqe_candidates_cost;
    for (const TQE& tqe : tqe_candidates) {
      tqe_candidates_cost.insert(
          {tqe,
           {default_pauliexp_tqe_cost(discount_rate, rotation_sets, rows, tqe),
            static_cast<double>(depth_tracker.gate_depth(
                std::get<1>(tqe), std::get<2>(tqe)))}});
    }
    // select the best one
    TQE selected_tqe = select_pauliexp_tqe(tqe_candidates_cost, depth_weight);
    // apply TQE
    apply_tqe_to_circ(selected_tqe, circ);
    depth_tracker.add_2q_gate(
        std::get<1>(selected_tqe), std::get<2>(selected_tqe));
    for (std::vector<PauliNode_ptr>& rotation_set : rotation_sets) {
      for (PauliNode_ptr& node : rotation_set) {
        node->update(selected_tqe);
      }
    }
    for (PauliNode_ptr& row : rows) {
      row->update(selected_tqe);
    }
  }
}

// convert a Pauli exponential to a PauliNode_ptr
static PauliNode_ptr get_node_from_exp(
    const std::vector<Pauli>& paulis, const Expr& theta,
    const qubit_vector_t& args, unsigned n) {
  // pad the Paulis
  std::vector<Pauli> string(n, Pauli::I);
  for (unsigned i = 0; i < args.size(); i++) {
    string[args[i].index()[0]] = paulis[i];
  }
  return std::make_shared<PauliRotation>(string, theta);
}

// convert a Clifford tableau to a vector of PauliNode_ptr
static std::vector<PauliNode_ptr> get_nodes_from_tableau(
    const UnitaryRevTableau& tab, unsigned n_qubits) {
  std::vector<PauliNode_ptr> rows;
  for (unsigned i = 0; i < n_qubits; i++) {
    Qubit q(i);
    SpPauliStabiliser z_stab = tab.get_zrow(q);
    SpPauliStabiliser x_stab = tab.get_xrow(q);
    bool z_sign = cast_coeff<quarter_turns_t, Complex>(z_stab.coeff) == 1.;
    bool x_sign = cast_coeff<quarter_turns_t, Complex>(x_stab.coeff) == 1.;
    TKET_ASSERT(z_stab.string.size() == n_qubits);
    std::vector<Pauli> z_string;
    std::vector<Pauli> x_string;
    for (unsigned j = 0; j < n_qubits; j++) {
      z_string.push_back(z_stab.string.at(Qubit(j)));
      x_string.push_back(x_stab.string.at(Qubit(j)));
    }
    rows.push_back(std::make_shared<PauliPropagation>(
        z_string, x_string, z_sign, x_sign, i));
  }
  return rows;
}

// detect trivial pauli exps, if true then return the global phase
static std::pair<bool, Expr> is_trivial_pauliexp(
    const std::vector<Pauli>& paulis, const Expr& theta) {
  if (static_cast<std::size_t>(std::count(
          paulis.begin(), paulis.end(), Pauli::I)) == paulis.size()) {
    // If all identity term
    return {true, -theta / 2};
  }
  if (equiv_0(theta, 2)) {
    if (equiv_0(theta, 4)) {
      return {true, 0};
    } else {
      return {true, -1};
    }
  }
  return {false, 0};
}

Circuit greedy_pauli_set_synthesis(
    const std::vector<SymPauliTensor>& unordered_set, double depth_weight) {
  if (unordered_set.size() == 0) {
    return Circuit();
  }
  unsigned n_qubits = unordered_set[0].string.size();

  Circuit c(n_qubits);
  std::vector<std::vector<PauliNode_ptr>> rotation_sets{{}};

  for (auto& pauli : unordered_set) {
    TKET_ASSERT(pauli.string.size() == n_qubits);
    rotation_sets[0].push_back(
        std::make_shared<PauliRotation>(pauli.string, pauli.coeff));
  }
  UnitaryRevTableau tab(n_qubits);
  std::vector<PauliNode_ptr> rows = get_nodes_from_tableau(tab, n_qubits);
  DepthTracker depth_tracker(n_qubits);
  // synthesise Pauli exps
  pauli_exps_synthesis(rotation_sets, rows, c, 0, depth_weight, depth_tracker);
  // synthesise the tableau
  tableau_row_nodes_synthesis(rows, c, depth_weight, depth_tracker);
  c.replace_SWAPs();
  return c;
}

Circuit greedy_pauli_graph_synthesis(
    const Circuit& circ, double discount_rate, double depth_weight) {
  // c is the circuit we are trying to build
  Circuit c(circ.all_qubits(), circ.all_bits());
  std::optional<std::string> name = circ.get_name();
  if (name != std::nullopt) {
    c.set_name(name.value());
  }
  c.add_phase(circ.get_phase());
  unit_map_t unit_map = c.flatten_registers();
  Circuit measure_circ(c.n_qubits(), c.n_bits());
  Circuit cliff(c.n_qubits());
  // circuit used to iterate the original commands with flattened registers
  Circuit circ_flat(circ);
  circ_flat.flatten_registers();
  std::vector<Command> commands = circ_flat.get_commands();
  // extract the final clifford and the measurement circuits
  for (const Command& cmd : commands) {
    OpType optype = cmd.get_op_ptr()->get_type();
    switch (optype) {
      case OpType::Measure: {
        measure_circ.add_op<UnitID>(OpType::Measure, cmd.get_args());
        break;
      }
      default: {
        if (optype == OpType::PauliExpBox ||
            optype == OpType::PauliExpPairBox ||
            optype == OpType::PauliExpCommutingSetBox)
          break;
        TKET_ASSERT(is_clifford_type(optype) && is_gate_type(optype));
        cliff.add_op<UnitID>(optype, cmd.get_args());
      }
    }
  }
  std::vector<std::vector<PauliNode_ptr>> rotation_sets;
  unsigned n_qubits = c.n_qubits();
  UnitaryRevTableau tab = circuit_to_unitary_rev_tableau(cliff);
  // convert the tableau into a set of nodes
  std::vector<PauliNode_ptr> rows = get_nodes_from_tableau(tab, n_qubits);
  // extract the Pauli exps
  for (const Command& cmd : commands) {
    OpType optype = cmd.get_op_ptr()->get_type();
    switch (optype) {
      case OpType::PauliExpBox: {
        const PauliExpBox& pbox =
            static_cast<const PauliExpBox&>(*cmd.get_op_ptr());
        const Expr phase = pbox.get_phase();
        const std::vector<Pauli> paulis = pbox.get_paulis();
        auto [trivial, global_phase] = is_trivial_pauliexp(paulis, phase);
        if (trivial) {
          c.add_phase(global_phase);
        } else {
          rotation_sets.push_back(
              {get_node_from_exp(paulis, phase, cmd.get_qubits(), n_qubits)});
        }
        break;
      }
      case OpType::PauliExpPairBox: {
        const PauliExpPairBox& pbox =
            static_cast<const PauliExpPairBox&>(*cmd.get_op_ptr());
        const auto [paulis1, paulis2] = pbox.get_paulis_pair();
        const auto [phase1, phase2] = pbox.get_phase_pair();
        auto [trivial1, global_phase1] = is_trivial_pauliexp(paulis1, phase1);
        auto [trivial2, global_phase2] = is_trivial_pauliexp(paulis2, phase2);
        std::vector<PauliNode_ptr> rotation_set;
        if (trivial1) {
          c.add_phase(global_phase1);
        } else {
          rotation_set.push_back(
              get_node_from_exp(paulis1, phase1, cmd.get_qubits(), n_qubits));
        }
        if (trivial2) {
          c.add_phase(global_phase2);
        } else {
          rotation_set.push_back(
              get_node_from_exp(paulis2, phase2, cmd.get_qubits(), n_qubits));
        }
        if (!rotation_set.empty()) {
          rotation_sets.push_back(rotation_set);
        }
        break;
      }
      case OpType::PauliExpCommutingSetBox: {
        const PauliExpCommutingSetBox& pbox =
            static_cast<const PauliExpCommutingSetBox&>(*cmd.get_op_ptr());
        const std::vector<SymPauliTensor> gadgets = pbox.get_pauli_gadgets();
        std::vector<PauliNode_ptr> rotation_set;
        for (const SymPauliTensor& pt : gadgets) {
          const std::vector<Pauli> paulis = pt.string;
          const Expr phase = pt.coeff;
          auto [trivial, global_phase] = is_trivial_pauliexp(paulis, phase);
          if (trivial) {
            c.add_phase(global_phase);
          } else {
            rotation_set.push_back(
                get_node_from_exp(paulis, phase, cmd.get_qubits(), n_qubits));
          }
        }
        if (rotation_set.size() > 0) {
          rotation_sets.push_back(rotation_set);
        }
        break;
      }
      default:
        break;
    }
  }

  DepthTracker depth_tracker(n_qubits);
  // synthesise Pauli exps
  pauli_exps_synthesis(
      rotation_sets, rows, c, discount_rate, depth_weight, depth_tracker);
  // synthesise the tableau
  tableau_row_nodes_synthesis(rows, c, depth_weight, depth_tracker);
  unit_map_t rev_unit_map;
  for (const auto& pair : unit_map) {
    rev_unit_map.insert({pair.second, pair.first});
  }
  c.append(measure_circ);
  c.rename_units(rev_unit_map);
  c.replace_SWAPs();
  return c;
}

}  // namespace GreedyPauliSimp

Transform greedy_pauli_optimisation(double discount_rate, double depth_weight) {
  return Transform([discount_rate, depth_weight](Circuit& circ) {
    synthesise_pauli_graph(PauliSynthStrat::Sets, CXConfigType::Snake)
        .apply(circ);
    circ = GreedyPauliSimp::greedy_pauli_graph_synthesis(
        circ, discount_rate, depth_weight);
    singleq_clifford_sweep().apply(circ);
    return true;
  });
}

}  // namespace Transforms

}  // namespace tket
