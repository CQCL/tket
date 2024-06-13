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
#include "tket/OpType/OpTypeInfo.hpp"
#include "tket/Ops/Op.hpp"
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

static void apply_tqe_to_tableau(const TQE& tqe, UnitaryRevTableau& tab) {
  auto [gate_type, a_int, b_int] = tqe;
  Qubit a(a_int);
  Qubit b(b_int);
  switch (gate_type) {
    case TQEType::XX:
      tab.apply_gate_at_end(OpType::H, {a});
      tab.apply_gate_at_end(OpType::CX, {a, b});
      tab.apply_gate_at_end(OpType::H, {a});
      break;
    case TQEType::XY:
      tab.apply_gate_at_end(OpType::H, {a});
      tab.apply_gate_at_end(OpType::CY, {a, b});
      tab.apply_gate_at_end(OpType::H, {a});
      break;
    case TQEType::XZ:
      tab.apply_gate_at_end(OpType::CX, {b, a});
      break;
    case TQEType::YX:
      tab.apply_gate_at_end(OpType::H, {b});
      tab.apply_gate_at_end(OpType::CY, {b, a});
      tab.apply_gate_at_end(OpType::H, {b});
      break;
    case TQEType::YY:
      tab.apply_gate_at_end(OpType::V, {a});
      tab.apply_gate_at_end(OpType::CY, {a, b});
      tab.apply_gate_at_end(OpType::Vdg, {a});
      break;
    case TQEType::YZ:
      tab.apply_gate_at_end(OpType::CY, {b, a});
      break;
    case TQEType::ZX:
      tab.apply_gate_at_end(OpType::CX, {a, b});
      break;
    case TQEType::ZY:
      tab.apply_gate_at_end(OpType::CY, {a, b});
      break;
    case TQEType::ZZ:
      tab.apply_gate_at_end(OpType::CZ, {a, b});
      break;
  }
}

PauliExpNode::PauliExpNode(std::vector<unsigned> support_vec, Expr theta)
    : support_vec_(support_vec), theta_(theta) {
  tqe_cost_ = support_vec_.size() -
              std::count(support_vec_.begin(), support_vec_.end(), 0) - 1;
}

pauli_letter_distances_t PauliExpNode::all_distances(
    const std::vector<unsigned>& support,
    std::shared_ptr<Architecture> architecture,
    const std::map<unsigned, Node>& node_mapping) const {
  pauli_letter_distances_t letter_distances(
      architecture->get_diameter() + 1, 0);
  for (unsigned i = 0; i < support.size() - 1; i++) {
    for (unsigned j = i + 1; j < support.size(); j++) {
      // pauli strings are detailed as unsigned ints where 0 => Identity
      auto it = node_mapping.find(i);
      TKET_ASSERT(it != node_mapping.end());
      auto jt = node_mapping.find(j);
      TKET_ASSERT(jt != node_mapping.end());
      if (support[i] > 0 && support[j] > 0)
        letter_distances[architecture->get_distance(it->second, jt->second)] +=
            1;
    }
  }
  return letter_distances;
}

int PauliExpNode::tqe_cost_increase(const TQE& tqe) const {
  unsigned supp0 = support_vec_[std::get<1>(tqe)];
  unsigned supp1 = support_vec_[std::get<2>(tqe)];
  unsigned new_supp0, new_supp1;
  std::tie(new_supp0, new_supp1) =
      SINGLET_PAIR_TRANSFORMATION_MAP.at({std::get<0>(tqe), supp0, supp1});

  return (new_supp0 > 0) + (new_supp1 > 0) - (supp0 > 0) - (supp1 > 0);
}

int PauliExpNode::aas_tqe_cost_increase(
    const TQE& tqe, std::shared_ptr<Architecture> architecture,
    const std::map<unsigned, Node>& node_mapping) const {
  std::vector<unsigned> comparison = support_vec_;
  unsigned index0 = std::get<1>(tqe);
  unsigned index1 = std::get<2>(tqe);
  TKET_ASSERT(index0 < support_vec_.size());
  TKET_ASSERT(index1 < support_vec_.size());
  unsigned supp0 = support_vec_[index0];
  unsigned supp1 = support_vec_[index1];
  unsigned new_supp0, new_supp1;
  std::tie(new_supp0, new_supp1) =
      SINGLET_PAIR_TRANSFORMATION_MAP.at({std::get<0>(tqe), supp0, supp1});
  comparison[index0] = new_supp0;
  comparison[index1] = new_supp1;

  // how do we put a number to this? minimum is better so can start easy
  // first get distances
  pauli_letter_distances_t old_distances =
      this->all_distances(support_vec_, architecture, node_mapping);
  pauli_letter_distances_t new_distances =
      this->all_distances(comparison, architecture, node_mapping);

  // for distance d, if old_distances[d] - new_distances[d] < 0, then that entry
  // has increased given this, increased distances at larger d add a larger
  // contribution & vice
  int cost = 0;
  TKET_ASSERT(old_distances.size() == new_distances.size());
  for (unsigned i = 0; i < old_distances.size(); i++) {
    cost += (i * (old_distances[i] - new_distances[i]));
  }
  return cost;
}

void PauliExpNode::update(const TQE& tqe) {
  unsigned a = std::get<1>(tqe);
  unsigned b = std::get<2>(tqe);
  unsigned supp0 = support_vec_[a];
  unsigned supp1 = support_vec_[b];
  unsigned new_supp0, new_supp1;
  std::tie(new_supp0, new_supp1) =
      SINGLET_PAIR_TRANSFORMATION_MAP.at({std::get<0>(tqe), supp0, supp1});
  support_vec_[a] = new_supp0;
  support_vec_[b] = new_supp1;
  tqe_cost_ += (new_supp0 > 0) + (new_supp1 > 0) - (supp0 > 0) - (supp1 > 0);
}

std::vector<TQE> PauliExpNode::reduction_tqes() const {
  std::vector<TQE> tqes;
  // qubits with support
  std::vector<unsigned> sqs;
  for (unsigned i = 0; i < support_vec_.size(); i++) {
    if (support_vec_[i] > 0) sqs.push_back(i);
  }
  for (unsigned i = 0; i < sqs.size() - 1; i++) {
    for (unsigned j = i + 1; j < sqs.size(); j++) {
      std::vector<TQEType> tqe_types = SINGLET_PAIR_REDUCTION_TQES.at(
          {support_vec_[sqs[i]], support_vec_[sqs[j]]});
      for (const TQEType& tt : tqe_types) {
        tqes.push_back({tt, sqs[i], sqs[j]});
      }
    }
  }
  return tqes;
}

std::vector<TQE> PauliExpNode::reduction_tqes_all_letters(
    std::shared_ptr<Architecture> architecture,
    const std::map<unsigned, Node>& node_mapping) const {
  std::vector<TQE> tqes;
  // qubits with support
  std::vector<unsigned> sqs;
  // First we try to find Architecture permitted options that
  // will convert a Pauli letter to an Identity
  for (unsigned i = 0; i < support_vec_.size(); i++) {
    if (support_vec_[i] > 0) sqs.push_back(i);
  }
  for (unsigned i = 0; i < sqs.size(); i++) {
    for (unsigned j = 0; j < sqs.size(); j++) {
      if (i == j) continue;
      unsigned index_i = sqs[i];
      unsigned index_j = sqs[j];
      auto it = node_mapping.find(index_i);
      auto jt = node_mapping.find(index_j);
      TKET_ASSERT(it != node_mapping.end());
      TKET_ASSERT(jt != node_mapping.end());
      Node node_i = it->second;
      Node node_j = jt->second;
      if (architecture->edge_exists(node_i, node_j) ||
          architecture->edge_exists(node_j, node_i)) {
        std::cout << "permitted indices: " << sqs[i] << " " << sqs[j]
                  << std::endl;
        std::vector<TQEType> tqe_types = ALL_SINGLET_PAIR_REDUCTION_TQES.at(
            {support_vec_[sqs[i]], support_vec_[sqs[j]]});
        for (const TQEType& tt : tqe_types) {
          tqes.push_back({tt, sqs[i], sqs[j]});
        }
      }
    }
  }
  if (!tqes.empty()) {
    return tqes;
  }
  for (unsigned i = 0; i < support_vec_.size() - 1; i++) {
    for (unsigned j = 0; j < support_vec_.size(); j++) {
      if (i == j) continue;
      if (!(support_vec_[i] > 0 || support_vec_[j] > 0)) continue;
      auto it = node_mapping.find(i);
      auto jt = node_mapping.find(j);
      TKET_ASSERT(it != node_mapping.end());
      TKET_ASSERT(jt != node_mapping.end());
      Node node_i = it->second;
      Node node_j = jt->second;
      if (architecture->edge_exists(node_i, node_j) ||
          architecture->edge_exists(node_j, node_i)) {
        std::cout << "Edge on following is permitted: " << support_vec_[i]
                  << " " << support_vec_[j] << " " << i << " " << j
                  << std::endl;
        std::vector<TQEType> tqe_types = ALL_SINGLET_PAIR_REDUCTION_TQES.at(
            {support_vec_[i], support_vec_[j]});
        TKET_ASSERT(tqe_types.size() > 0);
        for (const TQEType& tt : tqe_types) {
          std::cout << i << " " << sqs.size() << " " << sqs[i] << " " << sqs[j] << std::endl;
          tqes.push_back({tt, i, j});
        }
      }
    }
  }
  TKET_ASSERT(!tqes.empty());
  return tqes;
}

std::pair<unsigned, unsigned> PauliExpNode::first_support() const {
  for (unsigned i = 0; i < support_vec_.size(); i++) {
    if (support_vec_[i] > 0) {
      return {i, support_vec_[i]};
    }
  }
  // Should be impossible to reach here
  TKET_ASSERT(false);
}

void PauliExpNode::pad_support_vector(unsigned width) {
  while (support_vec_.size() < width) {
    support_vec_.push_back(0);
  }
}

TableauRowNode::TableauRowNode(std::vector<unsigned> support_vec)
    : support_vec_(support_vec) {
  n_weaks_ = 0;
  n_strongs_ = 0;
  for (const unsigned& supp : support_vec_) {
    SupportType st = FACTOR_WEAKNESS_MAP.at(supp);
    if (st == SupportType::Strong) {
      n_strongs_++;
    } else if (st == SupportType::Weak) {
      n_weaks_++;
    }
  }
  tqe_cost_ = static_cast<unsigned>(1.5 * (n_strongs_ - 1) + n_weaks_);
}

int TableauRowNode::tqe_cost_increase(const TQE& tqe) const {
  unsigned supp0 = support_vec_[std::get<1>(tqe)];
  unsigned supp1 = support_vec_[std::get<2>(tqe)];
  unsigned new_supp0, new_supp1;
  std::tie(new_supp0, new_supp1) =
      FACTOR_PAIR_TRANSFORMATION_MAP.at({std::get<0>(tqe), supp0, supp1});
  SupportType st_supp0 = FACTOR_WEAKNESS_MAP.at(supp0);
  SupportType st_supp1 = FACTOR_WEAKNESS_MAP.at(supp1);
  SupportType st_new_supp0 = FACTOR_WEAKNESS_MAP.at(new_supp0);
  SupportType st_new_supp1 = FACTOR_WEAKNESS_MAP.at(new_supp1);
  unsigned old_strongs =
      (st_supp0 == SupportType::Strong) + (st_supp1 == SupportType::Strong);
  unsigned old_weaks =
      (st_supp0 == SupportType::Weak) + (st_supp1 == SupportType::Weak);
  unsigned new_strongs = (st_new_supp0 == SupportType::Strong) +
                         (st_new_supp1 == SupportType::Strong);
  unsigned new_weaks =
      (st_new_supp0 == SupportType::Weak) + (st_new_supp1 == SupportType::Weak);
  int strong_increase = new_strongs - old_strongs;
  int weak_increase = new_weaks - old_weaks;
  return static_cast<int>(1.5 * strong_increase + weak_increase);
}

void TableauRowNode::update(const TQE& tqe) {
  unsigned a = std::get<1>(tqe);
  unsigned b = std::get<2>(tqe);
  unsigned supp0 = support_vec_[a];
  unsigned supp1 = support_vec_[b];
  unsigned new_supp0, new_supp1;
  std::tie(new_supp0, new_supp1) =
      FACTOR_PAIR_TRANSFORMATION_MAP.at({std::get<0>(tqe), supp0, supp1});
  support_vec_[a] = new_supp0;
  support_vec_[b] = new_supp1;
  SupportType st_supp0 = FACTOR_WEAKNESS_MAP.at(supp0);
  SupportType st_supp1 = FACTOR_WEAKNESS_MAP.at(supp1);
  SupportType st_new_supp0 = FACTOR_WEAKNESS_MAP.at(new_supp0);
  SupportType st_new_supp1 = FACTOR_WEAKNESS_MAP.at(new_supp1);
  unsigned old_strongs =
      (st_supp0 == SupportType::Strong) + (st_supp1 == SupportType::Strong);
  unsigned old_weaks =
      (st_supp0 == SupportType::Weak) + (st_supp1 == SupportType::Weak);
  unsigned new_strongs = (st_new_supp0 == SupportType::Strong) +
                         (st_new_supp1 == SupportType::Strong);
  unsigned new_weaks =
      (st_new_supp0 == SupportType::Weak) + (st_new_supp1 == SupportType::Weak);
  n_strongs_ += new_strongs - old_strongs;
  n_weaks_ += new_weaks - old_weaks;
  tqe_cost_ = static_cast<unsigned>(1.5 * (n_strongs_ - 1) + n_weaks_);
}

std::vector<TQE> TableauRowNode::reduction_tqes() const {
  std::vector<TQE> tqes;
  // qubits with support
  std::vector<unsigned> sqs;
  for (unsigned i = 0; i < support_vec_.size(); i++) {
    if (support_vec_[i] > 0) sqs.push_back(i);
  }
  for (unsigned i = 0; i < sqs.size() - 1; i++) {
    for (unsigned j = i + 1; j < sqs.size(); j++) {
      std::vector<TQEType> tqe_types;
      unsigned a = sqs[i];
      unsigned b = sqs[j];
      unsigned supp0 = support_vec_[a];
      unsigned supp1 = support_vec_[b];
      SupportType st_supp0 = FACTOR_WEAKNESS_MAP.at(supp0);
      SupportType st_supp1 = FACTOR_WEAKNESS_MAP.at(supp1);
      if (st_supp0 == SupportType::Strong) {
        if (st_supp1 == SupportType::Strong) {
          // TQEs that transform a SS pair to WW
          tqe_types = FACTOR_PAIR_SS_TO_WW_TQES.at({supp0, supp1});
        } else {
          // TQEs that transform a SW pair to a single strong
          tqe_types = FACTOR_PAIR_SW_TO_SN_TQES.at({supp0, supp1});
        }
      } else {
        if (st_supp1 == SupportType::Strong) {
          // TQEs that transform a WS pair to a single strong
          tqe_types = FACTOR_PAIR_SW_TO_SN_TQES.at({supp1, supp0});
          // flip qubits
          a = sqs[j];
          b = sqs[i];
        } else {
          // TQEs that transform a WW pair to a single weak, not always
          // possible
          tqe_types = FACTOR_PAIR_WW_TO_WN_OR_NW_TQES.at({supp0, supp1});
        }
      }
      for (const TQEType& tt : tqe_types) {
        tqes.push_back({tt, a, b});
      }
    }
  }
  return tqes;
}

std::pair<unsigned, unsigned> TableauRowNode::first_support() const {
  for (unsigned i = 0; i < support_vec_.size(); i++) {
    if (support_vec_[i] > 0) {
      return {i, support_vec_[i]};
    }
  }
  // Should be impossible to reach here
  TKET_ASSERT(false);
}

// return the sum of the cost increases on remaining tableau nodes
static double default_tableau_tqe_cost(
    const std::vector<TableauRowNode>& rows,
    const std::vector<unsigned>& remaining_indices, const TQE& tqe) {
  double cost = 0;
  for (const unsigned& index : remaining_indices) {
    cost += rows[index].tqe_cost_increase(tqe);
  }
  return cost;
}

// return the weighted sum of the cost increases on remaining nodes
// we discount the weight after each set
static double default_pauliexp_tqe_cost(
    const double discount_rate,
    const std::vector<std::vector<PauliExpNode>>& rotation_sets,
    const std::vector<TableauRowNode>& rows, const TQE& tqe) {
  double discount = 1 / (1 + discount_rate);
  double weight = 1;
  double exp_cost = 0;
  double tab_cost = 0;
  for (const std::vector<PauliExpNode>& rotation_set : rotation_sets) {
    for (const PauliExpNode& node : rotation_set) {
      exp_cost += weight * node.tqe_cost_increase(tqe);
    }
    weight *= discount;
  }
  for (const TableauRowNode& node : rows) {
    tab_cost += weight * node.tqe_cost_increase(tqe);
  }
  return exp_cost + tab_cost;
}

// return the weighted sum of the cost increases on remaining nodes
// we discount the weight after each set
static double aas_pauliexp_tqe_cost(
    const double discount_rate,
    const std::vector<std::vector<PauliExpNode>>& rotation_sets,
    const std::vector<TableauRowNode>& rows, const TQE& tqe,
    std::shared_ptr<Architecture> architecture,
    const std::map<unsigned, Node>& node_mapping) {
  double discount = 1 / (1 + discount_rate);
  double weight = 1;
  double exp_cost = 0;
  double tab_cost = 0;
  for (const std::vector<PauliExpNode>& rotation_set : rotation_sets) {
    for (const PauliExpNode& node : rotation_set) {
      exp_cost +=
          weight * node.aas_tqe_cost_increase(tqe, architecture, node_mapping);
    }
    weight *= discount;
  }
  for (const TableauRowNode& node : rows) {
    tab_cost += weight * node.tqe_cost_increase(tqe);
  }
  return exp_cost + tab_cost;
}

// given a map from TQE to a vector of costs, select the one with the minimum
// weighted sum of minmax-normalised costs
static TQE minmax_selection(
    const std::map<TQE, std::vector<double>>& tqe_candidates_cost,
    const std::vector<double>& weights) {
  std::cout << "tqe candidates cost size: " << tqe_candidates_cost.size()
            << std::endl;
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
 * @brief Given a tableau that is identity up to local Cliffords, qubit
 * permutation, and signs, transform it to exact identity and adding gates to
 * a circuit
 */
static void tableau_cleanup(
    std::vector<TableauRowNode>& rows, UnitaryRevTableau& tab, Circuit& circ) {
  // apply local Cliffords
  for (const TableauRowNode& node : rows) {
    unsigned q_index, supp;
    std::tie(q_index, supp) = node.first_support();
    Qubit q(q_index);
    std::vector<LocalCliffordType> local_cliffords =
        FACTOR_STRONG_TO_LOCALS.at(supp);
    for (const LocalCliffordType& lc : local_cliffords) {
      switch (lc) {
        case LocalCliffordType::H:
          tab.apply_gate_at_end(OpType::H, {q});
          circ.add_op<UnitID>(OpType::H, {q});
          break;
        case LocalCliffordType::S:
          tab.apply_gate_at_end(OpType::S, {q});
          circ.add_op<UnitID>(OpType::S, {q});
          break;
        case LocalCliffordType::V:
          tab.apply_gate_at_end(OpType::V, {q});
          circ.add_op<UnitID>(OpType::V, {q});
          break;
      }
    }
  }
  // remove signs
  for (const Qubit& q : circ.all_qubits()) {
    if (cast_coeff<quarter_turns_t, Complex>(tab.get_xrow(q).coeff) != 1.) {
      tab.apply_gate_at_end(OpType::Z, {q});
      circ.add_op<UnitID>(OpType::Z, {q});
    }
    if (cast_coeff<quarter_turns_t, Complex>(tab.get_zrow(q).coeff) != 1.) {
      tab.apply_gate_at_end(OpType::X, {q});
      circ.add_op<UnitID>(OpType::X, {q});
    }
  }
  // remove permutations
  // 1. find perm
  unsigned n_qubits = circ.n_qubits();
  std::vector<unsigned> perm(n_qubits);
  for (unsigned i = 0; i < n_qubits; i++) {
    QubitPauliMap z_row_string = tab.get_zrow(Qubit(i)).string;
    for (auto it = z_row_string.begin(); it != z_row_string.end(); it++) {
      if (it->second == Pauli::Z) {
        perm[it->first.index()[0]] = i;
        break;
      }
    }
  }
  // 2. traverse transpositions
  std::unordered_set<unsigned> done;
  for (unsigned k = 0; k < n_qubits; k++) {
    if (done.find(k) != done.end()) {
      continue;
    }
    unsigned head = k;
    unsigned current = k;
    unsigned next = perm[k];
    while (true) {
      if (next == head) {
        done.insert(current);
        break;
      }
      // the SWAP gates will be later converted to wire swaps
      tab.apply_gate_at_end(OpType::SWAP, {Qubit(current), Qubit(next)});
      circ.add_op<unsigned>(OpType::SWAP, {current, next});
      done.insert(current);
      current = next;
      next = perm[current];
    }
  }
}

/**
 * @brief Synthesise a vector of TableauRowNode
 */
static void tableau_row_nodes_synthesis(
    std::vector<TableauRowNode>& rows, UnitaryRevTableau& tab, Circuit& circ,
    double depth_weight, DepthTracker& depth_tracker) {
  // only consider nodes with a non-zero cost
  std::vector<unsigned> remaining_indices;
  for (unsigned i = 0; i < rows.size(); i++) {
    if (rows[i].tqe_cost() > 0) {
      remaining_indices.push_back(i);
    }
  }
  while (remaining_indices.size() != 0) {
    // get nodes with min cost
    std::vector<unsigned> min_nodes_indices = {remaining_indices[0]};
    unsigned min_cost = rows[remaining_indices[0]].tqe_cost();
    for (unsigned i = 1; i < remaining_indices.size(); i++) {
      unsigned node_cost = rows[remaining_indices[i]].tqe_cost();
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
      std::vector<TQE> node_reducing_tqes = rows[index].reduction_tqes();
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
    apply_tqe_to_tableau(selected_tqe, tab);
    // update depth tracker
    depth_tracker.add_2q_gate(
        std::get<1>(selected_tqe), std::get<2>(selected_tqe));
    // remove finished nodes
    for (unsigned i = remaining_indices.size(); i-- > 0;) {
      unsigned node_index = remaining_indices[i];
      rows[node_index].update(selected_tqe);
      if (rows[node_index].tqe_cost() == 0) {
        remaining_indices.erase(remaining_indices.begin() + i);
      }
    }
  }
  tableau_cleanup(rows, tab, circ);
}

/**
 * @brief Given a vector of sets of PauliExpNode, implement any node in the
 * first set where the tqe_cost is zero. Remove implemented nodes and the first
 * set if empty.
 *
 * @param rotation_sets
 * @param tab
 * @param circ
 * @return true if the first set is now empty and removed
 * @return false
 */
static bool consume_available_rotations(
    std::vector<std::vector<PauliExpNode>>& rotation_sets,
    UnitaryRevTableau& tab, Circuit& circ, DepthTracker& depth_tracker) {
  std::vector<unsigned> bin;
  if (rotation_sets.size() == 0) {
    return false;
  }
  std::vector<PauliExpNode>& first_set = rotation_sets[0];
  for (unsigned i = 0; i < first_set.size(); i++) {
    PauliExpNode& node = first_set[i];
    if (node.tqe_cost() > 0) continue;
    unsigned q_index, supp;
    std::tie(q_index, supp) = node.first_support();
    Qubit q(q_index);
    depth_tracker.add_1q_gate(q_index);
    switch (supp) {
      case 3: {
        // we apply S gate only to the frame, then check the sign, then Sdg
        // if + apply f.Sdg; circ.Ry(-a)
        // if - apply f.Sdg; circ.Ry(a)
        tab.apply_gate_at_end(OpType::S, {q});
        Complex x_coeff =
            cast_coeff<quarter_turns_t, Complex>(tab.get_xrow(q).coeff);
        tab.apply_gate_at_end(OpType::Sdg, {q});
        if (x_coeff == 1.) {
          circ.add_op<UnitID>(OpType::Ry, -node.theta(), {q});
        } else {
          circ.add_op<UnitID>(OpType::Ry, node.theta(), {q});
        }
        break;
      }
      case 1: {
        Complex z_coeff =
            cast_coeff<quarter_turns_t, Complex>(tab.get_zrow(q).coeff);
        if (z_coeff == 1.) {
          circ.add_op<UnitID>(OpType::Rz, node.theta(), {q});
        } else {
          circ.add_op<UnitID>(OpType::Rz, -node.theta(), {q});
        }
        break;
      }
      case 2: {
        Complex x_coeff =
            cast_coeff<quarter_turns_t, Complex>(tab.get_xrow(q).coeff);
        if (x_coeff == 1.) {
          circ.add_op<UnitID>(OpType::Rx, node.theta(), {q});
        } else {
          circ.add_op<UnitID>(OpType::Rx, -node.theta(), {q});
        }
        break;
      }
      default:
        // support can't be Pauli::I
        TKET_ASSERT(false);
    }
    bin.push_back(i);
  }
  if (bin.size() == 0) return false;
  // sort the bin so we remove elements from back to front
  std::sort(bin.begin(), bin.end(), std::greater<unsigned>());
  for (const unsigned& index : bin) {
    first_set.erase(first_set.begin() + index);
  }
  if (first_set.size() == 0) {
    rotation_sets.erase(rotation_sets.begin());
    return true;
  }
  return false;
}

/**
 * @brief Synthesise a vector of unordered rotation sets
 */
static void pauli_exps_synthesis(
    std::vector<std::vector<PauliExpNode>>& rotation_sets,
    std::vector<TableauRowNode>& rows, UnitaryRevTableau& tab, Circuit& circ,
    double discount_rate, double depth_weight, DepthTracker& depth_tracker) {
  while (true) {
    while (consume_available_rotations(
        rotation_sets, tab, circ, depth_tracker));  // do nothing
    if (rotation_sets.size() == 0) break;
    std::vector<PauliExpNode>& first_set = rotation_sets[0];
    // get nodes with min cost
    std::vector<unsigned> min_nodes_indices = {0};
    unsigned min_cost = first_set[0].tqe_cost();
    for (unsigned i = 1; i < first_set.size(); i++) {
      unsigned node_cost = first_set[i].tqe_cost();
      if (node_cost == min_cost) {
        min_nodes_indices.push_back(i);
      } else if (node_cost < min_cost) {
        min_nodes_indices = {i};
        min_cost = node_cost;
      }
    }
    std::set<TQE> tqe_candidates;
    for (const unsigned& index : min_nodes_indices) {
      std::vector<TQE> node_reducing_tqes = first_set[index].reduction_tqes();
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
    apply_tqe_to_tableau(selected_tqe, tab);
    depth_tracker.add_2q_gate(
        std::get<1>(selected_tqe), std::get<2>(selected_tqe));
    for (std::vector<PauliExpNode>& rotation_set : rotation_sets) {
      for (PauliExpNode& node : rotation_set) {
        node.update(selected_tqe);
      }
    }
    for (TableauRowNode& row : rows) {
      row.update(selected_tqe);
    }
  }
}

/**
 * @brief Synthesise a vector of unordered rotation sets
 */
static void aas_pauli_exps_synthesis(
    std::vector<std::vector<PauliExpNode>>& rotation_sets,
    std::vector<TableauRowNode>& rows, UnitaryRevTableau& tab, Circuit& circ,
    double discount_rate, double depth_weight, DepthTracker& depth_tracker,
    std::shared_ptr<Architecture> architecture,
    const std::map<unsigned, Node>& node_mapping) {
  while (true) {
    while (consume_available_rotations(
        rotation_sets, tab, circ, depth_tracker));  // do nothing
    if (rotation_sets.size() == 0) break;
    std::vector<PauliExpNode>& first_set = rotation_sets[0];
    // get nodes with min cost
    std::vector<unsigned> min_nodes_indices = {0};
    unsigned min_cost = first_set[0].tqe_cost();
    for (unsigned i = 1; i < first_set.size(); i++) {
      unsigned node_cost = first_set[i].tqe_cost();
      if (node_cost == min_cost) {
        min_nodes_indices.push_back(i);
      } else if (node_cost < min_cost) {
        min_nodes_indices = {i};
        min_cost = node_cost;
      }
    }
    std::set<TQE> tqe_candidates;
    for (const unsigned& index : min_nodes_indices) {
      std::vector<TQE> node_reducing_tqes =
          first_set[index].reduction_tqes_all_letters(
              architecture, node_mapping);
      std::cout << "n of cool new candidates produced: " << node_reducing_tqes.size() << std::endl;
      tqe_candidates.insert(
          node_reducing_tqes.begin(), node_reducing_tqes.end());
    }

    std::cout << "Number of candidates produced: " << tqe_candidates.size()
              << " | number of remaining rotation sets " << rotation_sets.size()
              << " | number of remaining exponetials in set: "
              << first_set.size() << std::endl;

    // for each tqe we compute costs which might subject to normalisation
    std::map<TQE, std::vector<double>> tqe_candidates_cost;
    for (const TQE& tqe : tqe_candidates) {
      // TODO: this assumes bidirectional edges, is this ok?
      auto it = node_mapping.find(std::get<1>(tqe));
      TKET_ASSERT(it != node_mapping.end());
      auto jt = node_mapping.find(std::get<2>(tqe));
      TKET_ASSERT(jt != node_mapping.end());
      if (architecture->edge_exists(it->second, jt->second)) {
        tqe_candidates_cost.insert(
            {tqe,
             {aas_pauliexp_tqe_cost(
                  discount_rate, rotation_sets, rows, tqe, architecture,
                  node_mapping),
              // static_cast<double>(depth_tracker.gate_depth(
              //     std::get<1>(tqe), std::get<2>(tqe)))}});
              static_cast<double>(1)}});
      }
    }
    // select the best one
    std::cout << "All costs: " << std::endl;
    for (auto c : tqe_candidates_cost) {
      std::cout << std::get<1>(c.first) << " " << std::get<2>(c.first) << " "
                << c.second[0] << " " << c.second[1] << std::endl;
    }
    TQE selected_tqe = select_pauliexp_tqe(tqe_candidates_cost, depth_weight);

    // std::cout << "What did it pick? " << std::get<1>(selected_tqe) << " " <<
    // std::get<2>(selected_tqe) << std::endl;

    // apply TQE
    apply_tqe_to_circ(selected_tqe, circ);
    apply_tqe_to_tableau(selected_tqe, tab);
    depth_tracker.add_2q_gate(
        std::get<1>(selected_tqe), std::get<2>(selected_tqe));
    for (std::vector<PauliExpNode>& rotation_set : rotation_sets) {
      for (PauliExpNode& node : rotation_set) {
        node.update(selected_tqe);
      }
    }
    for (TableauRowNode& row : rows) {
      row.update(selected_tqe);
    }
  }
}

// convert a Pauli exponential to a PauliExpNode
static PauliExpNode get_node_from_exp(
    const std::vector<Pauli>& paulis, const Expr& theta,
    const qubit_vector_t& args, unsigned n, const UnitaryTableau& forward_tab,
    const UnitaryRevTableau& tab) {
  std::map<Qubit, Pauli> pauli_map;
  for (unsigned i = 0; i < args.size(); i++) {
    pauli_map.insert({args[i], paulis[i]});
  }
  // this has the effect of bringing the final clifford
  // forward past the Pauli exponential
  SpPauliStabiliser pstab =
      forward_tab.get_row_product(SpPauliStabiliser(pauli_map));
  Complex sign = cast_coeff<quarter_turns_t, Complex>(pstab.coeff);

  std::vector<unsigned> support_vec;
  for (unsigned i = 0; i < n; i++) {
    SpPauliStabiliser zrow = tab.get_zrow(Qubit(i));
    SpPauliStabiliser xrow = tab.get_xrow(Qubit(i));
    bool z_supp = !zrow.commutes_with(pstab);
    bool x_supp = !xrow.commutes_with(pstab);
    if (!z_supp && !x_supp) {
      support_vec.push_back(0);
    } else if (!z_supp && x_supp) {
      support_vec.push_back(1);
    } else if (z_supp && !x_supp) {
      support_vec.push_back(2);
    } else if (z_supp && x_supp) {
      support_vec.push_back(3);
    }
  }
  return PauliExpNode(support_vec, sign.real() * theta);
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

Circuit greedy_pauli_graph_synthesis(
    const Circuit& circ, double discount_rate, double depth_weight,
    std::optional<std::shared_ptr<Architecture>> architecture) {
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
  std::vector<std::vector<PauliExpNode>> rotation_sets;
  std::vector<TableauRowNode> rows;
  // use forward Tableau to update the paulis by commuting the tableau to the
  // front
  UnitaryTableau forward_tab = circuit_to_unitary_tableau(cliff);
  // Tableau used for tracking Cliffords throughout the synthesis
  // TODO: this can be potentially made redundant
  UnitaryRevTableau tab = circuit_to_unitary_rev_tableau(cliff).dagger();
  unsigned n_qubits = c.n_qubits();
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
          rotation_sets.push_back({get_node_from_exp(
              paulis, phase, cmd.get_qubits(), n_qubits, forward_tab, tab)});
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
        std::vector<PauliExpNode> rotation_set;
        if (trivial1) {
          c.add_phase(global_phase1);
        } else {
          rotation_set.push_back(get_node_from_exp(
              paulis1, phase1, cmd.get_qubits(), n_qubits, forward_tab, tab));
        }
        if (trivial2) {
          c.add_phase(global_phase2);
        } else {
          rotation_set.push_back(get_node_from_exp(
              paulis2, phase2, cmd.get_qubits(), n_qubits, forward_tab, tab));
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
        std::vector<PauliExpNode> rotation_set;
        for (const SymPauliTensor& pt : gadgets) {
          const std::vector<Pauli> paulis = pt.string;
          const Expr phase = pt.coeff;
          auto [trivial, global_phase] = is_trivial_pauliexp(paulis, phase);
          if (trivial) {
            c.add_phase(global_phase);
          } else {
            rotation_set.push_back(get_node_from_exp(
                paulis, phase, cmd.get_qubits(), n_qubits, forward_tab, tab));
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
  // add identity TableauRowNodes
  // the tableau width is either set to the number of qubits or the number of
  // nodes depending on whether the architecture is present
  unsigned tableau_width = n_qubits;
  if (architecture) {
    tableau_width = (*architecture)->n_nodes();
  }
  for (unsigned i = 0; i < tableau_width; i++) {
    std::vector<unsigned> support_vec;
    // identity rows
    std::map<Qubit, Pauli> p;
    std::map<Qubit, Pauli> q;
    for (unsigned j = 0; j < tableau_width; j++) {
      if (j == i) {
        p.insert({Qubit(j), Pauli::Z});
        q.insert({Qubit(j), Pauli::X});
      } else {
        p.insert({Qubit(j), Pauli::I});
        q.insert({Qubit(j), Pauli::I});
      }
    }
    SpPauliStabiliser stab_p(p);
    SpPauliStabiliser stab_q(q);
    for (unsigned row_index = 0; row_index < tableau_width; row_index++) {
      SpPauliStabiliser zrow = tab.get_zrow(Qubit(row_index));
      SpPauliStabiliser xrow = tab.get_xrow(Qubit(row_index));
      bool lpx = !xrow.commutes_with(stab_p);
      bool lpz = !zrow.commutes_with(stab_p);
      bool lqx = !xrow.commutes_with(stab_q);
      bool lqz = !zrow.commutes_with(stab_q);
      support_vec.push_back(8 * lpx + 4 * lpz + 2 * lqx + lqz);
    }
    rows.push_back(TableauRowNode(support_vec));
  }
  DepthTracker depth_tracker(n_qubits);

  // here we also concern ourselves with architecture
  // essentially, if architecture contains an architecture we use it, else we
  // don't ...
  if (architecture) {
    // synthesise Pauli exps
    // TODO: Replace naive node mapping with something from graph placement ...
    std::map<unsigned, Node> node_mapping;
    std::shared_ptr<Architecture> arc = *architecture;
    TKET_ASSERT(arc->n_nodes() >= circ.n_qubits());
    std::vector<Node> nodes = arc->get_all_nodes_vec();
    for (unsigned i = 0; i < nodes.size(); i++) {
      node_mapping.insert({i, nodes[i]});
    }
    for (std::vector<PauliExpNode>& pexns : rotation_sets) {
      for (PauliExpNode& pexn : pexns) {
        pexn.pad_support_vector(nodes.size());
      }
    }
    aas_pauli_exps_synthesis(
        rotation_sets, rows, tab, c, discount_rate, depth_weight, depth_tracker,
        arc, node_mapping);
  } else {
    // synthesise Pauli exps
    pauli_exps_synthesis(
        rotation_sets, rows, tab, c, discount_rate, depth_weight,
        depth_tracker);
  }
  // synthesise the tableau
  tableau_row_nodes_synthesis(rows, tab, c, depth_weight, depth_tracker);
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
        circ, discount_rate, depth_weight, std::nullopt);
    singleq_clifford_sweep().apply(circ);
    return true;
  });
}

Transform aas_greedy_pauli_optimisation(
    std::shared_ptr<Architecture> architecture, double discount_rate,
    double depth_weight) {
  return Transform([architecture, discount_rate, depth_weight](Circuit& circ) {
    synthesise_pauli_graph(PauliSynthStrat::Sets, CXConfigType::Snake)
        .apply(circ);
    circ = GreedyPauliSimp::greedy_pauli_graph_synthesis(
        circ, discount_rate, depth_weight, architecture);
    singleq_clifford_sweep().apply(circ);
    return true;
  });
}

}  // namespace Transforms

}  // namespace tket
