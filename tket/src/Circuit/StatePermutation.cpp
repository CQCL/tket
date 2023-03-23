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

#include "Circuit/StatePermutation.hpp"

#include <boost/graph/max_cardinality_matching.hpp>

#include "Circuit/Circuit.hpp"
#include "Circuit/DiagonalBox.hpp"
#include "Circuit/Multiplexor.hpp"
#include "Gate/Rotation.hpp"
#include "Ops/OpJsonFactory.hpp"
#include "Utils/HelperFunctions.hpp"
#include "Utils/Json.hpp"

namespace tket {

StatePermutationBox::StatePermutationBox(const state_perm_t &permutation)
    : Box(OpType::StatePermutationBox), permutation_(permutation) {
  // verify the permutation is valid
  // we need to check every element has the same size, and every bitstrings
  // appears twice
  if (permutation.size() == 0) {
    throw std::invalid_argument("The permutation is empty.");
  }
  auto it = permutation.begin();
  unsigned n_qubits = (unsigned)it->first.size();
  std::unordered_set<std::vector<bool>> rhs_sates;
  if (permutation.size() != 1 << n_qubits) {
    throw std::invalid_argument(
        "The permutation doesn't contain all bitstrings.");
  }
  for (; it != permutation.end(); ++it) {
    if (it->first.size() != n_qubits || it->second.size() != n_qubits) {
      throw std::invalid_argument("Bitstrings don't have the same size.");
    }
    if (rhs_sates.find(it->second) != rhs_sates.end()) {
      throw std::invalid_argument("The permutation contains duplicate values.");
    }
    rhs_sates.insert(it->second);
  }
}

StatePermutationBox::StatePermutationBox(const StatePermutationBox &other)
    : Box(other), permutation_(other.permutation_) {}

Op_ptr StatePermutationBox::dagger() const {
  state_perm_t reverse_perm;
  for (auto it = permutation_.begin(); it != permutation_.end(); it++) {
    reverse_perm.insert({it->second, it->first});
  }
  return std::make_shared<StatePermutationBox>(reverse_perm);
}

Op_ptr StatePermutationBox::transpose() const { return dagger(); }

op_signature_t StatePermutationBox::get_signature() const {
  op_signature_t qubits(
      (unsigned)permutation_.begin()->first.size(), EdgeType::Quantum);
  return qubits;
}

state_perm_t StatePermutationBox::get_permutation() const {
  return permutation_;
}

void print_op_map(const ctrl_op_map_t& op_map) {
		for (auto it = op_map.begin(); it != op_map.end(); it++){
			std::cout<<"(";
			for(auto b : it->first) {
				if (b) {
					std::cout<<"1,";
				}
				else {
					std::cout<<"0,";
				}
			}
			std::cout<<")\n";
		}
}

void print_bool_vec(const std::vector<bool>& bool_vec) {
	std::cout<<"(";
	for(auto b : bool_vec) {
		if (b) {
			std::cout<<"1,";
		}
		else {
			std::cout<<"0,";
		}
	}
	std::cout<<")\n";
}
// big endian
static unsigned bin_to_dec(const std::vector<bool> &bin) {
  unsigned res = 0;
  for (unsigned i = 0; i < bin.size(); i++) {
    if (bin[i]) {
      res = res + (1 << (bin.size() - 1 - i));
    }
  }
  return res;
}

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
    cube_graph_t;
typedef boost::graph_traits<cube_graph_t>::vertex_descriptor CubeVert;

static std::vector<std::vector<bool>> distribute_along_qubit(
    const std::vector<bool> &prefix, unsigned n_qubits, unsigned partition_q,
    state_perm_t &perm) {
  // number of free qubits on the subcube defined by prefix and
  partition_q = 0;
  unsigned n_free_qubits = n_qubits - prefix.size() - 1;
  // construct the bipartie graph
  cube_graph_t g(1 << (n_free_qubits + 1));
  for (unsigned i = 0; i < (1 << n_free_qubits); i++) {
    std::vector<bool> postfix = dec_to_bin(i, n_free_qubits);
    std::vector<bool> upper_v(prefix);
    upper_v.push_back(0);
    upper_v.insert(upper_v.end(), postfix.begin(), postfix.end());
    std::vector<bool> lower_v(prefix);
    lower_v.push_back(1);
    lower_v.insert(lower_v.end(), postfix.begin(), postfix.end());
    std::vector<bool> upper_state_postfix(
        perm[upper_v].begin() + prefix.size() + 1, perm[upper_v].end());
    std::vector<bool> lower_state_postfix(
        perm[lower_v].begin() + prefix.size() + 1, perm[lower_v].end());
    boost::add_edge(
        i, bin_to_dec(upper_state_postfix) + (1 << n_free_qubits), g);
    if (lower_state_postfix != upper_state_postfix) {
      boost::add_edge(
          i, bin_to_dec(lower_state_postfix) + (1 << n_free_qubits), g);
    }
  }
  std::vector<CubeVert> mate(1 << (n_free_qubits + 1));
  edmonds_maximum_cardinality_matching(g, &mate[0]);
  std::vector<std::vector<bool>> swap_edges;
  for (unsigned i = 0; i < (1 << n_free_qubits); i++) {
    std::vector<bool> postfix = dec_to_bin(i, n_free_qubits);
    std::vector<bool> upper_v(prefix);
    upper_v.push_back(0);
    upper_v.insert(upper_v.end(), postfix.begin(), postfix.end());
    std::vector<bool> mapped_postfix = dec_to_bin(mate[i] - (1 << n_free_qubits), n_free_qubits);
    if (!std::equal(
            perm[upper_v].begin() + prefix.size() + 1, perm[upper_v].end(),
            mapped_postfix.begin())) {
      upper_v.erase(upper_v.begin() + partition_q);
      swap_edges.push_back(upper_v);
    }
  }
  return swap_edges;
}

// swap two states along qubit q
// assume two states only differ at q
// assume zflip_op to be one of {Rx(1), Ry(1)}
// update op_map, phases, and permutation
static void swap_states(
    const std::vector<bool> &edge, unsigned diff_q, const Op_ptr &zflip_op,
    ctrl_op_map_t &op_map, std::vector<Complex> &phases, state_perm_t &perm) {
  std::vector<bool> a(edge);
  std::vector<bool> b(edge);
  a.insert(a.begin() + diff_q, 0);
  b.insert(b.begin() + diff_q, 1);
  std::swap(perm.at(a), perm.at(b));
  unsigned a_bin = bin_to_dec(a);
  unsigned b_bin = bin_to_dec(b);
  if (zflip_op->get_type() == OpType::Rx &&
      equiv_val(zflip_op->get_params()[0], 1., 4)) {
    phases[a_bin] *= -i_;
    phases[b_bin] *= -i_;
  } else if (
      zflip_op->get_type() == OpType::Ry &&
      equiv_val(zflip_op->get_params()[0], 1., 4)) {
    phases[b_bin] *= -1;
  } else {
    throw std::invalid_argument("Unsupported rotation.");
  }
  std::swap(phases[a_bin], phases[b_bin]);
  op_map.insert({edge, zflip_op});
}

/**
 * @brief Recursively route the state permutation using Multiplexed Rx/Ry gates
 *
 * @param circ the circuit to implement the permutation
 * @param prev_partition_q defines the parallel subcubes partitioned by
 * [0,...,prev_partition_q]
 * @param phases keep track of the phases induced by Rx/Ry
 * @param perm permutation
 * @param zflip_op the base 1-q rotation for swap two states. Currently limited
 * to Rx(1) and Ry(1)
 *
 */
static void route_recursive(
    Circuit &circ, std::optional<unsigned> prev_partition_q,
    std::vector<Complex> &phases, state_perm_t &perm, const Op_ptr &zflip_op) {
  unsigned partition_q;
  if (prev_partition_q == std::nullopt) {
    partition_q = 0;
  } else {
    partition_q = prev_partition_q.value() + 1;
  }

  unsigned n_qubits = circ.n_qubits();
  std::vector<unsigned> multplx_args(n_qubits);
  std::iota(std::begin(multplx_args), std::end(multplx_args), 0);
  multplx_args.erase(multplx_args.begin() + partition_q);
  multplx_args.push_back(partition_q);
  // base case handling
  if (partition_q == n_qubits - 1) {
    ctrl_op_map_t op_map;
    for (unsigned i = 0; i < (1 << (n_qubits - 1)); i++) {
      std::vector<bool> edge = dec_to_bin(i, n_qubits - 1);
      edge.push_back(0);
      if (perm[edge][partition_q] != 0) {
        edge.pop_back();
        swap_states(edge, partition_q, zflip_op, op_map, phases, perm);
      }
    }
    if (!op_map.empty()) {
      circ.add_box(MultiplexedRotationBox(op_map), multplx_args);
    }
    return;
  }
  // distribute states for each parallel cubes defined by prev_partition_q
  ctrl_op_map_t op_map;
  if (prev_partition_q == std::nullopt) {
    std::vector<std::vector<bool>> swap_edges =
        distribute_along_qubit({}, circ.n_qubits(), partition_q, perm);
    for (auto e : swap_edges) {
      swap_states(e, partition_q, zflip_op, op_map, phases, perm);
    }
  } else {
    for (unsigned i = 0; i < (1 << (prev_partition_q.value() + 1)); i++) {
      std::vector<std::vector<bool>> swap_edges = distribute_along_qubit(
          dec_to_bin(i, prev_partition_q.value() + 1), circ.n_qubits(),
          partition_q, perm);
      for (auto e : swap_edges) {
        swap_states(e, partition_q, zflip_op, op_map, phases, perm);
      }
    }
  }
  if (!op_map.empty()) {
    circ.add_box(MultiplexedRotationBox(op_map), multplx_args);
  }

  route_recursive(circ, partition_q, phases, perm, zflip_op);

  ctrl_op_map_t last_op_map;
  for (unsigned i = 0; i < (1 << (n_qubits - 1)); i++) {
    std::vector<bool> edge = dec_to_bin(i, n_qubits - 1);
    edge.insert(edge.begin() + partition_q, 0);
    if (perm[edge][partition_q] != 0) {
      edge.erase(edge.begin() + partition_q);
      swap_states(edge, partition_q, zflip_op, last_op_map, phases, perm);
    }
  }
  if (!last_op_map.empty()) {
    circ.add_box(MultiplexedRotationBox(last_op_map), multplx_args);
  }
};

void StatePermutationBox::generate_circuit() const {
  unsigned n_qubits = (unsigned)permutation_.begin()->first.size();
  Circuit circ(n_qubits);
  std::vector<Complex> phases(1 << n_qubits, 1.);
  state_perm_t perm(permutation_);
  route_recursive(circ, std::nullopt, phases, perm, get_op_ptr(OpType::Ry, 1));
  Eigen::VectorXcd corrections(1 << n_qubits);
  for (unsigned i = 0; i < phases.size(); i++) {
    corrections[i] = 1. / phases[i];
  }
  DiagonalBox diag(corrections);
  circ.add_box(diag, circ.all_qubits());
  circ_ = std::make_shared<Circuit>(circ);
}

nlohmann::json StatePermutationBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const StatePermutationBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["permutation"] = box.get_permutation();
  return j;
}

Op_ptr StatePermutationBox::from_json(const nlohmann::json &j) {
  StatePermutationBox box =
      StatePermutationBox(j.at("permutation").get<state_perm_t>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(StatePermutationBox, StatePermutationBox)

}  // namespace tket