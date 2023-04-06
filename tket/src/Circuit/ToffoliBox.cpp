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

#include "Circuit/ToffoliBox.hpp"

#include <boost/graph/max_cardinality_matching.hpp>

#include "Circuit/Circuit.hpp"
#include "Circuit/DiagonalBox.hpp"
#include "Circuit/Multiplexor.hpp"
#include "Gate/Rotation.hpp"
#include "Ops/OpJsonFactory.hpp"
#include "Utils/HelperFunctions.hpp"
#include "Utils/Json.hpp"

namespace tket {

ToffoliBox::ToffoliBox(
    const state_perm_t &permutation, const OpType &rotation_axis)
    : Box(OpType::ToffoliBox),
      permutation_(permutation),
      rotation_axis_(rotation_axis) {
  // we need to check every element has the same size, and every bitstrings
  // appears twice
  if (permutation.size() == 0) {
    throw std::invalid_argument(
        "The permutation argument passed to ToffoliBox is empty.");
  }
  if (rotation_axis != OpType::Rx && rotation_axis != OpType::Ry) {
    throw std::invalid_argument(
        "The rotation_axis argument passed to ToffoliBox must be Rx or Ry.");
  }
  // verify the permutation is valid
  auto it = permutation.begin();
  if (it->first.size() > 32) {
    throw std::invalid_argument(
        "ToffoliBox only supports permutation up to 32 bits.");
  }
  unsigned n_qubits = (unsigned)it->first.size();
  std::set<std::vector<bool>> rhs_states;
  std::set<std::vector<bool>> lhs_states;
  for (; it != permutation.end(); ++it) {
    if (it->first.size() != n_qubits || it->second.size() != n_qubits) {
      throw std::invalid_argument(
          "The permutation argument passed to ToffoliBox contains bitstrings "
          "with different sizes.");
    }
    lhs_states.insert(it->first);
    rhs_states.insert(it->second);
  }
  if (lhs_states != rhs_states) {
    throw std::invalid_argument(
        "The permutation argument passed to ToffoliBox is not complete because "
        "some states aren't mapped.");
  }
}

ToffoliBox::ToffoliBox(const ToffoliBox &other)
    : Box(other),
      permutation_(other.permutation_),
      rotation_axis_(other.rotation_axis_) {}

Op_ptr ToffoliBox::dagger() const {
  state_perm_t reverse_perm;
  for (auto it = permutation_.begin(); it != permutation_.end(); it++) {
    reverse_perm.insert({it->second, it->first});
  }
  return std::make_shared<ToffoliBox>(reverse_perm, rotation_axis_);
}

Op_ptr ToffoliBox::transpose() const { return dagger(); }

op_signature_t ToffoliBox::get_signature() const {
  op_signature_t qubits(
      (unsigned)permutation_.begin()->first.size(), EdgeType::Quantum);
  return qubits;
}

state_perm_t ToffoliBox::get_permutation() const { return permutation_; }

OpType ToffoliBox::get_rotation_axis() const { return rotation_axis_; }

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
    cube_graph_t;
typedef boost::graph_traits<cube_graph_t>::vertex_descriptor perm_vert_t;

static std::vector<std::vector<bool>> rearrange_along_col(
    const std::vector<bool> &prefix, unsigned n_qubits, unsigned col_idx,
    state_perm_t &perm) {
  unsigned n_right_columns = n_qubits - prefix.size() - 1;
  // Given a group of rows identified by a shared prefix and col_idx in {0,1}
  // we rearrange these rows such that the bits(called postfix) after col_idx
  // in their indices span the entire n_right_columns-bit space
  // the rearrangement is done by swapping along col_idx

  // We only need to do the matching for the rows starting with prefix+0
  // the solution will also take care of the rows starting with prefix+1

  // Construct the bipartite graph with 2*2^n_right_columns vertices
  // that connects the postfix in the row entires and the postfix in the row
  // indices. First half of the vertices represents the row entries the second
  // half represents the row indices
  cube_graph_t g(1ULL << (n_right_columns + 1));
  for (unsigned long long postfix_dec = 0;
       postfix_dec < (1ULL << n_right_columns); postfix_dec++) {
    std::vector<bool> postfix = dec_to_bin(postfix_dec, n_right_columns);
    // row0 = prefix+0+postfix
    // row1 = prefix+1+postfix
    std::vector<bool> row0(prefix);
    std::vector<bool> row1(prefix);
    row0.push_back(0);
    row0.insert(row0.end(), postfix.begin(), postfix.end());
    row1.push_back(1);
    row1.insert(row1.end(), postfix.begin(), postfix.end());
    // the postfix of row0's index
    std::vector<bool> idx0_postfix(
        perm[row0].begin() + prefix.size() + 1, perm[row0].end());
    // the postfix of row1's index
    std::vector<bool> idx1_postfix(
        perm[row1].begin() + prefix.size() + 1, perm[row1].end());
    // row0 can stay where it is
    boost::add_edge(
        postfix_dec, bin_to_dec(idx0_postfix) + (1ULL << n_right_columns), g);
    if (idx0_postfix != idx1_postfix) {
      // row0 can also move to row1's index by a swap
      boost::add_edge(
          postfix_dec, bin_to_dec(idx1_postfix) + (1ULL << n_right_columns), g);
    }
  }
  // find a matching
  std::vector<perm_vert_t> match(1ULL << (n_right_columns + 1));
  // O(|V|*|E|*alpha(|V|, |E|)), alpha is a slow growing function
  // that never exceeds 4 for any realistic input.
  // TODO: use a specialised algo for bipartite matching
  // e.g. Hopcroft-Karp for O(|E| sqrt(|V|))
  edmonds_maximum_cardinality_matching(g, &match[0]);
  std::vector<std::vector<bool>> swap_pairs;
  for (unsigned long long postfix_dec = 0;
       postfix_dec < (1ULL << n_right_columns); postfix_dec++) {
    std::vector<bool> postfix = dec_to_bin(postfix_dec, n_right_columns);
    std::vector<bool> row0(prefix);
    row0.push_back(0);
    row0.insert(row0.end(), postfix.begin(), postfix.end());
    std::vector<bool> mapped_idx = dec_to_bin(
        match[postfix_dec] - (1ULL << n_right_columns), n_right_columns);
    if (!std::equal(
            perm[row0].begin() + prefix.size() + 1, perm[row0].end(),
            mapped_idx.begin())) {
      // row0 matched to row1's index
      row0.erase(row0.begin() + col_idx);
      swap_pairs.push_back(row0);
    }
  }
  return swap_pairs;
}

static void swap_rows(
    const std::vector<bool> &pair, unsigned col_idx, state_perm_t &perm,
    ctrl_op_map_t &op_map, std::vector<Complex> &phases,
    const Op_ptr &zflip_op) {
  // swap a pair of rows along col_idx
  // assume two rows only differ at q
  // update op_map, phases, and permutation
  std::vector<bool> row0(pair);
  std::vector<bool> row1(pair);
  row0.insert(row0.begin() + col_idx, 0);
  row1.insert(row1.begin() + col_idx, 1);
  std::swap(perm.at(row0), perm.at(row1));
  // phases is a vector, convert rows to decimal to track phases
  unsigned row0_dec = bin_to_dec(row0);
  unsigned row1_dec = bin_to_dec(row1);
  // we currently only support Rx(pi) and Ry(pi) for permuting states
  // we might introduce Rx(-pi) or Ry(-pi) in the future. Also, dynamically
  // choosing between these 4 might give possible phase cancellation
  // opportunities
  if (zflip_op->get_type() == OpType::Rx &&
      equiv_val(zflip_op->get_params()[0], 1., 4)) {
    // check the Rx angle is (1 mod 4)
    phases[row0_dec] *= -i_;
    phases[row1_dec] *= -i_;
  } else if (
      zflip_op->get_type() == OpType::Ry &&
      equiv_val(zflip_op->get_params()[0], 1., 4)) {
    // check the Ry angle is (1 mod 4)
    phases[row1_dec] *= -1;
  } else {
    // shouldn't happen since we assume the zflip_op must satisfy one of the
    // above conditions. Still throw an error in case something goes wrong.
    throw std::logic_error(
        "Attempt to perform state permutation in ToffoliBox with unsupported "
        "rotations.");
  }
  std::swap(phases[row0_dec], phases[row1_dec]);
  op_map.insert({pair, zflip_op});
}

static std::vector<unsigned> get_multiplexor_args(
    unsigned n_qubits, unsigned col_idx) {
  std::vector<unsigned> multplx_args(n_qubits);
  std::iota(std::begin(multplx_args), std::end(multplx_args), 0);
  multplx_args.erase(multplx_args.begin() + col_idx);
  multplx_args.push_back(col_idx);
  return multplx_args;
}

/**
 * @brief Construct a state permutation circuit
 *
 * @param perm permutation
 * @param n_qubits number of qubits
 * @param zflip_op the base 1-q rotation for swap two states. Currently limited
 * to Rx(1) and Ry(1)
 * @return Circuit
 */
static Circuit permute(
    state_perm_t &perm, unsigned n_qubits, const Op_ptr &zflip_op) {
  // Consider the permutation map as a boolean matrix with n columns and 2^n
  // rows. Row i contains current location of the coefficient that needs to be
  // permutated to the state |i>. We want to sort the rows such that the value
  // of each row matches the its index.
  //
  // The algorithm has two steps:
  // Step 1. We traverse columns to the right, finish at (n-2)th column.
  // At jth column, we first group the rows by their first j+1 entries,
  // so each group is identified by a bitstring [b_0,b_1,...,b_j]. We then
  // rearrange the rows such that for each group, the last n-j-1 bits of
  // its row indices span the entire (n-j-1)-bit space. Such rearrangement
  // can be accomplished with one multiplexor targeting q[j].
  //
  // Step 2. We traverse all columns from right to left.
  // At jth column j, for any row R whose jth bit doesn't match the jth bit of
  // its row index, swap it with row R' where R and R' only differ at the jth
  // bit. All the swaps can be done with one multiplexor targeting q[j].

  // In the implementation, we use a map from rows to rows indices to represent
  // the matrix defined in the algorithm.

  Circuit circ(n_qubits);
  // special case
  if (n_qubits == 1) {
    if (perm.begin()->first != perm.begin()->second) {
      circ.add_op<unsigned>(OpType::X, {0});
    }
    return circ;
  }
  // track relative phases caused by using pi rotations
  std::vector<Complex> phases(1ULL << n_qubits, 1.);
  // step 1
  for (unsigned col_idx = 0; col_idx < n_qubits - 1; col_idx++) {
    ctrl_op_map_t op_map;
    for (unsigned long long prefix = 0; prefix < (1ULL << col_idx); prefix++) {
      std::vector<bool> prefix_bin;
      if (col_idx != 0) {
        // if col_idx == 0, the prefix is empty, otherwise convert the decimal
        // prefix to its binary representation
        prefix_bin = dec_to_bin(prefix, col_idx);
      }
      std::vector<std::vector<bool>> swap_pairs =
          rearrange_along_col(prefix_bin, n_qubits, col_idx, perm);
      for (const std::vector<bool> &pair : swap_pairs) {
        // swap_rows mutates the perm (current permutation), phases (phases
        // accumulated by using SU2 gates) and also updates the op_map to
        // indicate which pairs of rows to swap
        swap_rows(pair, col_idx, perm, op_map, phases, zflip_op);
      }
    }
    if (!op_map.empty()) {
      std::vector<unsigned> multplx_args =
          get_multiplexor_args(n_qubits, col_idx);
      circ.add_box(MultiplexedRotationBox(op_map), multplx_args);
    }
  }
  // step 2
  for (unsigned col_idx = n_qubits; col_idx-- > 0;) {
    ctrl_op_map_t op_map;
    for (unsigned long long row = 0; row < (1ULL << (n_qubits - 1)); row++) {
      std::vector<bool> pair = dec_to_bin(row, n_qubits - 1);
      pair.insert(pair.begin() + col_idx, 0);
      if (perm[pair][col_idx] != 0) {
        pair.erase(pair.begin() + col_idx);
        swap_rows(pair, col_idx, perm, op_map, phases, zflip_op);
      }
    }
    if (!op_map.empty()) {
      std::vector<unsigned> multplx_args =
          get_multiplexor_args(n_qubits, col_idx);
      circ.add_box(MultiplexedRotationBox(op_map), multplx_args);
    }
  }
  // correct the phases with a diagonal operator
  Eigen::VectorXcd corrections(1ULL << n_qubits);
  for (unsigned i = 0; i < phases.size(); i++) {
    corrections[i] = 1. / phases[i];
  }
  DiagonalBox diag(corrections);
  circ.add_box(diag, circ.all_qubits());
  return circ;
}

void ToffoliBox::generate_circuit() const {
  unsigned n_qubits = (unsigned)permutation_.begin()->first.size();
  state_perm_t perm(permutation_);
  // fill the permutation with identities
  if (perm.size() != (1ULL << n_qubits)) {
    for (unsigned long long i = 0; i < (1ULL << n_qubits); i++) {
      std::vector<bool> state = dec_to_bin(i, n_qubits);
      if (perm.find(state) == perm.end()) {
        perm.insert({state, state});
      }
    }
  }
  circ_ = std::make_shared<Circuit>(
      permute(perm, n_qubits, get_op_ptr(rotation_axis_, 1)));
}

nlohmann::json ToffoliBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const ToffoliBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["permutation"] = box.get_permutation();
  j["rotation_axis"] = box.get_rotation_axis();
  return j;
}

Op_ptr ToffoliBox::from_json(const nlohmann::json &j) {
  ToffoliBox box = ToffoliBox(
      j.at("permutation").get<state_perm_t>(),
      j.at("rotation_axis").get<OpType>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(ToffoliBox, ToffoliBox)

}  // namespace tket