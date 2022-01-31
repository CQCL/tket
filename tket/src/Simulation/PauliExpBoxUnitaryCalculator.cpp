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

#include "PauliExpBoxUnitaryCalculator.hpp"

#include <algorithm>

#include "Circuit/Boxes.hpp"
#include "Utils/Assert.hpp"
#include "Utils/Exceptions.hpp"
#include "Utils/PauliStrings.hpp"

namespace tket {
namespace tket_sim {
namespace internal {

// A nonzero entry in the matrix containing only 1,0,-1 values,
// which we will build up for the tensor product
// (taking out a factor of i from Y).
typedef Eigen::Triplet<int> Entry;

static std::map<Pauli, std::array<Entry, 2>> get_pauli_map() {
  std::map<Pauli, std::array<Entry, 2>> result;
  {
    auto& i_matr = result[Pauli::I];
    i_matr[0] = Entry(0, 0, 1);
    i_matr[1] = Entry(1, 1, 1);
  }
  {
    auto& x_matr = result[Pauli::X];
    x_matr[0] = Entry(0, 1, 1);
    x_matr[1] = Entry(1, 0, 1);
  }
  {
    // Remove a factor of i.
    auto& y_matr = result[Pauli::Y];
    y_matr[0] = Entry(0, 1, -1);
    y_matr[1] = Entry(1, 0, 1);
  }
  {
    auto& z_matr = result[Pauli::Z];
    z_matr[0] = Entry(0, 0, 1);
    z_matr[1] = Entry(1, 1, -1);
  }
  return result;
}

// The tensor product of the left and right entries.
static Entry get_combined_entry(const Entry& left, const Entry& right) {
  return Entry(
      2 * left.row() + right.row(), 2 * left.col() + right.col(),
      left.value() * right.value());
}

namespace {
// To save repeated reallocation, build the data up inside here.
struct PauliExpBoxUnitaryCalculator {
  // The two nonzero entries of the 2x2 matrix.
  const std::map<Pauli, std::array<Entry, 2>> pauli_map;

  // The tensor product matrix, all factors of i removed.
  // Although QubitPauliString and QubitPauliTensor could probably
  // be made to do the work for us, they are more complicated
  // than we need.
  std::vector<Entry> sparse_matrix;

  // The number of Y which occurred.
  unsigned power_of_i;

  // A work vector to avoid repeated reallocation.
  std::vector<TripletCd> triplets;

  // We need to remember which diagonal entries were filled,
  // to construct the final triplets.
  std::vector<unsigned> set_diagonals;

  // In the tensor product of the given entry in "sparse_matrix"
  // with the new 2x2 Pauli matrix on the right,
  // add the new entries to "sparse_matrix"
  // (including changing the existing entry).
  // I.e., each entry in "sparse_matrix" is a 1 or -1,
  // which adds two more +/-1 values to the list of nonzero entries
  // (whilst being removed itself),
  // since we take out a factor of i from Y.
  void add_entries(unsigned sparse_matrix_index, Pauli pauli);

  // Call this to begin entering a new Pauli string "manually".
  void clear();

  // Add a single Pauli to the RIGHT of the current tensor product.
  void append(Pauli pauli);

  void fill_triplets(double phase);

  PauliExpBoxUnitaryCalculator();

  static PauliExpBoxUnitaryCalculator& get();
};

PauliExpBoxUnitaryCalculator& PauliExpBoxUnitaryCalculator::get() {
  static PauliExpBoxUnitaryCalculator calculator;
  return calculator;
}

PauliExpBoxUnitaryCalculator::PauliExpBoxUnitaryCalculator()
    : pauli_map(get_pauli_map()) {}

void PauliExpBoxUnitaryCalculator::clear() {
  sparse_matrix.resize(1);
  sparse_matrix[0] = Entry(0, 0, 1);
  power_of_i = 0;
  triplets.clear();
}

void PauliExpBoxUnitaryCalculator::add_entries(
    unsigned sparse_matrix_index, Pauli pauli) {
  TKET_ASSERT(sparse_matrix_index < sparse_matrix.size());
  const auto& single_pauli = pauli_map.at(pauli);
  sparse_matrix.push_back(
      get_combined_entry(sparse_matrix[sparse_matrix_index], single_pauli[0]));
  // Overwrite existing entry.
  // The reference is valid until after the new entry returns!
  auto& existing_entry = sparse_matrix[sparse_matrix_index];
  existing_entry = get_combined_entry(existing_entry, single_pauli[1]);
}

void PauliExpBoxUnitaryCalculator::append(Pauli pauli) {
  const auto current_size = sparse_matrix.size();
  if (pauli == Pauli::Y) {
    ++power_of_i;
  }
  for (unsigned ii = 0; ii < current_size; ++ii) {
    add_entries(ii, pauli);
  }
}

void PauliExpBoxUnitaryCalculator::fill_triplets(double phase) {
  triplets.reserve(sparse_matrix.size());
  triplets.clear();

  // If   M^2 = I   then exp(itM) = cos(t)I + i.sin(t)M.
  const double angle = -0.5 * PI * phase;
  const double cc = std::cos(angle);
  const double ss = std::sin(angle);

  // TODO: efficiency: (i) check if cos(t) or sin(t) is nearly zero;
  // (ii): fill the I coefficients FIRST, so we don't have to bother with
  //      "set_diagonals" (as we can overwrite diagonal entries directly).
  const std::complex<double> identity_coefficient(cc);

  // We took out factors of i from M (coming from Pauli Y),
  // so have to put them back in.
  power_of_i = power_of_i % 4;
  const auto matrix_coefficient =
      std::polar(1.0, 0.5 * (power_of_i + 1) * PI) * ss;

  // When we form the sparse representation of the matrix exp(itM),
  // we must check which diagonal entries were missed out,
  // coming from the cos(t)I term.
  set_diagonals.clear();

  for (const auto& entry : sparse_matrix) {
    // In Utils\PauliStrings.hpp, there is
    //  QubitPauliTensor operator*(Complex a, const QubitPauliTensor &qpt),
    // which messes up ordinary std statements like std::complex * int.
    // Not really a problem, but can be unexpected.
    std::complex<double> value =
        matrix_coefficient * std::complex<double>(entry.value());

    if (entry.row() == entry.col()) {
      set_diagonals.push_back(entry.row());
      value += identity_coefficient;
    }
    triplets.emplace_back(entry.row(), entry.col(), value);
  }
  // Complexity: let there be k diagonals set, out of N, and you want
  // the other N-k (remembering that N = 2^n can be large).
  // Here, we storing the diagonals we DON'T want in a std::vector, sort them,
  // then select the OTHER elements. This takes time
  // O(k.log k) + O(N.log k) time (with pretty small O constants,
  // since sorting a std::vector is fast).
  //
  // If instead we use std::set to store all diagonals,
  // remove each as it is set, and finally iterate through the
  // REMAINING elements, it would be
  // O(N.log N) + O((N-k).log N) + O(k.log k) time.
  // This is probably a bit worse.
  //
  // If we really wanted to speed it up, we'd use bit operations.
  std::sort(set_diagonals.begin(), set_diagonals.end());
  for (unsigned ii = 0; ii < sparse_matrix.size(); ++ii) {
    if (!std::binary_search(set_diagonals.cbegin(), set_diagonals.cend(), ii)) {
      triplets.emplace_back(ii, ii, identity_coefficient);
    }
  }
}
}  // namespace

std::vector<TripletCd> get_triplets(const PauliExpBox& box) {
  const auto phase_optional = eval_expr(box.get_phase());
  if (!phase_optional) {
    throw NotImplemented(
        "PauliExpBoxUnitaryCalculator called "
        "with symbolic phase parameter");
  }
  const auto pauli_string = box.get_paulis();
  auto& calculator = PauliExpBoxUnitaryCalculator::get();
  calculator.clear();
  for (auto pauli : pauli_string) {
    calculator.append(pauli);
  }
  calculator.fill_triplets(phase_optional.value());
  return calculator.triplets;
}

}  // namespace internal
}  // namespace tket_sim
}  // namespace tket
