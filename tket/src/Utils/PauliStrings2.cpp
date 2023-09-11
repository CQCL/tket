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

#include "tket/Utils/PauliStrings2.hpp"

namespace tket {

std::set<Qubit> common_qubits(
    const QubitPauliMap &first, const QubitPauliMap &second) {
  std::set<Qubit> common;
  for (const std::pair<const Qubit, Pauli> &p : first) {
    if (p.second == Pauli::I) continue;
    QubitPauliMap::const_iterator found = second.find(p.first);
    if (found != second.end() && found->second == p.second)
      common.insert(p.first);
  }
  return common;
}

std::set<unsigned> common_indices(
    const DensePauliMap &first, const DensePauliMap &second) {
  std::set<unsigned> common;
  unsigned min_size = std::min(first.size(), second.size());
  for (unsigned i = 0; i < min_size; ++i) {
    Pauli p = first.at(i);
    if (p != Pauli::I && p == second.at(i)) common.insert(i);
  }
  return common;
}

std::set<Qubit> own_qubits(
    const QubitPauliMap &first, const QubitPauliMap &second) {
  std::set<Qubit> own;
  for (const std::pair<const Qubit, Pauli> &p : first) {
    if (p.second == Pauli::I) continue;
    QubitPauliMap::const_iterator found = second.find(p.first);
    if (found == second.end() || found->second == Pauli::I) own.insert(p.first);
  }
  return own;
}

std::set<unsigned> own_indices(
    const DensePauliMap &first, const DensePauliMap &second) {
  std::set<unsigned> own;
  unsigned min_size = std::min(first.size(), second.size());
  for (unsigned i = 0; i < min_size; ++i) {
    if (first.at(i) != Pauli::I && second.at(i) == Pauli::I) own.insert(i);
  }
  for (unsigned i = min_size; i < first.size(); ++i) {
    if (first.at(i) != Pauli::I) own.insert(i);
  }
  return own;
}

std::set<Qubit> conflicting_qubits(
    const QubitPauliMap &first, const QubitPauliMap &second) {
  std::set<Qubit> conflicts;
  for (const std::pair<const Qubit, Pauli> &p : first) {
    if (p.second == Pauli::I) continue;
    QubitPauliMap::const_iterator found = second.find(p.first);
    if (found != second.end() && found->second != Pauli::I &&
        found->second != p.second)
      conflicts.insert(p.first);
  }
  return conflicts;
}

std::set<unsigned> conflicting_indices(
    const DensePauliMap &first, const DensePauliMap &second) {
  std::set<unsigned> conflicts;
  unsigned min_size = std::min(first.size(), second.size());
  for (unsigned i = 0; i < min_size; ++i) {
    Pauli p = first.at(i);
    Pauli p2 = second.at(i);
    if (p != Pauli::I && p2 != Pauli::I & p != p2) conflicts.insert(i);
  }
  return conflicts;
}

}  // namespace tket
