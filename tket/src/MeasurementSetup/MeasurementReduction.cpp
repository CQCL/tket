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

#include "MeasurementReduction.hpp"

namespace tket {

MeasurementSetup measurement_reduction(
    const std::list<QubitPauliString>& strings, PauliPartitionStrat strat,
    GraphColourMethod method, CXConfigType cx_config) {
  std::set<Qubit> qubits;
  for (const QubitPauliString& qpt : strings) {
    for (const std::pair<const Qubit, Pauli>& qb_p : qpt.map)
      qubits.insert(qb_p.first);
  }

  std::map<Qubit, unsigned> qb_location_map;
  unsigned u = 0;
  for (const Qubit& qb : qubits) {
    qb_location_map[qb] = u;
    ++u;
  }

  std::list<std::list<QubitPauliString>> all_terms =
      term_sequence(strings, strat, method);
  MeasurementSetup ms;
  unsigned i = 0;
  for (const std::list<QubitPauliString>& terms : all_terms) {
    std::list<std::pair<QubitPauliTensor, Expr>> gadgets;
    for (const QubitPauliString& string : terms) {
      QubitPauliTensor qps(string);
      gadgets.push_back({qps, 1.});
    }

    std::set<Qubit> mutable_qb_set(qubits);
    Circuit cliff_circ = mutual_diagonalise(gadgets, mutable_qb_set, cx_config);
    unsigned bit_count = 0;
    for (const Qubit& qb : cliff_circ.all_qubits()) {
      cliff_circ.add_bit(Bit(bit_count));
      cliff_circ.add_measure(qb, Bit(bit_count));
      ++bit_count;
    }
    ms.add_measurement_circuit(cliff_circ);

    std::list<std::pair<QubitPauliTensor, Expr>>::const_iterator gadgets_iter =
        gadgets.begin();
    for (const QubitPauliString& string : terms) {
      const std::pair<QubitPauliTensor, Expr>& new_gadget = *gadgets_iter;
      std::vector<unsigned> bits;
      for (const std::pair<const Qubit, Pauli>& qp_pair :
           new_gadget.first.string.map) {
        if (qp_pair.second == Pauli::Z)
          bits.push_back(qb_location_map.at(qp_pair.first));
      }
      bool invert = (std::abs(new_gadget.first.coeff + Complex(1)) < EPS);
      ms.add_result_for_term(string, {i, bits, invert});
      ++gadgets_iter;
    }
    ++i;
  }

  return ms;
}

}  // namespace tket
