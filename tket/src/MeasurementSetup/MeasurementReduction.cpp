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

#include "tket/MeasurementSetup/MeasurementReduction.hpp"

namespace tket {

MeasurementSetup measurement_reduction(
    const std::list<SpPauliString>& strings, PauliPartitionStrat strat,
    GraphColourMethod method, CXConfigType cx_config) {
  std::set<Qubit> qubits;
  for (const SpPauliString& qpt : strings) {
    for (const std::pair<const Qubit, Pauli>& qb_p : qpt.string)
      qubits.insert(qb_p.first);
  }

  std::map<Qubit, unsigned> qb_location_map;
  unsigned u = 0;
  for (const Qubit& qb : qubits) {
    qb_location_map[qb] = u;
    ++u;
  }

  std::list<std::list<SpPauliString>> all_terms =
      term_sequence(strings, strat, method);
  MeasurementSetup ms;
  unsigned i = 0;
  for (const std::list<SpPauliString>& terms : all_terms) {
    std::list<SpSymPauliTensor> gadgets;
    for (const SpPauliString& string : terms) {
      gadgets.push_back(string);
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

    std::list<SpSymPauliTensor>::const_iterator gadgets_iter = gadgets.begin();
    for (const SpPauliString& string : terms) {
      SpPauliStabiliser stab(*gadgets_iter);  // Force coeff to be real
      std::vector<unsigned> bits;
      for (const std::pair<const Qubit, Pauli>& qp_pair : stab.string) {
        if (qp_pair.second == Pauli::Z)
          bits.push_back(qb_location_map.at(qp_pair.first));
      }
      ms.add_result_for_term(string, {i, bits, stab.is_real_negative()});
      ++gadgets_iter;
    }
    ++i;
  }

  return ms;
}

}  // namespace tket
