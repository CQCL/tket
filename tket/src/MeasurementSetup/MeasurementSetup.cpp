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

#include "MeasurementSetup.hpp"

#include "Converters/Converters.hpp"

namespace tket {

std::string MeasurementSetup::MeasurementBitMap::to_str() const {
  std::stringstream ss;
  ss << "Circuit index: ";
  ss << circ_index;
  ss << "\nBits: ";
  for (unsigned i : bits) {
    ss << i;
    ss << " ";
  }
  ss << "\nInvert: ";
  if (invert)
    ss << "True";
  else
    ss << "False";
  return ss.str();
}

void MeasurementSetup::add_measurement_circuit(const Circuit &circ) {
  measurement_circs.push_back(circ);
}

void MeasurementSetup::add_result_for_term(
    const QubitPauliString &term, const MeasurementBitMap &result) {
  result_map[term].push_back(result);
}

void MeasurementSetup::add_result_for_term(
    const QubitPauliTensor &term, const MeasurementBitMap &result) {
  add_result_for_term(term.string, result);
}

bool MeasurementSetup::verify() const {
  std::map<std::pair<unsigned, unsigned>, QubitPauliTensor> pauli_map;
  // Identify Paulis measured onto each bit
  for (unsigned circ_id = 0; circ_id < measurement_circs.size(); ++circ_id) {
    Circuit circ = measurement_circs[circ_id];
    std::map<Qubit, unsigned> readout = circ.qubit_readout();
    // Remove Measures from circuit to allow CliffordTableau generation
    for (const Vertex &out : circ.q_outputs()) {
      Vertex pred = circ.get_predecessors(out).front();
      if (circ.get_OpType_from_Vertex(pred) == OpType::Measure) {
        circ.remove_vertex(
            pred, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);
      }
    }
    CliffTableau tab = circuit_to_tableau(circ);
    for (const Qubit &qb : tab.get_qubits()) {
      pauli_map.insert({{circ_id, readout[qb]}, tab.get_zpauli(qb)});
    }
  }
  for (const std::pair<const QubitPauliString, std::vector<MeasurementBitMap>>
           &term : result_map) {
    for (const MeasurementBitMap &bits : term.second) {
      QubitPauliTensor total;
      for (unsigned bit : bits.bits) {
        total = total * pauli_map[{bits.circ_index, bit}];
      }
      if (bits.invert) total.coeff *= -1.;
      QubitPauliTensor term_tensor(term.first);
      if (total != term_tensor) {
        std::stringstream out;
        out << "Invalid MeasurementSetup: expecting to measure "
            << term_tensor.to_str() << "; actually measured " << total.to_str();
        tket_log()->error(out.str());
        return false;
      }
    }
  }
  return true;
}

std::string MeasurementSetup::to_str() const {
  std::stringstream ss;
  ss << "Circuits: ";
  ss << measurement_circs.size();
  ss << "\n";
  for (const std::pair<const QubitPauliString, std::vector<MeasurementBitMap>>
           &tensor_map : result_map) {
    ss << "|| ";
    ss << tensor_map.first.to_str();
    ss << " ||\n";
    for (const MeasurementBitMap &mbm : tensor_map.second) {
      ss << mbm.to_str();
      ss << "\n";
    }
  }
  return ss.str();
}

void to_json(
    nlohmann::json &j, const MeasurementSetup::MeasurementBitMap &result) {
  j["circ_index"] = result.get_circ_index();
  j["bits"] = result.get_bits();
  j["invert"] = result.get_invert();
}

void from_json(
    const nlohmann::json &j, MeasurementSetup::MeasurementBitMap &result) {
  result = MeasurementSetup::MeasurementBitMap(
      j.at("circ_index").get<unsigned>(),
      j.at("bits").get<std::vector<unsigned>>(), j.at("invert").get<bool>());
}

void to_json(nlohmann::json &j, const MeasurementSetup &setup) {
  std::vector<std::pair<
      QubitPauliString, std::vector<MeasurementSetup::MeasurementBitMap>>>
      map_list;
  for (const std::pair<
           const QubitPauliString,
           std::vector<MeasurementSetup::MeasurementBitMap>> &tensor_map :
       setup.get_result_map()) {
    map_list.push_back(tensor_map);
  }
  // sort the list for consistent serialisation
  std::sort(map_list.begin(), map_list.end(), [](auto pair1, auto pair2) {
    return pair1.first < pair2.first;
  });
  j["result_map"] = map_list;
  j["circs"] = setup.get_circs();
}

void from_json(const nlohmann::json &j, MeasurementSetup &setup) {
  for (auto it = j["circs"].begin(); it != j["circs"].end(); ++it) {
    setup.add_measurement_circuit(it->get<Circuit>());
  }
  for (auto it = j["result_map"].begin(); it != j["result_map"].end(); ++it) {
    for (auto second_it = it->at(1).begin(); second_it != it->at(1).end();
         ++second_it) {
      setup.add_result_for_term(
          it->at(0).get<QubitPauliString>(),
          second_it->get<MeasurementSetup::MeasurementBitMap>());
    }
  }
}

}  // namespace tket
