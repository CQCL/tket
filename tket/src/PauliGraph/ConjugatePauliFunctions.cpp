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

#include "ConjugatePauliFunctions.hpp"

#include "PauliGraph/PauliGraph.hpp"
#include "Utils/Assert.hpp"
#include "Utils/Exceptions.hpp"

namespace tket {

std::pair<Pauli, bool> conjugate_Pauli(OpType op, Pauli p, bool reverse) {
  static const std::map<std::pair<OpType, Pauli>, std::pair<Pauli, bool>>
      conj_lut{
          {{OpType::H, Pauli::I}, {Pauli::I, false}},
          {{OpType::H, Pauli::X}, {Pauli::Z, false}},
          {{OpType::H, Pauli::Y}, {Pauli::Y, true}},
          {{OpType::H, Pauli::Z}, {Pauli::X, false}},
          {{OpType::S, Pauli::I}, {Pauli::I, false}},
          {{OpType::S, Pauli::X}, {Pauli::Y, true}},
          {{OpType::S, Pauli::Y}, {Pauli::X, false}},
          {{OpType::S, Pauli::Z}, {Pauli::Z, false}},
          {{OpType::Sdg, Pauli::I}, {Pauli::I, false}},
          {{OpType::Sdg, Pauli::X}, {Pauli::Y, false}},
          {{OpType::Sdg, Pauli::Y}, {Pauli::X, true}},
          {{OpType::Sdg, Pauli::Z}, {Pauli::Z, false}},
          {{OpType::V, Pauli::I}, {Pauli::I, false}},
          {{OpType::V, Pauli::X}, {Pauli::X, false}},
          {{OpType::V, Pauli::Y}, {Pauli::Z, true}},
          {{OpType::V, Pauli::Z}, {Pauli::Y, false}},
          {{OpType::Vdg, Pauli::I}, {Pauli::I, false}},
          {{OpType::Vdg, Pauli::X}, {Pauli::X, false}},
          {{OpType::Vdg, Pauli::Y}, {Pauli::Z, false}},
          {{OpType::Vdg, Pauli::Z}, {Pauli::Y, true}},
          {{OpType::X, Pauli::I}, {Pauli::I, false}},
          {{OpType::X, Pauli::X}, {Pauli::X, false}},
          {{OpType::X, Pauli::Y}, {Pauli::Y, true}},
          {{OpType::X, Pauli::Z}, {Pauli::Z, true}},
          {{OpType::Y, Pauli::I}, {Pauli::I, false}},
          {{OpType::Y, Pauli::X}, {Pauli::X, true}},
          {{OpType::Y, Pauli::Y}, {Pauli::Y, false}},
          {{OpType::Y, Pauli::Z}, {Pauli::Z, true}},
          {{OpType::Z, Pauli::I}, {Pauli::I, false}},
          {{OpType::Z, Pauli::X}, {Pauli::X, true}},
          {{OpType::Z, Pauli::Y}, {Pauli::Y, true}},
          {{OpType::Z, Pauli::Z}, {Pauli::Z, false}}};

  static const std::map<std::pair<OpType, Pauli>, std::pair<Pauli, bool>>
      rev_conj_lut{
          {{OpType::H, Pauli::I}, {Pauli::I, false}},
          {{OpType::H, Pauli::X}, {Pauli::Z, false}},
          {{OpType::H, Pauli::Y}, {Pauli::Y, true}},
          {{OpType::H, Pauli::Z}, {Pauli::X, false}},
          {{OpType::S, Pauli::I}, {Pauli::I, false}},
          {{OpType::S, Pauli::X}, {Pauli::Y, false}},
          {{OpType::S, Pauli::Y}, {Pauli::X, true}},
          {{OpType::S, Pauli::Z}, {Pauli::Z, false}},
          {{OpType::Sdg, Pauli::I}, {Pauli::I, false}},
          {{OpType::Sdg, Pauli::X}, {Pauli::Y, true}},
          {{OpType::Sdg, Pauli::Y}, {Pauli::X, false}},
          {{OpType::Sdg, Pauli::Z}, {Pauli::Z, false}},
          {{OpType::V, Pauli::I}, {Pauli::I, false}},
          {{OpType::V, Pauli::X}, {Pauli::X, false}},
          {{OpType::V, Pauli::Y}, {Pauli::Z, false}},
          {{OpType::V, Pauli::Z}, {Pauli::Y, true}},
          {{OpType::Vdg, Pauli::I}, {Pauli::I, false}},
          {{OpType::Vdg, Pauli::X}, {Pauli::X, false}},
          {{OpType::Vdg, Pauli::Y}, {Pauli::Z, true}},
          {{OpType::Vdg, Pauli::Z}, {Pauli::Y, false}},
          {{OpType::X, Pauli::I}, {Pauli::I, false}},
          {{OpType::X, Pauli::X}, {Pauli::X, false}},
          {{OpType::X, Pauli::Y}, {Pauli::Y, true}},
          {{OpType::X, Pauli::Z}, {Pauli::Z, true}},
          {{OpType::Y, Pauli::I}, {Pauli::I, false}},
          {{OpType::Y, Pauli::X}, {Pauli::X, true}},
          {{OpType::Y, Pauli::Y}, {Pauli::Y, false}},
          {{OpType::Y, Pauli::Z}, {Pauli::Z, true}},
          {{OpType::Z, Pauli::I}, {Pauli::I, false}},
          {{OpType::Z, Pauli::X}, {Pauli::X, true}},
          {{OpType::Z, Pauli::Y}, {Pauli::Y, true}},
          {{OpType::Z, Pauli::Z}, {Pauli::Z, false}}};

  if (reverse) {
    return rev_conj_lut.at({op, p});
  }
  return conj_lut.at({op, p});
}

void conjugate_PauliTensor(
    QubitPauliTensor& qpt, OpType op, const Qubit& q, bool reverse) {
  QubitPauliMap::iterator it = qpt.string.map.find(q);
  if (it == qpt.string.map.end()) {
    return;
  }
  std::pair<Pauli, bool> conj = conjugate_Pauli(op, it->second, reverse);
  it->second = conj.first;
  if (conj.second) {
    qpt.coeff *= -1;
  }
}

void conjugate_PauliTensor(
    QubitPauliTensor& qpt, OpType op, const Qubit& q0, const Qubit& q1) {
  static const std::map<std::pair<Pauli, Pauli>, std::tuple<Pauli, Pauli, bool>>
      cx_conj_lut{
          {{Pauli::I, Pauli::I}, {Pauli::I, Pauli::I, false}},
          {{Pauli::I, Pauli::X}, {Pauli::I, Pauli::X, false}},
          {{Pauli::I, Pauli::Y}, {Pauli::Z, Pauli::Y, false}},
          {{Pauli::I, Pauli::Z}, {Pauli::Z, Pauli::Z, false}},
          {{Pauli::X, Pauli::I}, {Pauli::X, Pauli::X, false}},
          {{Pauli::X, Pauli::X}, {Pauli::X, Pauli::I, false}},
          {{Pauli::X, Pauli::Y}, {Pauli::Y, Pauli::Z, false}},
          {{Pauli::X, Pauli::Z}, {Pauli::Y, Pauli::Y, true}},
          {{Pauli::Y, Pauli::I}, {Pauli::Y, Pauli::X, false}},
          {{Pauli::Y, Pauli::X}, {Pauli::Y, Pauli::I, false}},
          {{Pauli::Y, Pauli::Y}, {Pauli::X, Pauli::Z, true}},
          {{Pauli::Y, Pauli::Z}, {Pauli::X, Pauli::Y, false}},
          {{Pauli::Z, Pauli::I}, {Pauli::Z, Pauli::I, false}},
          {{Pauli::Z, Pauli::X}, {Pauli::Z, Pauli::X, false}},
          {{Pauli::Z, Pauli::Y}, {Pauli::I, Pauli::Y, false}},
          {{Pauli::Z, Pauli::Z}, {Pauli::I, Pauli::Z, false}},
      };
  if (op != OpType::CX) {
    throw NotImplemented("Conjugations of Pauli strings only defined for CXs");
  }
  QubitPauliMap::iterator it0 = qpt.string.map.find(q0);
  QubitPauliMap::iterator it1 = qpt.string.map.find(q1);
  Pauli p0, p1;
  if (it0 == qpt.string.map.end()) {
    p0 = Pauli::I;
  } else {
    p0 = it0->second;
  }
  if (it1 == qpt.string.map.end()) {
    p1 = Pauli::I;
  } else {
    p1 = it1->second;
  }
  std::tuple<Pauli, Pauli, bool> conj = cx_conj_lut.at({p0, p1});
  qpt.string.map[q0] = std::get<0>(conj);
  qpt.string.map[q1] = std::get<1>(conj);
  if (std::get<2>(conj)) {
    qpt.coeff *= -1;
  }
}

void conjugate_PauliTensor(
    QubitPauliTensor& qpt, OpType op, const Qubit& q0, const Qubit& q1,
    const Qubit& q2) {
  /* XXPhase3 gates used for conjugations always implicitly use angle π/2
   * i.e. XXPhase3(1/2). Note that up to phase the 3-qb gate is self-inverse:
   *            XXPhase3(1/2) == exp(π/2i) XXPhase3(-1/2).
   *
   * We conjugate XXPhase3(1/2) by conjugating its CX-equivalent circuit
   *    __________
   * --|          |--               ------- X--H--C--H--X--
   *   | XXPhase3 |                         |     |
   * --|  (1/2)   |-- = exp(-3π/4i) --H--C--C--H--+-----X--
   *   |          |                      |        |
   * --|__________|--               -----X--------X-----X--
   *
   *  We can safely ignore phase differences as for each conjugated XXPhase3,
   *  we insert XXPhase3(1/2) and its dagger XXPhase3(-1/2), so that the phases
   *  always cancel out.
   */
  if (op != OpType::XXPhase3) {
    throw NotImplemented(
        "3qb-Conjugations of Pauli strings only defined for XXPhase3");
  }
  Conjugations equiv = {{OpType::H, {q1}},      {OpType::CX, {q1, q2}},
                        {OpType::CX, {q1, q0}}, {OpType::H, {q0}},
                        {OpType::H, {q1}},      {OpType::CX, {q0, q2}},
                        {OpType::H, {q0}},      {OpType::X, {q0}},
                        {OpType::X, {q1}},      {OpType::X, {q2}}};

  for (auto [op, qbs] : equiv) {
    if (qbs.size() == 1) {
      conjugate_PauliTensor(qpt, op, qbs[0]);
    } else {
      TKET_ASSERT(qbs.size() == 2);
      conjugate_PauliTensor(qpt, op, qbs[0], qbs[1]);
    }
  }
}

}  // namespace tket
