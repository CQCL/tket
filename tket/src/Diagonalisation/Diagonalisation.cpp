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

#include "tket/Diagonalisation/Diagonalisation.hpp"

#include <tkassert/Assert.hpp>

#include "tket/Ops/Op.hpp"
#include "tket/PauliGraph/ConjugatePauliFunctions.hpp"
#include "tket/Utils/UnitID.hpp"

namespace tket {

Circuit mutual_diagonalise(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType cx_config, DiagonalisationMethod diag_meth) {
  switch (diag_meth) {
    case DiagonalisationMethod::Greedy: {
      return mutual_diagonalise_greedy(gadgets, qubits, cx_config);
    }
    case DiagonalisationMethod::JGM: {
      return mutual_diagonalise_JGM(gadgets, qubits, cx_config);
    }
    case DiagonalisationMethod::vdBT_PE: {
      return mutual_diagonalise_vdBT_PE(gadgets, qubits, cx_config);
    }
    case DiagonalisationMethod::vdBT_CX: {
      return mutual_diagonalise_vdBT_CX(gadgets, qubits, cx_config);
    }
    case DiagonalisationMethod::vdBT_greedy1: {
      return mutual_diagonalise_vdBT_greedy(gadgets, qubits, cx_config, false);
    }
    case DiagonalisationMethod::vdBT_greedy2: {
      return mutual_diagonalise_vdBT_greedy(gadgets, qubits, cx_config, true);
    }
    case DiagonalisationMethod::CSW_CZ: {
      return mutual_diagonalise_CSW_CZ(gadgets, qubits, cx_config);
    }
    case DiagonalisationMethod::CSW_CX: {
      return mutual_diagonalise_CSW_CX(gadgets, qubits, cx_config);
    }
    default: {
      throw std::logic_error("Unrecognised diagonalisation method");
    }
  }
}

void check_easy_diagonalise(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> &qubits,
    Circuit &circ) {
  Conjugations conjugations;
  std::set<Qubit>::iterator qb_iter = qubits.begin();
  for (std::set<Qubit>::iterator next = qb_iter; qb_iter != qubits.end();
       qb_iter = next) {
    ++next;
    Pauli p1 = Pauli::I;
    bool remove_qb = true;
    for (const SpSymPauliTensor &gadget : gadgets) {
      Pauli p2 = gadget.get(*qb_iter);
      if (p2 == Pauli::I) continue;
      if (p1 == Pauli::I) {
        p1 = p2;
      } else if (p1 != p2) {
        remove_qb = false;
        break;
      }
    }
    if (remove_qb) {
      switch (p1) {
        case Pauli::I:
          break;
        case Pauli::X:
          conjugations.push_back({OpType::H, {*qb_iter}});
          circ.add_op<Qubit>(OpType::H, {*qb_iter});
          break;
        case Pauli::Y:
          conjugations.push_back({OpType::Vdg, {*qb_iter}});
          circ.add_op<Qubit>(OpType::V, {*qb_iter});
          break;
        case Pauli::Z:
          break;
        default:
          throw std::logic_error(
              "Unknown Pauli encountered in checking diagonalisation");
      }
      qubits.erase(qb_iter);
    }
  }
  for (SpSymPauliTensor &gadget : gadgets) {
    apply_conjugations(gadget, conjugations);
  }
}

std::optional<std::pair<Pauli, Pauli>> check_pair_compatibility(
    const Qubit &qb1, const Qubit &qb2,
    const std::list<SpSymPauliTensor> &gadgets) {
  if (qb1 == qb2) return std::nullopt;

  /* Do exhaustive search for a Pauli A and Pauli B that
  satisfy Theorem */
  std::list<Pauli> paulis{Pauli::Z, Pauli::X, Pauli::Y};
  for (Pauli pauli1 : paulis) {
    for (Pauli pauli2 : paulis) {
      bool found_pair = true;
      for (const SpSymPauliTensor &gadget : gadgets) {
        Pauli inner_p_1 = gadget.get(qb1);
        Pauli inner_p_2 = gadget.get(qb2);

        if (inner_p_1 == Pauli::I || inner_p_1 == pauli1) {
          if (!(inner_p_2 == Pauli::I || inner_p_2 == pauli2)) {
            found_pair = false;
            break;
          }
        }
        if (inner_p_2 == Pauli::I || inner_p_2 == pauli2) {
          if (!(inner_p_1 == Pauli::I || inner_p_1 == pauli1)) {
            found_pair = false;
            break;
          }
        }
      }
      if (found_pair) {
        return std::make_pair(pauli1, pauli2);
      }
    }
  }
  return std::nullopt;
}

void greedy_diagonalise(
    const std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> &qubits,
    Conjugations &conjugations, Circuit &circ, CXConfigType cx_config) {
  unsigned total_counter = UINT_MAX;
  QubitPauliMap to_diag;
  for (const SpSymPauliTensor &gadget : gadgets) {
    unsigned support_counter = 0;
    QubitPauliMap to_diag_candidates;
    for (const Qubit &qb : qubits) {
      Pauli p = gadget.get(qb);
      if (p != Pauli::I) {
        ++support_counter;
        to_diag_candidates.insert({qb, p});
      }
    }
    if (support_counter < total_counter && support_counter > 1) {
      total_counter = support_counter;
      to_diag = to_diag_candidates;
    }
  }
  if (to_diag.empty()) {
    throw std::logic_error("Brute Force Diagonalise can't find a candidate!");
  }
  for (const std::pair<const Qubit, Pauli> &qp : to_diag) {
    const Qubit &qb = qp.first;
    switch (qp.second) {
      case Pauli::X: {
        conjugations.push_back({OpType::H, {qb}});
        circ.add_op<Qubit>(OpType::H, {qb});
        break;
      }
      case Pauli::Y: {
        conjugations.push_back({OpType::Vdg, {qb}});
        circ.add_op<Qubit>(OpType::V, {qb});
        break;
      }
      case Pauli::Z: {
        break;
      }
      case Pauli::I:
      default:
        throw std::logic_error("Unknown Pauli in greedy diagonalisation.");
    }
  }
  qubit_vector_t diag_qubits;
  for (const auto &[k, v] : to_diag) diag_qubits.push_back(k);
  unsigned n_qubits = diag_qubits.size();
  Qubit first_qb = diag_qubits[0];

  switch (cx_config) {
    case CXConfigType::Snake: {
      for (unsigned i = n_qubits - 1; i > 0; --i) {
        Qubit qb = diag_qubits[i], before = diag_qubits[i - 1];
        conjugations.push_back({OpType::CX, {qb, before}});
        circ.add_op<Qubit>(OpType::CX, {qb, before});
      }
      break;
    }
    case CXConfigType::Star: {
      for (unsigned i = 1; i < n_qubits; ++i) {
        Qubit qb = diag_qubits[i];
        conjugations.push_back({OpType::CX, {qb, first_qb}});
        circ.add_op<Qubit>(OpType::CX, {qb, first_qb});
      }
      break;
    }
    case CXConfigType::Tree: {
      unsigned complete_layers = floor(log2(n_qubits));
      unsigned dense_end = pow(2, complete_layers);
      for (unsigned j = 0; j < n_qubits - dense_end; j++) {
        circ.add_op<Qubit>(
            OpType::CX,
            {diag_qubits[dense_end + j], diag_qubits[dense_end - 1 - j]});
        conjugations.push_back(
            {OpType::CX,
             {diag_qubits[dense_end + j], diag_qubits[dense_end - 1 - j]}});
      }
      for (unsigned step_size = 1; step_size < dense_end; step_size *= 2) {
        for (unsigned j = 0; j < dense_end; j += 2 * step_size) {
          circ.add_op<Qubit>(
              OpType::CX, {diag_qubits[j + step_size], diag_qubits[j]});
          conjugations.push_back(
              {OpType::CX, {diag_qubits[j + step_size], diag_qubits[j]}});
        }
      }
      break;
    }
    case CXConfigType::MultiQGate: {
      int sign_correction = 1;
      for (int q = n_qubits - 1; q > 0; q -= 2) {
        Qubit qb = diag_qubits[q];
        if (q - 1 > 0) {
          Qubit before = diag_qubits[q - 1];
          circ.add_op<Qubit>(OpType::H, {qb});
          circ.add_op<Qubit>(OpType::H, {before});
          circ.add_op<Qubit>(OpType::XXPhase3, 0.5, {qb, before, first_qb});
          conjugations.push_back({OpType::H, {qb}});
          conjugations.push_back({OpType::H, {before}});
          conjugations.push_back({OpType::XXPhase3, {qb, before, first_qb}});
          sign_correction *= -1;
        } else {
          circ.add_op<Qubit>(OpType::CX, {qb, first_qb});
          conjugations.push_back({OpType::CX, {qb, first_qb}});
        }
      }
      if (sign_correction < 0) {
        circ.add_op<Qubit>(OpType::X, {first_qb});
        conjugations.push_back({OpType::X, {first_qb}});
      }
      break;
    }
    default:
      throw std::logic_error("Unknown CXConfigType in greedy diagonalisation.");
  }
  qubits.erase(first_qb);
}

/* Diagonalise a set of Pauli Gadgets simultaneously using Cliffords*/
Circuit mutual_diagonalise_greedy(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType cx_config) {
  Circuit cliff_circ;
  for (const Qubit &qb : qubits) {
    cliff_circ.add_qubit(qb);
  }
  check_easy_diagonalise(gadgets, qubits, cliff_circ);
  while (!qubits.empty()) {
    Conjugations conjugations;
    Qubit qb_a;
    Qubit qb_b;
    std::pair<Pauli, Pauli> pauli_pair;
    bool found_match = false;
    /* First, try to find some qubits we can reduce using only 1 CX */
    for (const Qubit &qb1 : qubits) {
      for (const Qubit &qb2 : qubits) {
        std::optional<std::pair<Pauli, Pauli>> compatible =
            check_pair_compatibility(qb1, qb2, gadgets);
        if (compatible.has_value()) {
          pauli_pair = *compatible;
          qb_a = qb1;
          qb_b = qb2;
          found_match = true;
          break;
        }
      }
      if (found_match) break;
    }
    if (found_match) {
      Pauli p1 = pauli_pair.first;
      Pauli p2 = pauli_pair.second;
      switch (p1) {
        case Pauli::X: {
          conjugations.push_back({OpType::H, {qb_a}});
          cliff_circ.add_op<Qubit>(OpType::H, {qb_a});
          break;
        }
        case Pauli::Y: {
          conjugations.push_back({OpType::Vdg, {qb_a}});
          cliff_circ.add_op<Qubit>(OpType::V, {qb_a});
          break;
        }
        case Pauli::Z: {
          break;
        }
        case Pauli::I:
        default:
          throw std::logic_error("Unknown Pauli in mutual diagonalisation.");
      }
      switch (p2) {
        case Pauli::X: {
          conjugations.push_back({OpType::H, {qb_b}});
          cliff_circ.add_op<Qubit>(OpType::H, {qb_b});
          break;
        }
        case Pauli::Y: {
          conjugations.push_back({OpType::Vdg, {qb_b}});
          cliff_circ.add_op<Qubit>(OpType::V, {qb_b});
          break;
        }
        case Pauli::Z: {
          break;
        }
        case Pauli::I:
        default:
          throw std::logic_error("Unknown Pauli in mutual diagonalisation.");
      }
      conjugations.push_back({OpType::CX, {qb_a, qb_b}});
      cliff_circ.add_op<Qubit>(OpType::CX, {qb_a, qb_b});
      qubits.erase(qb_b);
    }
    /* If we can't, do it with `n-1` CXs, where `n` := no. of undiagonalised
     * qubits */
    if (!found_match) {
      greedy_diagonalise(gadgets, qubits, conjugations, cliff_circ, cx_config);
    }
    for (SpSymPauliTensor &gadget : gadgets) {
      apply_conjugations(gadget, conjugations);
    }
    // we may have made some easy-to-remove qubits
    check_easy_diagonalise(gadgets, qubits, cliff_circ);
  }
  return cliff_circ;
}

void apply_conjugations(
    SpSymPauliTensor &qps, const Conjugations &conjugations) {
  SpPauliStabiliser stab(qps.string);
  for (const auto &optype_qubit_pair : conjugations) {
    OpType ot = optype_qubit_pair.first;
    const qubit_vector_t &qbs = optype_qubit_pair.second;
    if (!optypeinfo().at(ot).signature ||
        optypeinfo().at(ot).signature->size() != qbs.size())
      throw std::logic_error("Incompatible qubit count for conjugations");
    switch (ot) {
      case OpType::H:
      case OpType::S:
      case OpType::Sdg:
      case OpType::V:
      case OpType::Vdg:
      case OpType::X:
      case OpType::Z:
        conjugate_PauliTensor(stab, ot, qbs[0]);
        break;
      case OpType::CX:
        conjugate_PauliTensor(stab, ot, qbs[0], qbs[1]);
        break;
      case OpType::XXPhase3:
        conjugate_PauliTensor(stab, ot, qbs[0], qbs[1], qbs[2]);
        break;
      default:
        throw std::logic_error(
            "Unknown OpType received when applying conjugations.");
    }
  }
  qps.string = stab.string;
  qps.coeff *= cast_coeff<quarter_turns_t, Expr>(stab.coeff);
}

}  // namespace tket
