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

#include "tket/PauliGraph/ConjugatePauliFunctions.hpp"
#include "tket/Utils/UnitID.hpp"

namespace tket {

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
Circuit mutual_diagonalise(
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

std::pair<Circuit, Qubit> reduce_pauli_to_z(
    const SpPauliStabiliser &pauli, CXConfigType cx_config) {
  Circuit circ;
  qubit_vector_t qubits;
  for (const std::pair<const Qubit, Pauli> &qp : pauli.string) {
    circ.add_qubit(qp.first);
    if (qp.second != Pauli::I) qubits.push_back(qp.first);
    switch (qp.second) {
      case Pauli::X: {
        circ.add_op<Qubit>(OpType::H, {qp.first});
        break;
      }
      case Pauli::Y: {
        circ.add_op<Qubit>(OpType::V, {qp.first});
        break;
      }
      default: {
        break;
      }
    }
  }
  unsigned n_qubits = qubits.size();
  if (n_qubits == 0) throw std::logic_error("Cannot reduce identity to Z");
  switch (cx_config) {
    case CXConfigType::Snake: {
      for (unsigned i = n_qubits - 1; i != 0; --i) {
        circ.add_op<Qubit>(OpType::CX, {qubits.at(i), qubits.at(i - 1)});
      }
      break;
    }
    case CXConfigType::Star: {
      for (unsigned i = n_qubits - 1; i != 0; --i) {
        circ.add_op<Qubit>(OpType::CX, {qubits.at(i), qubits.front()});
      }
      break;
    }
    case CXConfigType::Tree: {
      for (unsigned step_size = 1; step_size < n_qubits; step_size *= 2) {
        for (unsigned i = 0; step_size + i < n_qubits; i += 2 * step_size) {
          circ.add_op<Qubit>(
              OpType::CX, {qubits.at(step_size + i), qubits.at(i)});
        }
      }
      break;
    }
    case CXConfigType::MultiQGate: {
      bool flip_phase = false;
      for (unsigned i = n_qubits - 1; i != 0; --i) {
        if (i == 1) {
          circ.add_op<Qubit>(OpType::CX, {qubits.at(i), qubits.front()});
        } else {
          /**
           * This is only equal to the CX decompositions above up to phase,
           * but phase differences are cancelled out by its dagger
           */
          circ.add_op<Qubit>(OpType::H, {qubits.at(i)});
          circ.add_op<Qubit>(OpType::H, {qubits.at(i - 1)});
          circ.add_op<Qubit>(
              OpType::XXPhase3, 0.5,
              {qubits.at(i), qubits.at(i - 1), qubits.front()});
          --i;
          flip_phase = !flip_phase;
        }
      }
      if (flip_phase) circ.add_op<Qubit>(OpType::X, {qubits.front()});
      break;
    }
  }
  return {circ, qubits.front()};
}

static void reduce_shared_qs_by_CX_snake(
    Circuit &circ, std::set<Qubit> &match, SpPauliStabiliser &pauli0,
    SpPauliStabiliser &pauli1) {
  unsigned match_size = match.size();
  while (match_size > 1) {  // We allow one match left over
    auto it = --match.end();
    Qubit to_eliminate = *it;
    match.erase(it);
    Qubit helper = *match.rbegin();
    // extend CX snake
    circ.add_op<Qubit>(OpType::CX, {to_eliminate, helper});
    pauli0.string.erase(to_eliminate);
    pauli1.string.erase(to_eliminate);
    match_size--;
  }
}

static void reduce_shared_qs_by_CX_star(
    Circuit &circ, std::set<Qubit> &match, SpPauliStabiliser &pauli0,
    SpPauliStabiliser &pauli1) {
  std::set<Qubit>::iterator iter = match.begin();
  for (std::set<Qubit>::iterator next = match.begin(); match.size() > 1;
       iter = next) {
    ++next;
    Qubit to_eliminate = *iter;
    circ.add_op<Qubit>(OpType::CX, {to_eliminate, *match.rbegin()});
    pauli0.string.erase(to_eliminate);
    pauli1.string.erase(to_eliminate);
    match.erase(iter);
  }
}

static void reduce_shared_qs_by_CX_tree(
    Circuit &circ, std::set<Qubit> &match, SpPauliStabiliser &pauli0,
    SpPauliStabiliser &pauli1) {
  while (match.size() > 1) {
    std::set<Qubit> remaining;
    std::set<Qubit>::iterator it = match.begin();
    while (it != match.end()) {
      Qubit maintained = *it;
      it++;
      remaining.insert(maintained);
      if (it != match.end()) {
        Qubit to_eliminate = *it;
        it++;
        circ.add_op<Qubit>(OpType::CX, {to_eliminate, maintained});
        pauli0.string.erase(to_eliminate);
        pauli1.string.erase(to_eliminate);
      }
    }
    match = remaining;
  }
}

static void reduce_shared_qs_by_CX_multiqgate(
    Circuit &circ, std::set<Qubit> &match, SpPauliStabiliser &pauli0,
    SpPauliStabiliser &pauli1) {
  if (match.size() <= 1) {
    return;
  }
  // last qubit is target
  Qubit target = *match.rbegin();
  while (match.size() > 1) {
    std::set<Qubit>::iterator iter = match.begin();
    if (match.size() == 2) {
      // use CX
      Qubit to_eliminate = *iter;
      match.erase(iter);
      pauli0.string.erase(to_eliminate);
      pauli1.string.erase(to_eliminate);

      circ.add_op<Qubit>(OpType::CX, {to_eliminate, target});
    } else {
      // use XXPhase3
      Qubit to_eliminate1 = *iter;
      match.erase(iter++);
      pauli0.string.erase(to_eliminate1);
      pauli1.string.erase(to_eliminate1);

      Qubit to_eliminate2 = *iter;
      match.erase(iter);
      pauli0.string.erase(to_eliminate2);
      pauli1.string.erase(to_eliminate2);

      circ.add_op<Qubit>(OpType::H, {to_eliminate1});
      circ.add_op<Qubit>(OpType::H, {to_eliminate2});
      circ.add_op<Qubit>(
          OpType::XXPhase3, 0.5, {to_eliminate1, to_eliminate2, target});
      circ.add_op<Qubit>(OpType::X, {target});
    }
  }
}

std::pair<Circuit, std::optional<Qubit>> reduce_overlap_of_paulis(
    SpPauliStabiliser &pauli0, SpPauliStabiliser &pauli1,
    CXConfigType cx_config, bool allow_matching_final) {
  /*
   * Cowtan, Dilkes, Duncan, Simmons, Sivarajah: Phase Gadget Synthesis for
   * Shallow Circuits, Lemma 4.9
   * Let s and t be Pauli strings; then there exists a Clifford unitary U such
   * that
   * P(a, s) . P(b, t) = U . P(a, s') . P(b, t') . U^\dagger
   * where s' and t' are Pauli strings with intersection at most 1.
   *
   * Follows the procedure to reduce the intersection of the gadgets and then
   * synthesises the remainder individually.
   */

  /*
   * Step 1: Identify qubits in both pauli0 and pauli1 which either match or
   * don't
   */
  std::set<Qubit> match = pauli0.common_qubits(pauli1);
  std::set<Qubit> mismatch = pauli0.conflicting_qubits(pauli1);

  /*
   * Step 2: Build the unitary U that minimises the intersection of the gadgets.
   */
  Circuit u;
  for (const std::pair<const Qubit, Pauli> &qp : pauli0.string)
    u.add_qubit(qp.first);
  for (const std::pair<const Qubit, Pauli> &qp : pauli1.string) {
    if (!u.contains_unit(qp.first)) u.add_qubit(qp.first);
  }

  /*
   * Step 2.i: Remove (almost) all matches by converting to Z basis and applying
   * CXs
   */
  for (const Qubit &qb : match) {
    switch (pauli0.get(qb)) {
      case Pauli::X:
        u.add_op<Qubit>(OpType::H, {qb});
        pauli0.set(qb, Pauli::Z);
        pauli1.set(qb, Pauli::Z);
        break;
      case Pauli::Y:
        u.add_op<Qubit>(OpType::V, {qb});
        pauli0.set(qb, Pauli::Z);
        pauli1.set(qb, Pauli::Z);
        break;
      default:
        break;
    }
  }
  switch (cx_config) {
    case CXConfigType::Snake: {
      reduce_shared_qs_by_CX_snake(u, match, pauli0, pauli1);
      break;
    }
    case CXConfigType::Star: {
      reduce_shared_qs_by_CX_star(u, match, pauli0, pauli1);
      break;
    }
    case CXConfigType::Tree: {
      reduce_shared_qs_by_CX_tree(u, match, pauli0, pauli1);
      break;
    }
    case CXConfigType::MultiQGate: {
      reduce_shared_qs_by_CX_multiqgate(u, match, pauli0, pauli1);
      break;
    }
    default:
      throw std::logic_error(
          "Unknown CXConfigType received when decomposing gadget.");
  }
  /*
   * Step 2.ii: Convert mismatches to Z in pauli0 and X in pauli1
   */
  for (const Qubit &qb : mismatch) {
    switch (pauli0.get(qb)) {
      case Pauli::X: {
        switch (pauli1.get(qb)) {
          case Pauli::Y:
            u.add_op<Qubit>(OpType::Sdg, {qb});
            u.add_op<Qubit>(OpType::Vdg, {qb});
            break;
          case Pauli::Z:
            u.add_op<Qubit>(OpType::H, {qb});
            break;
          default:
            TKET_ASSERT(false);
            break;  // Cannot hit this case
        }
        break;
      }
      case Pauli::Y: {
        switch (pauli1.get(qb)) {
          case Pauli::X:
            u.add_op<Qubit>(OpType::V, {qb});
            break;
          case Pauli::Z:
            u.add_op<Qubit>(OpType::V, {qb});
            u.add_op<Qubit>(OpType::S, {qb});
            break;
          default:
            TKET_ASSERT(false);
            break;  // Cannot hit this case
        }
        break;
      }
      default: {  // Necessarily Z
        if (pauli1.get(qb) == Pauli::Y) u.add_op<Qubit>(OpType::Sdg, {qb});
        // No need to act if already X
      }
    }
    pauli0.set(qb, Pauli::Z);
    pauli1.set(qb, Pauli::X);
  }

  /*
   * Step 2.iii: Remove the final matching qubit against a mismatch if one
   * exists, otherwise remove into another qubit or allow final overlap to be
   * matching
   */
  std::optional<Qubit> last_overlap = std::nullopt;
  if (!match.empty()) {
    Qubit last_match = *match.begin();
    if (!mismatch.empty()) {
      Qubit mismatch_used =
          *mismatch.rbegin();  // Prefer to use the one that may be left over
                               // after reducing pairs
      u.add_op<Qubit>(OpType::S, {mismatch_used});
      u.add_op<Qubit>(OpType::CX, {last_match, mismatch_used});
      u.add_op<Qubit>(OpType::Sdg, {mismatch_used});
      pauli0.string.erase(last_match);
      pauli1.string.erase(last_match);
    } else if (!allow_matching_final) {
      std::optional<std::pair<Qubit, Pauli>> other;
      for (const std::pair<const Qubit, Pauli> &qp : pauli0.string) {
        if (qp.first != last_match && qp.second != Pauli::I) {
          other = qp;
          pauli0.string.erase(last_match);
          break;
        }
      }
      if (!other) {
        for (const std::pair<const Qubit, Pauli> &qp : pauli1.string) {
          if (qp.first != last_match && qp.second != Pauli::I) {
            other = qp;
            pauli1.string.erase(last_match);
            break;
          }
        }
        if (!other)
          throw std::logic_error(
              "Cannot reduce identical Paulis to different qubits");
      }
      if (other->second == Pauli::X) {
        u.add_op<Qubit>(OpType::H, {other->first});
        u.add_op<Qubit>(OpType::CX, {last_match, other->first});
        u.add_op<Qubit>(OpType::H, {other->first});
      } else {
        u.add_op<Qubit>(OpType::CX, {last_match, other->first});
      }
    } else {
      last_overlap = last_match;
    }
  }

  /*
   * Step 2.iv: Reduce pairs of mismatches to different qubits.
   * Allow both gadgets to build a remaining qubit if it exists.
   */
  std::set<Qubit>::iterator mis_it = mismatch.begin();
  while (mis_it != mismatch.end()) {
    Qubit z_in_0 = *mis_it;
    mis_it++;
    if (mis_it != mismatch.end()) {
      Qubit x_in_1 = *mis_it;
      u.add_op<Qubit>(OpType::CX, {x_in_1, z_in_0});
      pauli0.string.erase(x_in_1);
      pauli1.string.erase(z_in_0);
      mis_it++;
    } else {
      last_overlap = z_in_0;
    }
  }

  return {u, last_overlap};
}

std::pair<Circuit, Qubit> reduce_anticommuting_paulis_to_z_x(
    SpPauliStabiliser pauli0, SpPauliStabiliser pauli1,
    CXConfigType cx_config) {
  std::pair<Circuit, std::optional<Qubit>> reduced_overlap =
      reduce_overlap_of_paulis(pauli0, pauli1, cx_config);
  Circuit &u = reduced_overlap.first;
  if (!reduced_overlap.second)
    throw std::logic_error("No overlap for anti-commuting paulis");
  Qubit &last_mismatch = *reduced_overlap.second;

  /**
   * Reduce each remaining Pauli to the shared mismatching qubit.
   * Since reduce_pauli_to_Z does not allow us to pick the final qubit, we
   * reserve the mismatching qubit, call reduce_pauli_to_Z on the rest, and add
   * a CX.
   */
  pauli0.string.erase(last_mismatch);
  pauli0.compress();
  if (!pauli0.string.empty()) {
    std::pair<Circuit, Qubit> diag0 = reduce_pauli_to_z(pauli0, cx_config);
    u.append(diag0.first);
    u.add_op<Qubit>(OpType::CX, {diag0.second, last_mismatch});
  }
  pauli1.compress();
  pauli1.string.erase(last_mismatch);
  if (!pauli1.string.empty()) {
    std::pair<Circuit, Qubit> diag1 = reduce_pauli_to_z(pauli1, cx_config);
    u.append(diag1.first);
    u.add_op<Qubit>(OpType::H, {last_mismatch});
    u.add_op<Qubit>(OpType::CX, {diag1.second, last_mismatch});
    u.add_op<Qubit>(OpType::H, {last_mismatch});
  }

  return {u, last_mismatch};
}

std::tuple<Circuit, Qubit, Qubit> reduce_commuting_paulis_to_zi_iz(
    SpPauliStabiliser pauli0, SpPauliStabiliser pauli1,
    CXConfigType cx_config) {
  std::pair<Circuit, std::optional<Qubit>> reduced_overlap =
      reduce_overlap_of_paulis(pauli0, pauli1, cx_config);
  Circuit &u = reduced_overlap.first;
  if (reduced_overlap.second)
    throw std::logic_error("Overlap remaining for commuting paulis");

  /**
   * Reduce each remaining Pauli to a single qubit.
   */
  std::pair<Circuit, Qubit> diag0 = reduce_pauli_to_z(pauli0, cx_config);
  u.append(diag0.first);
  std::pair<Circuit, Qubit> diag1 = reduce_pauli_to_z(pauli1, cx_config);
  u.append(diag1.first);

  return {u, diag0.second, diag1.second};
}

}  // namespace tket
