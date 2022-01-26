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

#pragma once

#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Architecture/Architecture.hpp"
#include "Characterisation/ErrorTypes.hpp"
#include "Circuit/Circuit.hpp"
#include "Utils/Json.hpp"
#include "Utils/PauliStrings.hpp"

namespace tket {

/* Dictates whether synthesis of a PauliGraph should
    be done on the Paulis individually, making use of the pairwise
    interactions or collecting into mutually commuting sets. */
enum class PauliSynthStrat { Individual, Pairwise, Sets };

NLOHMANN_JSON_SERIALIZE_ENUM(
    PauliSynthStrat, {{PauliSynthStrat::Individual, "Individual"},
                      {PauliSynthStrat::Pairwise, "Pairwise"},
                      {PauliSynthStrat::Sets, "Sets"}});

class Transform {
 public:
  typedef std::function<bool(Circuit&)> Transformation;
  typedef std::function<unsigned(const Circuit&)> Metric;

  // the actual transformation to be applied
  // performs transformation in place and returns true iff made some change
  Transformation apply;  // this would ideally be `const`, but that deletes the
                         // copy assignment operator for Transform.

  explicit Transform(const Transformation& trans) : apply(trans) {}

  static const Transform id;  // identity Transform (does nothing to Circuit)

  friend Transform operator>>(const Transform& lhs, const Transform& rhs);

  /**
   * Squash sequences of 3-qubit instructions into their canonical 20-CX form.
   *
   * The circuit should comprise only CX and single-qubit gates. The transform
   * may perform a combination of 2-qubit (KAK) and 3-qubit decompositions of
   * subcircuits, but only does so if this reduces the CX count.
   *
   * @return Transform implementing the squash
   */
  static Transform three_qubit_squash();

  ////////////////////////
  // Contextual Reduction//
  ////////////////////////

  /**
   * Remove all operations that have no @ref OpType::Output or
   * @ref OpType::ClOutput in their causal future.
   */
  static Transform remove_discarded_ops();

  /** Allow insertion of classical operations when simplifying? */
  enum class AllowClassical { Yes, No };

  /** Automatically create all qubits in zero state when simplifying? */
  enum class CreateAllQubits { Yes, No };

  /**
   * Simplify the circuit where it acts on known basis states.
   *
   * Whenever a gate transforms a known basis state to another known basis
   * state, remove it, inserting X gates where necessary to achieve the same
   * state. Beginning with Create and Reset vertices, move forward through the
   * circuit as far as we can in this way.
   *
   * Global phase is ignored.
   *
   * @param allow_classical if allowed, insert SetBits operations in place of
   *  Measure operations that act on known basis states
   * @param create_all_qubits if enabled, annotate all qubits as initialized
   *  to zero as part of the transform, before applying simplification
   * @param xcirc 1-qubit circuit implementing an X gate (if null, an X gate
   *  is used)
   */
  static Transform simplify_initial(
      AllowClassical allow_classical = AllowClassical::Yes,
      CreateAllQubits create_all_qubits = CreateAllQubits::No,
      std::shared_ptr<const Circuit> = 0);

  /**
   * Commute classical maps through measurements.
   *
   * A "classical map" is a pure quantum operation that acts as a permutation
   * of the computational basis states composed with a diagonal operator.
   *
   * A "measured classical map" is a classical map whose succeeding vertices
   * in the DAG are all @ref OpType::Measure.
   *
   * This transform replaces all measured classical maps that are followed by
   * Measure operations whose quantum output is discarded with classical
   * operations following the Measure.
   *
   * The process is repeated until no such replacements are possible.
   *
   * Global phase is not preserved.
   */
  static Transform simplify_measured();
};

inline const Transform Transform::id = Transform([](const Circuit&) {
  return false;
});  // returns `false` as it does not change the Circuit in any way

}  // namespace tket
