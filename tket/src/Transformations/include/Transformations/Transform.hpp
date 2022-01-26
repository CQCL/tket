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

  /////////////////////////
  // Clifford Optimisation//
  /////////////////////////

  // These apply some of the Clifford rules in the paper "Optimising Clifford
  // Circuits with Quantomatic" All of these expect and produce CX, Z, X, S, V,
  // and other single qubit gates

  static Transform singleq_clifford_sweep();

  // mutli-qubit patterns that decrease the CX count
  // inserting swaps can sometimes cause errors elsewhere (e.g. routing), so
  // they can be turned off
  static Transform multiq_clifford_replacement(bool allow_swaps = false);

  static Transform clifford_reduction(bool allow_swaps = false);

  // copies Z through the target of a CX and X through the control
  static Transform copy_pi_through_CX();

  //////////////////////
  // Basic Optimisation//
  //////////////////////

  // removes gate-inverse pairs, merges rotations, removes identity rotations,
  // and removes redundant gates before measure Expects: Any gates Produces: The
  // same gate set
  static Transform remove_redundancies();

  /**
   * general u_squash by converting any chains of p, q gates (p, q in
   * {Rx,Ry,Rz}) to triples -p-q-p- or -q-p-q-
   *
   * if strict = false (default), then the chains will be squashed to either
   * -p-q-p- or -q-p-q-, and the third rotation will be commuted through the
   * next gate whenever possible if strict = true, then all chains will be
   * squashed to -p-q-p-
   *
   * Expects: p, q, and any multi-qubit gates
   * Produces: p, q, and any multi-qubit gates
   */
  static Transform squash_1qb_to_pqp(
      const OpType& q, const OpType& p, bool strict = false);

  // 1qb squashing into -Rz-Rx-Rz- or -Rx-Rz-Rx- form
  // Expects: Rx, Rz, and any multi-qubit gates
  // Produces: Rx, Rz, and any multi-qubit gates
  static Transform reduce_XZ_chains();

  /**
   * Squash all single-qubit gates to TK1.
   */
  static Transform squash_1qb_to_tk1();

  // identifies single-qubit chains and squashes them in the target gate set
  // Expects: any gates
  // Produces: singleqs and any multi-qubit gates
  static Transform squash_factory(
      const OpTypeSet& singleqs,
      const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
          tk1_replacement);

  // moves single qubit operations past multiqubit operations they commute with,
  // towards front of circuit (hardcoded)
  // Expects: Any gates
  // Produces: Any gates
  static Transform commute_through_multis();

  // replaces every non-global PhasedX gate with two global PhasedX gates,
  // using:
  //
  //         PhX(α, β) = PhX(-1/2, β + 1/2) Rz(α) PhX(1/2, β + 1/2)
  //
  // Setting α = 0 on non-targeted qubits makes the RHS the identity
  static Transform globalise_phasedx();

  // commutes Rz gates through ZZMax, and combines adjacent ZZMax gates
  // Expects: ZZMax, Rz, Rx
  // Produces: ZZMax, Rz, Rx
  static Transform commute_and_combine_HQS2();

  // commutes single qubit gates through SWAP gates, leaving them on the
  // PhysicalQubit with best fidelity for the given OP Expects: any single qubit
  // gates, SWAP gates Produces: any single qubit gates, SWAP gates
  static Transform commute_SQ_gates_through_SWAPS(
      const avg_node_errors_t& node_errors);
  static Transform commute_SQ_gates_through_SWAPS(
      const op_node_errors_t& node_errors);

  // squashes each sequence of two-qubit instructions into their canonical 3-CX
  // form (Cartan/KAK decomposition) The optional CX gate fidelity cx_fidelity
  // is used to produce circuit decompositions that maximise expected circuit
  // fidelity, assuming perfect local operations (i.e. with fidelity 1.).
  // Expects: CX and any single qubit gates
  // Produces: CX and single TK1 qubit gates
  static Transform two_qubit_squash(double cx_fidelity = 1.);

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

  // helper class subcircuits representing 2qb interactions
  struct Interaction {
    Interaction(const Qubit& _q0, const Qubit& _q1) : q0(_q0), q1(_q1) {}
    Qubit q0;  // Qubit numbers
    Qubit q1;
    Edge e0;  // In edges starting interaction
    Edge e1;
    unsigned count;      // Number of two qubit gates in interaction
    VertexSet vertices;  // Vertices in interaction subcircuit
  };
};

inline const Transform Transform::id = Transform([](const Circuit&) {
  return false;
});  // returns `false` as it does not change the Circuit in any way

}  // namespace tket
