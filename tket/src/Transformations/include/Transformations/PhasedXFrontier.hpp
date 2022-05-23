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

#include <cmath>
#include <set>

#include "Circuit/Circuit.hpp"
#include "Ops/OpPtr.hpp"
#include "SingleQubitSquash.hpp"

namespace tket {

namespace Transforms {

// type aliases
using OptVertex = std::optional<Vertex>;
using OptVertexVec = std::vector<OptVertex>;
using OptEdge = std::optional<Edge>;
using OptEdgeVec = std::vector<OptEdge>;

/**
 * @brief Frontier and circuit manipulation for `globalise_PhasedX`
 *
 * This helper class for the \ref globalise_PhasedX transform stores the circuit
 * frontier as useful for that specific pass. It also provides a simple
 * interface for all the circuit manipulations that are required within that
 * pass.
 *
 * The transform itself therefore only needs to worry about the logic of the
 * transformation and none of the bookkeeping and circuit substitutions.
 *
 * All gates of `circ` must be either PhasedX, NPhasedX, Rz or multi-qubit
 * gates.
 *
 * ## Intervals
 * For each qubit, the frontier stores the current interval. An interval is a
 * single-qubit subcircuit of `circ` forming a sequence of single-qubit gates
 * between two multi-qubit gates (or non-gate vertices) of the original
 * circuit. Each interval is defined by start and end edges, which are edges of
 * the circuit DAG ("end" is in the future of "start"). The source of the start
 * edge and the target of the end edge are either multi-qubit gates or non-gate
 * vertices such as input/output vertices, create/discard vertices etc.
 *
 * ## Initialisation and moving forward
 * The frontier is initialised at the beginning of the circuit, at which stage
 * an interval is created for each qubit of `circ`.
 * The frontier then moves forward in jumps, where the new start of the interval
 * will be the edge after the former end edge. The new end edge will be moved
 * forward along the qubit until its target is a multi-qubit gate or the Output
 * vertex.
 *
 * One can move forward an individual interval by calling `next_interval`.
 * Careful, however: all intervals are assumed to be space-like separated (i.e.
 * there cannot be a path between an interval and another at any time).
 * Moving arbitrary intervals forward may break this invariant.
 * To avoid this issue, the `next_multiqb` method is provided. It takes a multi-
 * qb vertex `v` as argument and will move forward all the qubits whose end edge
 * has target `v`. By calling `next_multiqb` on each multi-qubit gate in a
 * causal order, the space-like separation invariant will be maintained.
 *
 * Within each interval, edges pointing to PhasedX (or Rx) gates are called
 * beta edges: the angle parameter of these gates correspond to the β angle for
 * the next global PhasedX(β, α) gate. When squashed to normal form, each
 * interval will have at most one beta edge.
 * For non-normalised intervals, a beta edge can also be "shadowed" if the gate
 * it points to is a multi-qubit NPhasedX gate that is not entirely in the
 * current frontier, i.e. there is at least one qubit that needs to process
 * another gate before getting to the NPhasedX gate.
 * We cannot transform that gate until all other gates in its causal past
 * have been processed, and we therefore say ignore that beta edge.
 *
 * ## Operations on intervals
 * ### Squashing
 * Each interval can be squashed into a normal form `PhasedX + Rz` using
 * `squash_interval`, or `squash_intervals` to squash all current intervals.
 * It is recommended to do this before acting on the intervals, as otherwise
 * multiple single-qb operations can accumulate and lead to inefficient
 * decompositions. This is what `globalise_PhasedX(squash=true)` does
 * (the default). This is however not necessary for the frontier to be well-
 * behaved
 *
 * ### Get beta angles
 * When squashed, each interval is associated with a unique "beta angle", given
 * by the formula U_2 = Rz(α)Rx(β)Rz(ɣ) = PhasedX(β, α)Rz(α + ɣ).
 * These can be obtained for each interval using `get_all_betas`, while
 * `get_beta_edge` and `get_beta_edges` will return the edge pointing to the
 * gate that contains the β-angle (ie the PhasedX gate).
 *
 * ### Insert global NPhasedX gates
 * `insert_1_phasedx` and `insert_2_phasedx` will insert 1x resp 2x global
 * NPhasedX gates at the beginning of the current frontier. The resuling
 * circuit will always be equivalent to the original one.
 */
class PhasedXFrontier {
 public:
  /**
   * @brief Construct a new Phased X Frontier object.
   *
   * @param circ The circuit to traverse.
   */
  PhasedXFrontier(Circuit& circ);

  /**
   * @brief Squash the current intervals on each qubit
   */
  void squash_intervals();

  /**
   * @brief Move frontier forward on qubit `i`.
   *
   * If there is no further interval on qubit `i`, the interval will be the
   * empty interval `(last_e, last_e)` where `last_e` is the last edge on qubit
   * `i`. Calling `next_interval` again on `i` will not do anything.
   *
   * @param i The qubit index whose interval will be moved.
   */
  void next_interval(unsigned i);

  /**
   * @brief Move frontier on all qubits with end target `v`.
   *
   * This is the recommended way to move the frontier forward. Make successive
   * calls to this method by looping through all multi-qubit gates of `circ` in
   * topological ordering.
   *
   * @param v A multi-qubit gate in `circ`.
   */
  void next_multiqb(const Vertex& v);

  /**
   * @brief Insert a single global NPhasedX at the current frontier.
   *
   * This will insert a NPhasedX gate with β-angle given by the β-angle of qubit
   * `i`. This will get rid of the PhasedX on qubit `i` (and any other qubits
   * with the same β-angle), but might introduce new garbage on the remaining
   * qubits.
   *
   * @param i The qubit index whose PhasedX will be replaced by a global gate.
   */
  void insert_1_phasedx(unsigned i);

  /**
   * @brief Insert two global NPhasedX gates at the current frontier.
   *
   * This will replace all current PhasedX gates with two global NPhasedX gates
   * using the equation
   *
   *    PhX(α, β) = Rz(β + 1/2) PhX(1/2, 0) Rz(α) PhX(-1/2, 0) Rz(-β - 1/2)
   *
   * ie the PhX gates can be made global.
   */
  void insert_2_phasedx();

  /**
   * @brief The qubits whose end edge's target is `v`.
   *
   * @param v  A multi-qubit gate in `circ`.
   * @return std::set<unsigned> Indices of the qubits whose end target is `v`.
   */
  std::set<unsigned> qubits_ending_in(const Vertex& v) const;

  // =============== beta edges, beta vertices and beta angles ===============
  // Any SU(2) operation is written as Rz(α)Rx(β)Rz(ɣ). We care mostly about
  // the β-angle as we can always add Rz gates left and right to fix α and ɣ.
  // =========================================================================

  /**
   * @brief Get all β-angles.
   *
   * If a qubit has no β-angle gate, its β-angle is 0. Note that some gates
   * might be "shadowed" and thus appear as 0 even though a β-angle gate exists
   * (see `get_all_beta_edges`)l
   *
   * @return std::vector<Expr> The vector of β-angles.
   */
  std::vector<Expr> get_all_betas() const;

  /**
   * @brief Get vertices of all β-angle gates.
   *
   * If a qubit has no β-angle gate (or it is shadowed), return std::nullopt.
   *
   * @return OptVertexVec The vertices of all β-angle gates.
   */
  OptVertexVec get_all_beta_vertices() const;

  /**
   * @brief Get the edge with target the β-angle gate on qubit `i`.

   * If qubit `i` has no β-angle gate, return std::nullopt.
   *
   * @param i The qubit index.
   * @return OptEdge The edge with target the current β-angle gate.
   */
  OptEdge get_beta_edge(unsigned i) const;

  /**
   * @brief Get edges with target the β-angle gates.
   *
   * An edge will be std::nullopt if there is no PhasedX gate on that qubit
   * or if it is shadowed. Shadowed gates are NPhasedX gates that are not
   * entirely in the current frontier, i.e. there is at least one qubit that
   * needs to process another gate before getting to the NPhasedX gate.
   * We cannot transform that gate until all other gates in its causal past
   * have been processed, and we therefore ignore it.
   *
   * To get around this issue, we recommend using the `decompose_nphasedx`
   * transform first, so that there are no NPhasedX gates left in the input
   * circuit. This of course is not possible if one wants to preserve some of
   * the existing structure of the circuit.
   *
   * @return OptEdgeVec The edges with target the β-angle gates.
   */
  OptEdgeVec get_all_beta_edges() const;

  /**
   * @brief Whether there are further non-trivial intervals ahead.
   *
   * Useful to know if a certain computation can be deferred until later.
   */
  bool are_phasedx_left() const;

  /**
   * @brief Shorten current intervals by moving the start edge past `n` global
   * gates.
   *
   * @param n The number of global gates to move past.
   */
  void skip_global_gates(unsigned n = 1);

  /**
   * @brief Whether the frontier has reached the end of the circuit.
   */
  bool is_finished();

  /**
   * @brief Whether the vertex is within a single-qb interval
   *
   * Checks if v is a multi-qubit gate or a final optype.
   *
   * @param op The operation.
   */
  static bool is_interval_boundary(Op_ptr op);

 private:
  // squash gates in the current interval on qubit i to ZXZ form
  void squash_interval(unsigned i);

  // forwards to public static method
  bool is_interval_boundary(Vertex v) const;

  // given an edge within an interval (eg the first), advances through the
  // circuit until hitting a interval boundary vertex
  Edge get_interval_end(Edge e) const;

  // given the last edge of an interval, returns the beginning of the next
  // interval
  Edge get_interval_start(Edge e) const;

  struct BackupIntervals {
    std::vector<VertPort> start;
    std::vector<VertPort> end;
  };

  // backup current frontier, so that you can insert at the current hole
  BackupIntervals backup_intervals() const;

  // using backup, restore old frontier
  void restore_intervals(const BackupIntervals& b);

  // for each qubit: first and last edge of current interval
  std::vector<std::pair<Edge, Edge>> intervals_;

  // a reference to the circuit
  Circuit& circ_;

  // squasher to squash gates
  SingleQubitSquash squasher_;
};

/**
 * @brief Whether all elements of `vec` are std::nullopt.
 */
bool all_nullopt(const OptVertexVec& vec);

}  // namespace Transforms

}  // namespace tket
