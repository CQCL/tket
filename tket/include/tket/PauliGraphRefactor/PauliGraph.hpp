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

#pragma once

#include "PGOp.hpp"
#include "tket/Utils/GraphHeaders.hpp"

namespace tket {
namespace pg {

/**
 * PauliGraph
 *
 * This data structure provides a balance between the simple rewriting of an
 * instruction graph (with arcs between operations sharing the same physical
 * resource, e.g. Circuit) and the abstraction of a dependency DAG (abstracts
 * away all commutations).
 *
 * We attribute each instruction to a small number of Pauli strings, with the
 * guarantee that if each string from A commutes with each string of B then A
 * and B commute (this is a safe under-approximation of commutativity - there
 * may be commutations this doesn't identify). Rewriting requires us to update
 * the Pauli strings and the relation of anticommutations between the strings.
 *
 * We separately use a true dependency DAG for the classical dependencies (i.e.
 * there is a single edge between two operations if reordering them would cause
 * a RAW, WAR, or WAW hazard).
 *
 * We intend to support the following rewrites during optimisation:
 * - Reordering commuting operations
 * - Pauli reorder rules (just updating phases of strings)
 * - Clifford reorder rules (updating Pauli strings by multiplication)
 * - Merging compatible vertices (rotations, measurements, discards, etc.)
 * - "Product Rotation Lemma" actions (multiplies a Pauli string by a
 * stabilizer; see Simmons 2021)
 * - Deletion of identity vertices
 * - Deletions of vertices at start and end
 * - Absorbing Cliffords into the start and end tableaux
 * - Changing vertex types (e.g. continuously-parameterised rotation to discrete
 * Clifford rotation, reset expansion)
 *
 * Each operation corresponds to exactly one node in the classical graph but may
 * use multiple Pauli strings, so we attach operation details to the vertices of
 * the classical graph. The heterogeneity of contents for different kinds of
 * operations encourages an object-oriented structure for node contents, similar
 * to Ops in Circuits. Unlike Ops, the large variability in Pauli strings means
 * we won't benefit significantly from reusing immutable objects, so we instead
 * store separate objects for each vertex and allow them to be mutable to update
 * in-place where possible.
 *
 * Few rewrites will update the classical data so maintaining the classical
 * dependency for fast lookup is best (as opposed to maintaining a candidate
 * temporal ordering of the operations and determining classical dependencies on
 * the fly). Dependencies are typically sparse, so a directed adjacency list is
 * suitable.
 *
 * Some additional lookup maps maintain the most recent reads and writes to each
 * classical Bit to aid vertex insertion. These will be largely unimportant when
 * it comes to rewriting though.
 *
 * We store the anticommutation between the Pauli strings of different
 * operations to save recalculating them a lot on the fly. We specifically store
 * a directed form of the anticommutation relation that also factors in the
 * ordering of the operations, i.e. (P, Q) means both P and Q anticommute and
 * P's operation occurs after Q's. This can be a relatively dense relation and
 * updates due to multiplying strings involve taking XOR or symmetric difference
 * between the ancestors/descendants, so we store it as a Binary matrix for easy
 * updating via row/column updates. Row i indicates the anticommuting ancestors
 * (earlier in the circuit) of Pauli i, and column i indicates the anticommuting
 * descendants (later in the circuit).
 *
 * During rewrites, once we have decided on a vertex to rewrite around, we will
 * need to both find the rows/columns in the anticommutation matrix
 * corresponding to a particular vertex. Often the entries in the matrix will
 * then inform which other vertices need to be rewritten, e.g. when moving a
 * Clifford instruction to the start of the circuit, the positive indices in its
 * row give the ancestors that need to be updated, so we also need a reverse
 * lookup from the table indices. It is easiest to maintain this mapping as a
 * multi-indexed container, allowing other data to also be attached if needed in
 * the future.
 *
 * Each Pauli string within the PauliGraph can be uniquely identified either by
 * its index in the anticommutation matrix, or by a combination of the vertex
 * and index of the PauliString within the PGOp, referred to as its port. The
 * number of ports and their ordering/interpretation is fixed based on the
 * PGOpType/subclass of PGOp.
 *
 * During rewrites which eliminate vertices, we leave unused rows/columns in the
 * anticommutation matrix rather than attempt to reduce it at every opportunity.
 * A cleanup method can be written if we wish to run this occasionally during
 * long rewrite procedures.
 *
 * Whilst previous iterations of PauliGraph contained an explicit Clifford
 * tableau at the start or end of the circuit, we choose to represent these
 * within the graph itself, since including them in the anticommutation matrix
 * allows for easy identification of opportunities for eliminating instructions
 * around discards or stabilizers, or applying PRL actions. In the case where we
 * need to relate Pauli strings to inputs or outputs, we follow the style of
 * ChoiMixedTableau in describing pairs of related Pauli strings over the inputs
 * and interior or over the interior and outputs. However, we only care about
 * the interior Pauli strings in the anticommutation matrix. If they are not
 * provided explicitly, they are assumed to be identity circuits.
 *
 * When a vertex may contain multiple ports, such as InputTableau and
 * OutputTableau, we view the actions on the ports as happening simultaneously,
 * so the anticommutation matrix will read false in the corresponding entries
 * even if the Pauli strings anticommute.
 */

typedef boost::adjacency_list<
    boost::setS, boost::listS, boost::bidirectionalS, PGOp_ptr>
    PGClassicalGraph;
typedef PGClassicalGraph::vertex_descriptor PGVert;

struct PGPauli {
  unsigned index;
  PGVert vert;
  unsigned port;
};

struct TagID {};
struct TagOp {};

typedef boost::multi_index::multi_index_container<
    PGPauli,
    boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
            boost::multi_index::tag<TagID>,
            boost::multi_index::member<PGPauli, unsigned, &PGPauli::index>>,
        boost::multi_index::hashed_non_unique<
            boost::multi_index::tag<TagOp>,
            boost::multi_index::member<PGPauli, PGVert, &PGPauli::vert>>>>
    PGIndex;

class PauliGraph {
 public:
  /**
   * Construct an empty PauliGraph with no Qubits or Bits.
   */
  explicit PauliGraph();

  /**
   * Construct an empty PauliGraph representing the identity over some defined
   * set of Qubits and Bits. This will initially lack any PGInputTableau or
   * PGOutputTableau, so these should be added explicitly if they wish to be
   * used.
   */
  explicit PauliGraph(
      const std::set<Qubit>& qubits, const std::set<Bit>& bits = {});

  /**
   * Get a reference to the set of all qubits used in the Circuit captured. Such
   * qubits may not be open boundaries, as they may be initialised and discarded
   * in the input and output tableaux.
   */
  const std::set<Qubit>& get_qubits() const;

  /**
   * Get a reference to the set of all classical bits used in the Circuit
   * captured.
   */
  const std::set<Bit>& get_bits() const;

  /**
   * Get the vertex of the unique PGInputTableau. If no such vertex exists, it
   * is interpreted as an identity process, and this method returns
   * std::nullopt.
   */
  std::optional<PGVert> get_input_tableau() const;

  /**
   * Get the vertex of the unique PGOutputTableau. If no such vertex exists, it
   * is interpreted as an identity process, and this method returns
   * std::nullopt.
   */
  std::optional<PGVert> get_output_tableau() const;

  /**
   * Given a PGVert within the PauliGraph, looks up the PGOp_ptr stored there.
   * This does not actively verify that the PGVert belongs to this PauliGraph
   * (errors such as segmentation faults may occur if misused). The PGOp_ptr is
   * a non-const shared pointer to the internal data, so even though this method
   * is marked const it is possible to update internal data of the PauliGraph by
   * modifying the PGOp through this pointer.
   */
  PGOp_ptr get_vertex_PGOp_ptr(const PGVert& v) const;

  /**
   * Writes a graphviz representation of the PauliGraph to a stream. Use this
   * for visualisation. Each vertex in the PauliGraph is represented as a
   * cluster of graphviz vertices (one per active Pauli). Classical dependencies
   * are drawn as edges between clusters and the anti-commutation dependencies
   * between Paulis are drawn as edges between the corresponding vertices.
   */
  void to_graphviz(std::ostream& out) const;

  /**
   * Inserts a new vertex at the end of the PauliGraph. Throws an exception if a
   * PGInitialTableau is inserted after other vertices or if any vertex is
   * inserted after a PGOutputTableau.
   */
  PGVert add_vertex_at_end(PGOp_ptr op);

  /**
   * Verification of validity of the data structure. This is computationally
   * expensive so it is intended for use in debugging and tests, but not live
   * code.
   */
  void verify() const;

  /**
   * Returns all PGOps in a valid topological sort of the diagram. The exact
   * order depends on the internal order of vertices in c_graph_.
   */
  std::list<PGOp_ptr> pgop_sequence() const;

  std::list<std::list<PGOp_ptr>> pgop_commuting_sets() const;

  /**
   * Symbolic substitution: replaces each PGOp in the PauliGraph according to
   * the substitution map sending some set of symbols (not necessarily the same
   * as free_symbols()) to some other expressions.
   */
  void symbol_substitution(const SymEngine::map_basic_basic& sub_map);

  // Set of all free symbols occurring in operation parameters
  SymSet free_symbols() const;

  // Whether the PauliGraph's operations contain any symbolic parameters
  bool is_symbolic() const;

 private:
  MatrixXb pauli_ac_;
  PGIndex pauli_index_;
  PGClassicalGraph c_graph_;
  std::set<Qubit> qubits_;
  std::set<Bit> bits_;

  // Helper variables for tracking previous reads from and writes to each bit to
  // simplify adding dependencies in add_vertex_at_end.
  std::map<Bit, PGVert> last_writes_;
  std::map<Bit, std::unordered_set<PGVert>> last_reads_;

  std::optional<PGVert> initial_tableau_;
  std::optional<PGVert> final_tableau_;

  /**
   * Replaces the QubitPauliString of row \p target_r with i^{ \p coeff } *
   * source * target and updates pauli_ac_ accordingly.
   */
  void multiply_strings(
      unsigned source_r, unsigned target_r, quarter_turns_t coeff = 0);
};

}  // namespace pg
}  // namespace tket
