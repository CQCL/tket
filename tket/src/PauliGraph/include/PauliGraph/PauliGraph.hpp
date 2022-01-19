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

#include <fstream>

#include "Clifford/CliffTableau.hpp"
#include "Utils/Expression.hpp"
#include "Utils/GraphHeaders.hpp"
#include "Utils/PauliStrings.hpp"
#include "Utils/SequencedContainers.hpp"

namespace tket {
class Gate;

struct PauliGadgetProperties {
  QubitPauliTensor tensor_;
  Expr angle_;
};

struct DependencyEdgeProperties {};

typedef boost::adjacency_list<
    boost::listS, boost::listS, boost::bidirectionalS, PauliGadgetProperties,
    DependencyEdgeProperties>
    PauliDAG;
typedef boost::graph_traits<PauliDAG>::vertex_descriptor PauliVert;
typedef boost::graph_traits<PauliDAG>::edge_descriptor PauliEdge;

typedef sequence_set_t<PauliVert> PauliVertSet;
typedef sequence_set_t<PauliEdge> PauliEdgeSet;
typedef std::list<std::pair<OpType, qubit_vector_t>> Conjugations;

class Circuit;

/**
 * Dependency graph of a circuit wrt Pauli Gadgets.
 * Constructed by effectively commuting all non-Clifford gates to the front
 * of the circuit and determining their dependencies based on commutation
 * of the Pauli strings.
 * The Clifford effect of a circuit is maintained as a tableau, thought of as
 * being applied after all of the gadgets.
 */
class PauliGraph {
 public:
  /** Construct an empty dependency graph for the identity over n qubits. */
  explicit PauliGraph(unsigned n);

  /** Construct an empty dependency graph for the identity over given qubits. */
  explicit PauliGraph(const qubit_vector_t &qbs, const bit_vector_t &bits = {});

  /**
   * Applies the given gate to the end of the circuit.
   * Clifford gates transform the tableau.
   * Non-Clifford gates are transformed into gadgets by the tableau and added
   * to the graph.
   */
  void apply_gate_at_end(const Gate &gate, const unit_vector_t &args);

  /** Visualisation of the dependency graph */
  void to_graphviz_file(const std::string &filename) const;
  void to_graphviz(std::ostream &out) const;

  const CliffTableau &get_clifford_ref() { return cliff_; }
  unsigned n_vertices() const { return boost::num_vertices(this->graph_); }

  friend PauliGraph circuit_to_pauli_graph(const Circuit &circ);
  friend Circuit pauli_graph_to_circuit_individually(
      const PauliGraph &pg, CXConfigType cx_config);
  friend Circuit pauli_graph_to_circuit_pairwise(
      const PauliGraph &pg, CXConfigType cx_config);
  friend Circuit pauli_graph_to_circuit_sets(
      const PauliGraph &pg, CXConfigType cx_config);

 private:
  /** The dependency graph of Pauli gadgets */
  PauliDAG graph_;
  /** The tableau of the Clifford effect of the circuit */
  CliffTableau cliff_;
  /** The record of measurements at the very end of the circuit */
  boost::bimap<Qubit, Bit> measures_;
  bit_vector_t bits_;
  /**
   * Caches of the set of Pauli gadgets that can be commuted to the start
   * or the end of the circuit.
   */
  PauliVertSet start_line_;
  PauliVertSet end_line_;

  PauliVertSet get_successors(const PauliVert &vert) const;
  PauliVertSet get_predecessors(const PauliVert &vert) const;
  PauliEdgeSet get_in_edges(const PauliVert &vert) const;
  PauliEdgeSet get_out_edges(const PauliVert &vert) const;
  PauliVert source(const PauliEdge &edge) const;
  PauliVert target(const PauliEdge &edge) const;

  /**
   * Appends a pauli gadget at the end of the dependency graph.
   * Assumes this is the result AFTER pushing it through the Clifford
   * tableau.
   */
  void apply_pauli_gadget_at_end(
      const QubitPauliTensor &pauli, const Expr &angle);

  /**
   * Iterates through the vertices of a PauliGraph in a topological ordering.
   * When there are multiple commuting vertices that could be emitted, this
   * selects the one with the lowest lexicographic ordering on the Pauli string.
   */
  class TopSortIterator {
   public:
    TopSortIterator();
    explicit TopSortIterator(const PauliGraph &pg);

    const PauliVert &operator*() const;
    const PauliVert *operator->() const;
    bool operator==(const TopSortIterator &other) const;
    bool operator!=(const TopSortIterator &other) const;

    TopSortIterator operator++(int);
    TopSortIterator &operator++();

   private:
    const PauliGraph *pg_;
    PauliVert current_vert_;
    std::set<std::pair<QubitPauliTensor, PauliVert>>
        search_set_;  // Use pair to force ordering by string
    std::unordered_set<PauliVert> visited_;
  };

  TopSortIterator begin() const;
  TopSortIterator end() const;
};

}  // namespace tket
