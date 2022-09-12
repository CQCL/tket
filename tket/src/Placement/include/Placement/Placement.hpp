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

#include "Architecture/Architecture.hpp"
#include "Circuit/Circuit.hpp"
#include "Graphs/QubitGraph.hpp"

namespace tket {

class Placement {
 public:
  typedef std::shared_ptr<Placement> Ptr;

  struct Frontier {
    // set of 2qb vertices which need to be solved for
    std::shared_ptr<Slice> slice;
    // Quantum Edges coming in to vertices in slice, indexed by qubit
    std::shared_ptr<unit_frontier_t> quantum_in_edges;
    // Quantum Edges leaving vertices in slice, indexed by qubit
    std::shared_ptr<unit_frontier_t> quantum_out_edges;
    // Boolean edges coming in to vertices in slice. Guarantees that all edges
    // into every vertex in slice is represented in next_cut
    std::shared_ptr<b_frontier_t> boolean_in_edges;

    // reference to circuit that it acts on
    const Circuit& circ;

    // initialise at front of circuit
    explicit Frontier(const Circuit& _circ);
    // move to next slice
    void next_slicefrontier();
  };

  /**
   *
   */
  explicit Placement(const Architecture& _arc) : arc_(_arc) {}

  /**
   *
   */
  Placement(){};

  /**
   * Reassigns some UnitID in circ_ as UnitID arc_
   *
   * @param circ_ Circuit to be relabelled
   *
   */
  bool place(Circuit& circ_) const;

  /**
   *
   * For some Circuit, returns a map between Circuit UnitID and
   * Architecture UnitID that can be used for reassigning UnitID in
   * Circuit. Map is expected to give best performance for given method.
   *
   * @param circ_ Circuit relabelling map is constructed from
   *
   * @return Map between Circuit and Architecture UnitID
   */
  std::map<UnitID, UnitID> get_placement_map(const Circuit& circ_) const;

  /**
   *
   * For some Circuit, returns maps between Circuit UnitID and
   * Architecture UnitID that can be used for reassigning UnitID in
   * Circuit. Maps expected to give similiar performance for given method.
   *
   * @param circ_ Circuit relabelling map is constructed from
   *
   * @return Map between Circuit and Architecture UnitID
   */
  virtual std::vector<std::map<UnitID, UnitID>> get_all_placement_maps(
      const Circuit& circ_) const;

 protected:
  Architecture arc_;
};

JSON_DECL(Placement::Ptr);

class GraphPlacement : public Placement {
 public:
  /**
   * Holds information for constructing a weighted edge in a QubitGraph.
   * @param node0 UnitID for first node in edge
   * @param node1 UnitID for second node in edge
   * @param weight Unsigned giving a weight for implied edge
   */
  struct WeightedEdge {
    UnitID node0;
    UnitID node1;
    unsigned weight;
  };

  static const std::vector<WeightedEdge> default_weighting(
      const Circuit& circuit);

  explicit GraphPlacement(
      const Architecture& _arc, unsigned _maximum_matches, unsigned _timeout,
      const std::function<std::vector<WeightedEdge>(const Circuit&)>
          _weight_pattern_graph = default_weighting)
      : weight_pattern_graph_(_weight_pattern_graph),
        maximum_matches_(_maximum_matches),
        timeout_(_timeout) {
    arc_ = _arc;
  }

  /**
   *
   * For some Circuit, returns maps between Circuit UnitID and
   * Architecture UnitID that can be used for reassigning UnitID in
   * Circuit. Maps are constructed by running a Weighted Subgraph Monomorphism
   * for the given problem and returning maximum_matches_ number of
   * potential solutions, ranked.
   *
   * @param circ_ Circuit relabelling map is constructed from
   *
   * @return Map between Circuit and Architecture UnitID
   */
  std::vector<std::map<UnitID, UnitID>> get_all_placement_maps(
      const Circuit& circ_) const override;

 protected:
  const std::function<std::vector<WeightedEdge>(const Circuit&)>
      weight_pattern_graph_;
  unsigned maximum_matches_;
  unsigned timeout_;

  QubitGraph construct_pattern_graph(
      const std::vector<WeightedEdge>& edges) const;

  /** Solves the pure unweighted subgraph monomorphism problem, trying
   * to embed the pattern graph into the target graph.
   * Note that graph edge weights are IGNORED by this function.
   */
  std::vector<boost::bimap<Qubit, Node>> get_weighted_subgraph_monomorphisms(
      const QubitGraph::UndirectedConnGraph& pattern_graph,
      const Architecture::UndirectedConnGraph& target_graph,
      unsigned max_matches, unsigned timeout_ms) const;
};

}  // namespace tket