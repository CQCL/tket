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

  /**
   *
   */
  explicit Placement(const Architecture& _architecture)
      : architecture_(_architecture) {}

  /**
   *
   */
  Placement(){};

  /**
   * Reassigns some UnitID in circ_ as UnitID in architecture_
   *
   * @param circ_ Circuit to be relabelled
   * @param compilation_map For tracking changes during compilation
   *
   * @return true iff circuit or maps are modified
   */
  bool place(
      Circuit& circ_,
      std::shared_ptr<unit_bimaps_t> compilation_map = nullptr) const;

  /**
   * Reassigns some UnitID in circ_ as UnitID in architecture_, according to
   * given map.
   *
   * @param circ_ Circuit to be relabelled
   * @param map_ relabelling
   * @param compilation_map For tracking changes during compilation
   *
   * @return true iff circuit or maps were modified
   */
  static bool place_with_map(
      Circuit& circ, std::map<Qubit, Node>& map_,
      std::shared_ptr<unit_bimaps_t> compilation_map = nullptr);

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
  std::map<Qubit, Node> get_placement_map(const Circuit& circ_) const;

  /**
   *
   * For some Circuit, returns maps between Circuit UnitID and
   * Architecture UnitID that can be used for reassigning UnitID in
   * Circuit. Maps expected to give similiar performance for given method.
   * For Placement this naively assigns every Qubit to some Node.
   *
   * @param circ_ Circuit relabelling map is constructed from
   *
   * @return Map between Circuit and Architecture UnitID
   */
  virtual std::vector<std::map<Qubit, Node>> get_all_placement_maps(
      const Circuit& circ_) const;

  /**
   * Returns a reference to held Architecture.
   * Used to know Architecture properties to set predicates
   * during compilation.
   *
   * @return Architecture
   */
  const Architecture& get_architecture_ref() { return architecture_; }

  virtual ~Placement(){};

 protected:
  Architecture architecture_;
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

  /**
   * Holds information for slice wise iteration of Circuit
   * @param _circ Circuit to iterate through
   */
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

  static const std::vector<WeightedEdge> default_pattern_weighting(
      const Circuit& circuit);
  static const std::vector<WeightedEdge> default_target_weighting(
      Architecture& passed_architecture);

  explicit GraphPlacement(
      const Architecture& _architecture, unsigned _maximum_matches = 10000000,
      unsigned _timeout = 100,
      const std::function<std::vector<WeightedEdge>(const Circuit&)>
          _weight_pattern_graph = default_pattern_weighting,
      const std::function<std::vector<WeightedEdge>(Architecture&)>
          _weight_target_graph = default_target_weighting)
      : weight_pattern_graph_(_weight_pattern_graph),
        weight_target_graph_(_weight_target_graph),
        maximum_matches_(_maximum_matches),
        timeout_(_timeout) {
    // TODO: weight_target_graph is not const as it caches all distances
    // This is arguably beneficial as this will be done at some point
    // However it means we need to copy _architecture before it's used.
    // And additionally as we are finding new weights we don't just mutate
    // architecture_ but assign a new object to it.
    // Meaning, this logic seems contrived but maybe it's fine? A second opinion
    // is welcomed. We're going to copy _architecture either way.
    architecture_ = _architecture;
    this->weighted_target_edges =
        this->weight_target_graph_(architecture_);
    // architecture_ = this->construct_target_graph(weighted_target_edges);
  }

  /**
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
  std::vector<std::map<Qubit, Node>> get_all_placement_maps(
      const Circuit& circ_) const override;

  /**
   * @return maximum matches found during placement
   */
  unsigned get_maximum_matches() const { return this->maximum_matches_; }

  /**
   * @return maximum time (ms)
   */
  unsigned get_timeout() const { return this->timeout_; }

 protected:
  const std::function<std::vector<WeightedEdge>(const Circuit&)>
      weight_pattern_graph_;
  const std::function<std::vector<WeightedEdge>(Architecture&)>
      weight_target_graph_;
  unsigned maximum_matches_;
  unsigned timeout_;
  std::vector<WeightedEdge> weighted_target_edges;

  QubitGraph construct_pattern_graph(
      const std::vector<WeightedEdge>& edges) const;

  Architecture construct_target_graph(
      const std::vector<WeightedEdge>& edges) const;
};

/** Solves the pure unweighted subgraph monomorphism problem, trying
 * to embed the pattern graph into the target graph.
 * Note that graph edge weights are IGNORED by this function.
 */
std::vector<boost::bimap<Qubit, Node>> get_weighted_subgraph_monomorphisms(
QubitGraph::UndirectedConnGraph& pattern_graph,
    Architecture::UndirectedConnGraph& target_graph, unsigned max_matches,
    unsigned timeout_ms);

void to_json(nlohmann::json& j, const Placement::Ptr& placement_ptr);
void from_json(const nlohmann::json& j, Placement::Ptr& placement_ptr);

}  // namespace tket