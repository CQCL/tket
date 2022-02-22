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

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Architecture/Architecture.hpp"
#include "Characterisation/DeviceCharacterisation.hpp"
#include "Circuit/Circuit.hpp"
#include "Graphs/Utils.hpp"
#include "Utils/BiMapHeaders.hpp"
#include "Utils/GraphHeaders.hpp"
#include "Utils/Json.hpp"

namespace tket {

extern template class graphs::DirectedGraphBase<Qubit>;
extern template class graphs::DirectedGraph<Qubit>;

struct QubitWeight;
struct InteractionWeight;
class Placement;

typedef std::map<Qubit, Node> qubit_mapping_t;
typedef boost::bimap<Qubit, Node> qubit_bimap_t;
typedef qubit_bimap_t::left_map::const_iterator l_const_iterator_t;
typedef qubit_bimap_t::right_map::const_iterator r_const_iterator_t;
// Adjacent elements in a QubitLine interact in some timesteps such that
// these qubits do not need to be moved to be executed
typedef qubit_vector_t QubitLine;
typedef std::vector<QubitLine>
    QubitLineList;  // Used in placement of qubits methods

typedef std::shared_ptr<Placement> PlacementPtr;

JSON_DECL(PlacementPtr)

class QubitGraphInvalidity : public std::logic_error {
 public:
  explicit QubitGraphInvalidity(const std::string& message)
      : std::logic_error(message) {}
};

class PlacementError : public std::logic_error {
 public:
  explicit PlacementError(const std::string& message)
      : std::logic_error(message) {}
};

/** Print a map to stdout. */
template <class MapType>
void print_map(const MapType& m) {
  typedef typename MapType::const_iterator const_iterator;
  for (const_iterator iter = m.begin(), iend = m.end(); iter != iend; ++iter) {
    std::cout << iter->first << "-->" << iter->second << std::endl;
  }
}

// structure of configuration parameters for placement
struct PlacementConfig {
  // circuit look ahead limit
  unsigned depth_limit;
  // max edges in interaction graph
  unsigned max_interaction_edges;
  // max number of matches from monomorphism calculator
  unsigned vf2_max_matches = 1000;
  /*
  value of num_gates/num_qubits above which to contract architecture before
  placement for high values of this ratio it is assumed swap count is more
  critical than initial noise minimisation for which is architecture contraction
  to most mst highly connected subgraph is critical
  */
  unsigned arc_contraction_ratio = 10;
  // Timeout, corresponds to milliseconds. Default 1 minute.
  unsigned timeout = 60000;

  PlacementConfig(){};

  PlacementConfig(
      unsigned _depth_limit, unsigned _max_interaction_edges,
      unsigned _vf2_max_matches = 1000, unsigned _arc_contraction_ratio = 10,
      unsigned _timeout = 60000);

  bool operator==(const PlacementConfig& other) const;
};

JSON_DECL(PlacementConfig)

// stores and tracks the points of the circuit up to which has been solved
struct PlacementFrontier {
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
  explicit PlacementFrontier(const Circuit& _circ);
  // move to next slice
  void next_slicefrontier();
};

// Class for storing interaction graph.
// Interacting qubits have an edge between them.
class QubitGraph : public graphs::DirectedGraph<Qubit> {
 private:
  using Base = graphs::DirectedGraph<Qubit>;

 public:
  QubitGraph() : Base() {}
  explicit QubitGraph(const qubit_vector_t& _qubits) : Base(_qubits) {}
};

/* ACTUALLY PLACEMENT METHODS */

// generate interaction graph of circuit
QubitGraph generate_interaction_graph(
    const Circuit& circ, unsigned depth_limit = 10);
// generate lines of interacting qubits
QubitLineList qubit_lines(const Circuit& circ);
// generate mapping of qubit lines to lines on architecture. n_qubits is total
// number of qubits in circuit
qubit_mapping_t lines_on_arc(
    Architecture arc, QubitLineList qb_lines, unsigned n_qubits);

// build interaction graph of circ, max_edges in graph, depth_limit sets how far
// to look ahead
QubitGraph monomorph_interaction_graph(
    const Circuit& circ, const unsigned max_edges, unsigned depth_limit);

/**
 * Search for embeddings of the qubit graph in the architecture graph.
 *
 * @param arc architecture
 * @param q_graph qubit graph
 * @param max_matches maximum number of matches to find
 * @param timeout timeout in milliseconds
 *
 * @return vector of matches found, sorted in canonical order
 */
std::vector<qubit_bimap_t> monomorphism_edge_break(
    const Architecture& arc, const QubitGraph& q_graph, unsigned max_matches,
    unsigned timeout);

node_set_t best_nodes(Architecture& arc, unsigned n_remove);

class PatternError : public std::logic_error {
 public:
  explicit PatternError(const std::string& message)
      : std::logic_error(message) {}
};

// Note: boost's VF2 standard subgraph monomorphism uses the `always_equivalent`
// predicate, so it will match any vertex with any vertex
// and any edge with any edge, regardless of their bundled properties. This is
// very good for us!
struct QubitWeight {
  QubitWeight() : val(0.) {}
  explicit QubitWeight(const boost::no_property) : val(0.) {}
  explicit QubitWeight(double d) : val(d) {}
  double val;
};

struct InteractionWeight {
  InteractionWeight() : val(0.) {}
  explicit InteractionWeight(const boost::no_property) : val(0.) {}
  explicit InteractionWeight(double d) : val(d) {}
  template <typename Property>
  explicit InteractionWeight(Property p)
      : val(static_cast<double>(p.m_value)) {}
  double val;
};

// structure of qubit mapping with associated cost
struct MapCost {
  qubit_mapping_t map;
  double cost;
  bool operator<(const MapCost& other) const { return this->cost < other.cost; }
  bool operator>(const MapCost& other) const { return this->cost > other.cost; }
};

template <typename GraphP, typename GraphT>
class vf2_match_add_callback {
  using qubit_grapht_bimap_t =
      boost::bimap<Node, graphs::utils::vertex<GraphT>>;
  using qubit_graphp_bimap_t =
      boost::bimap<Qubit, graphs::utils::vertex<GraphP>>;

 public:
  vf2_match_add_callback(
      std::vector<qubit_bimap_t>& all_maps, const GraphP& pattern_graph,
      const GraphT& target_graph, unsigned _max)
      : n_maps_(all_maps),
        max(_max),
        pattern_graph_(pattern_graph),
        target_graph_(target_graph) {}

  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(const CorrespondenceMap1To2& f, const CorrespondenceMap2To1&);

  std::vector<qubit_bimap_t>& n_maps_;
  const unsigned max;

 private:
  const GraphP& pattern_graph_;
  const GraphT& target_graph_;
};

class Placement {
 public:
  explicit Placement(const Architecture& _arc) : arc_(_arc) {}
  Placement(){};

  /**
   * Modify qubits in place.
   *
   * @return true iff circuit or maps are modified
   */
  bool place(
      Circuit& circ_, std::shared_ptr<unit_bimaps_t> maps = nullptr) const;

  /**
   * Relabel circuit qubits to device nodes according to given map.
   *
   * @return true iff circuit or maps were modified
   */
  static bool place_with_map(
      Circuit& circ, qubit_mapping_t& map_,
      std::shared_ptr<unit_bimaps_t> maps = nullptr);

  virtual qubit_mapping_t get_placement_map(const Circuit& circ_) const;

  // methods that return maps, for base this returns one empty map
  virtual std::vector<qubit_mapping_t> get_all_placement_maps(
      const Circuit& circ_) const;

  static const std::string& unplaced_reg();

  const Architecture& get_architecture_ref() { return arc_; }
  virtual ~Placement(){};

 protected:
  Architecture arc_;
};

/**
 * NaivePlacement class provides methods for relabelling any
 * Qubit objects in some Circuit to Node objects in some Architecture
 * given the constraint that only Qubit that are not already labelled
 * as some Node can be relabelled, and only to Architecture Node
 * that are not already in the Circuit.
 */
class NaivePlacement : public Placement {
 public:
  /**
   * @param _arc Architecture object later relabellings are produced for
   */
  explicit NaivePlacement(const Architecture& _arc) { arc_ = _arc; }
  /**
   * Given some circuit, returns a map between Qubit which defines some
   * relabelling of some Circuit qubits to Architecture qubits
   *
   * @param circ_ Circuit map relabelling is defined for
   *
   * @return Map defining relabelling for circuit Qubit objects
   */
  qubit_mapping_t get_placement_map(const Circuit& circ_) const override;

  /**
   * Given some circuit, returns a single map for relabelling
   * in a vector.
   *
   * @param circ_ Circuit map relabelling is defined for
   *
   * @return Vector of a single Map defining relabelling for Circuit
   * Qubit objects.
   */
  std::vector<qubit_mapping_t> get_all_placement_maps(
      const Circuit& circ_) const override;
};

class LinePlacement : public Placement {
 public:
  explicit LinePlacement(const Architecture& _arc) { arc_ = _arc; }

  qubit_mapping_t get_placement_map(const Circuit& circ_) const override;

  // methods that return maps, for base this returns one empty map
  std::vector<qubit_mapping_t> get_all_placement_maps(
      const Circuit& circ_) const override;
};

class GraphPlacement : public Placement {
 public:
  explicit GraphPlacement(const Architecture& _arc) {
    arc_ = _arc;
    config_.depth_limit = 5;
    config_.max_interaction_edges = arc_.n_connections();
    config_.vf2_max_matches = 10000;
    config_.arc_contraction_ratio = 10;
  }

  explicit GraphPlacement(
      const Architecture& _arc, const PlacementConfig& _config)
      : Placement(_arc), config_(_config) {}

  explicit GraphPlacement(const PlacementConfig& _config) : config_(_config) {}

  PlacementConfig get_config() { return config_; }
  void set_config(const PlacementConfig& new_config) { config_ = new_config; }

  qubit_mapping_t get_placement_map(const Circuit& circ_) const override;
  // methods that return maps, for base this returns one empty map
  std::vector<qubit_mapping_t> get_all_placement_maps(
      const Circuit& circ_) const override;

 private:
  PlacementConfig config_;
};

///////////////////////////////
//   NOISE-AWARE PLACEMENT   //
///////////////////////////////

// Class for performing noise-aware placement using graph monomorphism
class Monomorpher {
 public:
  Monomorpher(
      const Circuit& _circ, const Architecture& _arc,
      const DeviceCharacterisation& _characterisation,
      const PlacementConfig& _config)
      : circ(_circ),
        arc(_arc),
        characterisation(_characterisation),
        config(_config) {
    q_graph = monomorph_interaction_graph(
        circ, config.max_interaction_edges, config.depth_limit);
  }

  // return best maps, up to max_return in number, unsorted
  std::vector<MapCost> place(unsigned max_return);
  // calculate cost of map
  double map_cost(const qubit_bimap_t& n_map);

 private:
  const Circuit& circ;
  Architecture arc;
  DeviceCharacterisation characterisation;
  PlacementConfig config;
  QubitGraph q_graph;
};

class NoiseAwarePlacement : public Placement {
 public:
  NoiseAwarePlacement(
      const Architecture& _arc,
      std::optional<avg_node_errors_t> _node_errors = std::nullopt,
      std::optional<avg_link_errors_t> _link_errors = std::nullopt,
      std::optional<avg_readout_errors_t> _readout_errors = std::nullopt) {
    arc_ = _arc;
    characterisation_ = {
        _node_errors ? *_node_errors : avg_node_errors_t(),
        _link_errors ? *_link_errors : avg_link_errors_t(),
        _readout_errors ? *_readout_errors : avg_readout_errors_t()};
    config_.depth_limit = 5;
    config_.max_interaction_edges = arc_.n_connections();
    config_.vf2_max_matches = 10000;
    config_.arc_contraction_ratio = 10;
    config_.timeout = 60000;
  }
  explicit NoiseAwarePlacement(
      const Architecture& _arc, const PlacementConfig& _config,
      std::optional<avg_node_errors_t> _node_errors = std::nullopt,
      std::optional<avg_link_errors_t> _link_errors = std::nullopt,
      std::optional<avg_readout_errors_t> _readout_errors = std::nullopt)
      : Placement(_arc), config_(_config) {
    characterisation_ = {
        _node_errors ? *_node_errors : avg_node_errors_t(),
        _link_errors ? *_link_errors : avg_link_errors_t(),
        _readout_errors ? *_readout_errors : avg_readout_errors_t()};
  }
  explicit NoiseAwarePlacement(const PlacementConfig& _config)
      : config_(_config) {}
  PlacementConfig get_config() { return config_; }
  void set_config(const PlacementConfig& new_config) { config_ = new_config; }
  qubit_mapping_t get_placement_map(const Circuit& circ_) const override;
  // methods that return maps, for base this returns one empty map
  std::vector<qubit_mapping_t> get_all_placement_maps(
      const Circuit& circ_) const override;

 private:
  friend void to_json(nlohmann::json& j, const PlacementPtr& placement_ptr);
  friend void from_json(const nlohmann::json& j, PlacementPtr& placement_ptr);

  PlacementConfig config_;
  DeviceCharacterisation characterisation_;
};

}  // namespace tket
