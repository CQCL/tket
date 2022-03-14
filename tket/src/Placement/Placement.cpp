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

#include "Placement.hpp"

#include "Utils/HelperFunctions.hpp"
#include "Utils/Json.hpp"

namespace tket {

// qubit_mapping_t -> std::map<Qubit, Node>
// qubit_map_t -> std::map<Qubit, Qubit>
// qubit_bimap_t -> boost::bimap<Qubit ArcVertex>

PlacementConfig::PlacementConfig(
    unsigned _depth_limit, unsigned _max_interaction_edges,
    unsigned _vf2_max_matches, unsigned _arc_contraction_ratio,
    unsigned _timeout)
    : depth_limit(_depth_limit),
      max_interaction_edges(_max_interaction_edges),
      vf2_max_matches(_vf2_max_matches),
      arc_contraction_ratio(_arc_contraction_ratio),
      timeout(_timeout) {}

bool PlacementConfig::operator==(const PlacementConfig &other) const {
  return (this->depth_limit == other.depth_limit) &&
         (this->max_interaction_edges == other.max_interaction_edges) &&
         (this->vf2_max_matches == other.vf2_max_matches) &&
         (this->arc_contraction_ratio == other.arc_contraction_ratio) &&
         (this->timeout == other.timeout);
}

void to_json(nlohmann::json &j, const PlacementPtr &placement_ptr) {
  j["architecture"] = placement_ptr->get_architecture_ref();
  if (std::shared_ptr<GraphPlacement> cast_placer =
          std::dynamic_pointer_cast<GraphPlacement>(placement_ptr)) {
    j["type"] = "GraphPlacement";
    j["config"] = cast_placer->get_config();
  } else if (
      std::shared_ptr<NoiseAwarePlacement> cast_placer =
          std::dynamic_pointer_cast<NoiseAwarePlacement>(placement_ptr)) {
    j["type"] = "NoiseAwarePlacement";
    j["config"] = cast_placer->get_config();
    j["characterisation"] = cast_placer->characterisation_;
  } else if (
      std::shared_ptr<LinePlacement> cast_placer =
          std::dynamic_pointer_cast<LinePlacement>(placement_ptr)) {
    j["type"] = "LinePlacement";
  } else {
    j["type"] = "Placement";
  }
}

void from_json(const nlohmann::json &j, PlacementPtr &placement_ptr) {
  std::string classname = j.at("type").get<std::string>();
  Architecture arc = j.at("architecture").get<Architecture>();
  if (classname == "GraphPlacement") {
    PlacementConfig config = j.at("config").get<PlacementConfig>();
    placement_ptr = std::make_shared<GraphPlacement>(arc, config);
  } else if (classname == "NoiseAwarePlacement") {
    PlacementConfig config = j.at("config").get<PlacementConfig>();
    DeviceCharacterisation ch =
        j.at("characterisation").get<DeviceCharacterisation>();
    std::shared_ptr<NoiseAwarePlacement> nap =
        std::make_shared<NoiseAwarePlacement>(arc, config);
    nap->characterisation_ = ch;
    placement_ptr = nap;
  } else if (classname == "LinePlacement") {
    placement_ptr = std::make_shared<LinePlacement>(arc);
  } else {
    placement_ptr = std::make_shared<Placement>(arc);
  }
}

void to_json(nlohmann::json &j, const PlacementConfig &config) {
  j["depth_limit"] = config.depth_limit;
  j["max_interaction_edges"] = config.max_interaction_edges;
  j["vf2_max_matches"] = config.vf2_max_matches;
  j["arc_contraction_ratio"] = config.arc_contraction_ratio;
  j["timeout"] = config.timeout;
}

void from_json(const nlohmann::json &j, PlacementConfig &config) {
  config.depth_limit = j.at("depth_limit").get<unsigned>();
  config.max_interaction_edges = j.at("max_interaction_edges").get<unsigned>();
  config.vf2_max_matches = j.at("vf2_max_matches").get<unsigned>();
  config.arc_contraction_ratio = j.at("arc_contraction_ratio").get<unsigned>();
  config.timeout = j.at("timeout").get<unsigned>();
}

const std::string &Placement::unplaced_reg() {
  static std::unique_ptr<const std::string> regname =
      std::make_unique<const std::string>("unplaced");
  return *regname;
}

void fill_partial_mapping(
    const qubit_vector_t &current_qubits, qubit_mapping_t &partial_mapping) {
  unsigned up_nu = 0;
  for (Qubit q : current_qubits) {
    if (partial_mapping.find(q) == partial_mapping.end()) {
      partial_mapping.insert({q, Node(Placement::unplaced_reg(), up_nu)});
      up_nu++;
    }
  }

  if (std::any_of(
          partial_mapping.begin(), partial_mapping.end(),
          [&current_qubits](auto x) {
            return std::find(
                       current_qubits.begin(), current_qubits.end(), x.first) ==
                   current_qubits.end();
          })) {
    tket_log()->warn("mapping contains qubits not present in circuit");
  }
}

// Default placement methods
bool Placement::place(
    Circuit &circ_, std::shared_ptr<unit_bimaps_t> maps) const {
  qubit_mapping_t map_ = get_placement_map(circ_);
  return place_with_map(circ_, map_, maps);
}

bool Placement::place_with_map(
    Circuit &circ_, qubit_mapping_t &map_,
    std::shared_ptr<unit_bimaps_t> maps) {
  qubit_vector_t circ_qbs = circ_.all_qubits();
  fill_partial_mapping(circ_qbs, map_);
  bool changed = circ_.rename_units(map_);
  changed |= update_maps(maps, map_, map_);
  return changed;
}

qubit_mapping_t Placement::get_placement_map(const Circuit &circ_) const {
  qubit_mapping_t out_map;
  fill_partial_mapping(circ_.all_qubits(), out_map);
  return out_map;
}

std::vector<qubit_mapping_t> Placement::get_all_placement_maps(
    const Circuit &circ_) const {
  return {get_placement_map(circ_)};
}

qubit_mapping_t NaivePlacement::get_placement_map(const Circuit &circ_) const {
  return get_all_placement_maps(circ_).at(0);
}

std::vector<qubit_mapping_t> NaivePlacement::get_all_placement_maps(
    const Circuit &circ_) const {
  qubit_mapping_t placement;
  qubit_vector_t to_place;
  std::vector<Node> placed;

  // Find which/if any qubits need placing
  for (const Qubit &q : circ_.all_qubits()) {
    Node n(q);
    if (!this->arc_.node_exists(n)) {
      to_place.push_back(n);
    } else {
      placed.push_back(n);
      // if already placed, make sure qubit retains placement
      placement.insert({n, n});
    }
  }
  // avoid doing std::set_difference unless qubits need to be placed
  unsigned n_placed = to_place.size();
  if (n_placed > 0) {
    std::vector<Node> difference,
        architecture_nodes = this->arc_.get_all_nodes_vec();
    std::set_difference(
        architecture_nodes.begin(), architecture_nodes.end(), placed.begin(),
        placed.end(), std::inserter(difference, difference.begin()));
    // should always be enough remaining qubits to assign unplaced qubits to
    TKET_ASSERT(difference.size() >= n_placed);
    for (unsigned i = 0; i < n_placed; i++) {
      // naively assign each qubit to some free node
      placement.insert({to_place[i], difference[i]});
    }
  }
  return {placement};
}

qubit_mapping_t LinePlacement::get_placement_map(const Circuit &circ_) const {
  return get_all_placement_maps(circ_).at(0);
}

std::vector<qubit_mapping_t> LinePlacement::get_all_placement_maps(
    const Circuit &circ_) const {
  qubit_mapping_t partial_map;
  QubitLineList qb_lines = qubit_lines(circ_);
  if (!qb_lines.empty()) {
    partial_map = lines_on_arc(arc_, qb_lines, circ_.n_qubits());
  }
  fill_partial_mapping(circ_.all_qubits(), partial_map);
  return {partial_map};
}

qubit_mapping_t GraphPlacement::get_placement_map(const Circuit &circ_) const {
  QubitGraph q_graph = monomorph_interaction_graph(
      circ_, arc_.n_connections(), config_.depth_limit);
  std::vector<qubit_bimap_t> all_bimaps = monomorphism_edge_break(
      arc_, q_graph, config_.vf2_max_matches, config_.timeout);
  qubit_mapping_t out_map = bimap_to_map(all_bimaps[0].left);
  fill_partial_mapping(circ_.all_qubits(), out_map);
  return out_map;
}

std::vector<qubit_mapping_t> GraphPlacement::get_all_placement_maps(
    const Circuit &circ_) const {
  QubitGraph q_graph = monomorph_interaction_graph(
      circ_, arc_.n_connections(), config_.depth_limit);
  std::vector<qubit_bimap_t> all_bimaps = monomorphism_edge_break(
      arc_, q_graph, config_.vf2_max_matches, config_.timeout);
  std::vector<qubit_mapping_t> all_qmaps;
  qubit_vector_t all_qbs = circ_.all_qubits();
  for (qubit_bimap_t bm : all_bimaps) {
    qubit_mapping_t qm = bimap_to_map(bm.left);
    fill_partial_mapping(all_qbs, qm);
    all_qmaps.push_back(qm);
  }
  return all_qmaps;
}

qubit_mapping_t NoiseAwarePlacement::get_placement_map(
    const Circuit &circ_) const {
  return get_all_placement_maps(circ_)[0];
}

std::vector<qubit_mapping_t> NoiseAwarePlacement::get_all_placement_maps(
    const Circuit &circ_) const {
  // Set up default placements configuation
  // Place according to configuration
  Monomorpher placer(circ_, arc_, characterisation_, config_);
  std::vector<MapCost> results = placer.place(config_.depth_limit * 2);
  std::sort(results.begin(), results.end());
  std::vector<qubit_mapping_t> output;
  qubit_vector_t all_qbs = circ_.all_qubits();
  for (const MapCost &map_c : results) {
    qubit_mapping_t map_ = map_c.map;
    fill_partial_mapping(all_qbs, map_);
    output.push_back(map_c.map);
  }
  return output;
}

}  // namespace tket
