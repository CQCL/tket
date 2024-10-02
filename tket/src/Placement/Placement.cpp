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

#include "tket/Placement/Placement.hpp"

namespace tket {

const std::string& Placement::unplaced_reg() {
  static std::unique_ptr<const std::string> regname =
      std::make_unique<const std::string>("unplaced");
  return *regname;
}

void fill_partial_mapping(
    const qubit_vector_t& current_qubits,
    std::map<Qubit, Node>& partial_mapping) {
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
    tket_log()->warn(
        "Placement map has some Qubit not present in the Circuit.");
  }
}

Placement::Placement(const Architecture& _architecture)
    : architecture_(_architecture) {};

bool Placement::place(
    Circuit& circ_, std::shared_ptr<unit_bimaps_t> compilation_map) const {
  if (circ_.n_qubits() > this->architecture_.n_nodes()) {
    throw std::invalid_argument(
        "Circuit has more qubits than Architecture has nodes.");
  }
  std::map<Qubit, Node> map_ = this->get_placement_map(circ_);
  return this->place_with_map(circ_, map_, compilation_map);
}

bool Placement::place_with_map(
    Circuit& circ_, std::map<Qubit, Node>& map_,
    std::shared_ptr<unit_bimaps_t> compilation_map) {
  bool changed = circ_.rename_units(map_);
  changed |= update_maps(compilation_map, map_, map_);
  return changed;
}

std::map<Qubit, Node> Placement::get_placement_map(const Circuit& circ_) const {
  std::vector<std::map<Qubit, Node>> all_maps =
      this->get_all_placement_maps(circ_, 1);

  // basic handling to avoid segmentation faults, as placement method may not
  // return any valid map
  auto it = all_maps.begin();
  if (it != all_maps.end()) {
    std::map<Qubit, Node> map = *it;
    fill_partial_mapping(circ_.all_qubits(), map);
    return map;
  } else {
    return {};
  }
}

std::vector<std::map<Qubit, Node>> Placement::get_all_placement_maps(
    const Circuit& circ_, unsigned /*matches*/) const {
  std::map<Qubit, Node> placement;
  qubit_vector_t to_place;
  std::vector<Node> placed;

  // Find which/if any qubits need placing
  for (const Qubit& q : circ_.all_qubits()) {
    Node n(q);
    if (!this->architecture_.node_exists(n)) {
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
        architecture_nodes = this->architecture_.get_all_nodes_vec();
    std::set_difference(
        architecture_nodes.begin(), architecture_nodes.end(), placed.begin(),
        placed.end(), std::inserter(difference, difference.begin()));
    // should always be enough remaining qubits to assign unplaced qubits to
    if (difference.size() < n_placed) {
      throw std::invalid_argument(
          "There are more unplaced Qubit in Circuit then there are free Nodes "
          "in Architecture.");
    }
    for (unsigned i = 0; i < n_placed; i++) {
      // naively assign each qubit to some free node
      placement.insert({to_place[i], difference[i]});
    }
  }
  return {placement};
}

void to_json(nlohmann::json& j, const Placement::Ptr& placement_ptr) {
  // n.b. due to inheritance NoiseAwarePlacement and LinePlacement
  // need to be cast before GraphPlacement
  j["architecture"] = placement_ptr->get_architecture_ref();
  if (std::shared_ptr<LinePlacement> cast_placer =
          std::dynamic_pointer_cast<LinePlacement>(placement_ptr)) {
    j["type"] = "LinePlacement";
    j["maximum_pattern_gates"] = cast_placer->get_maximum_pattern_gates();
    j["maximum_pattern_depth"] = cast_placer->get_maximum_pattern_depth();
  } else if (
      std::shared_ptr<NoiseAwarePlacement> cast_placer =
          std::dynamic_pointer_cast<NoiseAwarePlacement>(placement_ptr)) {
    j["type"] = "NoiseAwarePlacement";
    j["matches"] = cast_placer->get_maximum_matches();
    j["timeout"] = cast_placer->get_timeout();
    j["maximum_pattern_gates"] = cast_placer->get_maximum_pattern_gates();
    j["maximum_pattern_depth"] = cast_placer->get_maximum_pattern_depth();
    j["characterisation"] = cast_placer->get_characterisation();
  } else if (
      std::shared_ptr<GraphPlacement> cast_placer =
          std::dynamic_pointer_cast<GraphPlacement>(placement_ptr)) {
    j["type"] = "GraphPlacement";
    j["matches"] = cast_placer->get_maximum_matches();
    j["timeout"] = cast_placer->get_timeout();
    j["maximum_pattern_gates"] = cast_placer->get_maximum_pattern_gates();
    j["maximum_pattern_depth"] = cast_placer->get_maximum_pattern_depth();
  } else {
    j["type"] = "Placement";
  }
}

void from_json(const nlohmann::json& j, Placement::Ptr& placement_ptr) {
  std::string classname = j.at("type").get<std::string>();
  Architecture arc = j.at("architecture").get<Architecture>();
  if (classname == "GraphPlacement") {
    unsigned matches = j.at("matches").get<unsigned>();
    unsigned timeout = j.at("timeout").get<unsigned>();
    unsigned max_pattern_gates = j.at("maximum_pattern_gates").get<unsigned>();
    unsigned max_pattern_depth = j.at("maximum_pattern_depth").get<unsigned>();
    placement_ptr = std::make_shared<GraphPlacement>(
        arc, matches, timeout, max_pattern_gates, max_pattern_depth);
  } else if (classname == "LinePlacement") {
    unsigned max_pattern_gates = j.at("maximum_pattern_gates").get<unsigned>();
    unsigned max_pattern_depth = j.at("maximum_pattern_depth").get<unsigned>();
    std::shared_ptr<LinePlacement> lp = std::make_shared<LinePlacement>(
        arc, max_pattern_gates, max_pattern_depth);
    placement_ptr = lp;
  } else if (classname == "NoiseAwarePlacement") {
    unsigned matches = j.at("matches").get<unsigned>();
    unsigned timeout = j.at("timeout").get<unsigned>();
    unsigned max_pattern_gates = j.at("maximum_pattern_gates").get<unsigned>();
    unsigned max_pattern_depth = j.at("maximum_pattern_depth").get<unsigned>();
    DeviceCharacterisation characterisation =
        j.at("characterisation").get<DeviceCharacterisation>();
    avg_node_errors_t empty_node_errors = {};
    avg_readout_errors_t empty_readout_errors = {};
    avg_link_errors_t empty_link_errors = {};
    std::shared_ptr<NoiseAwarePlacement> nap =
        std::make_shared<NoiseAwarePlacement>(
            arc, empty_node_errors, empty_link_errors, empty_readout_errors,
            matches, timeout, max_pattern_gates, max_pattern_depth);
    nap->set_characterisation(characterisation);
    placement_ptr = nap;
  } else {
    placement_ptr = std::make_shared<Placement>(arc);
  }
}
}  // namespace tket