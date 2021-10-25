#include "main_entry_functions.hpp"

#include <sstream>
#include <stdexcept>

#include "BestFullTsa.hpp"
#include "TSAUtils/VertexMappingFunctions.hpp"
#include "Utils/Assert.hpp"

namespace tket {

using namespace tsa_internal;

std::vector<std::pair<Node, Node>> get_swaps(
    const Architecture& architecture, const NodeMapping& node_mapping) {
  std::vector<std::pair<Node, Node>> swaps;
  // Before all the conversion and object construction,
  // doesn't take long to check if it's actually trivial
  bool trivial = true;
  for (const auto& entry : node_mapping) {
    if (entry.first != entry.second) {
      trivial = false;
      break;
    }
  }
  if (trivial) {
    return swaps;
  }
  // Now convert the Nodes into raw vertices for use in TSA objects.
  const ArchitectureMapping arch_mapping(architecture);
  VertexMapping vertex_mapping;
  for (const auto& node_entry : node_mapping) {
    vertex_mapping[arch_mapping.get_vertex(node_entry.first)] =
        arch_mapping.get_vertex(node_entry.second);
  }
  TKET_ASSERT(vertex_mapping.size() == node_mapping.size());
  check_mapping(vertex_mapping);

  SwapList raw_swap_list;
  BestFullTsa().append_partial_solution(
      raw_swap_list, vertex_mapping, arch_mapping);

  // Finally, convert the raw swaps back to nodes.
  swaps.reserve(raw_swap_list.size());
  for (auto id_opt = raw_swap_list.front_id(); id_opt;
       id_opt = raw_swap_list.next(id_opt.value())) {
    const auto& raw_swap = raw_swap_list.at(id_opt.value());
    swaps.emplace_back(std::make_pair(
        arch_mapping.get_node(raw_swap.first),
        arch_mapping.get_node(raw_swap.second)));
  }
  return swaps;
}

std::tuple<Circuit, unit_map_t, unit_map_t> get_swaps(
    const Architecture& architecture,
    const unit_map_t& initial_logical_to_physical_map,
    const unit_map_t& desired_logical_to_physical_map) {
  // The physical qubits are nodes inside the architecture.
  // Some Node <--> UnitID conversion is unavoidable with the current design,
  // since Architecture uses Node objects, rather than UnitID objects,
  // and types like vector<Node> and vector<UnitID> cannot be converted
  // to each other without copying, even though each Node is just
  // a UnitID with no extra data (C++ containers are not "covariant").
  NodeMapping node_mapping;
  for (const std::pair<const UnitID, UnitID>& initial_entry :
       initial_logical_to_physical_map) {
    const auto citer =
        desired_logical_to_physical_map.find(initial_entry.first);
    if (citer == desired_logical_to_physical_map.cend()) {
      std::stringstream ss;
      ss << "Logical qubit " << initial_entry.first.repr()
         << " is present in the initial logical->physical map, but not in the "
            "target logical->physical map";
      throw std::runtime_error(ss.str());
    }
    const Node source_physical_node(initial_entry.second);
    const Node target_physical_node(citer->second);
    node_mapping[source_physical_node] = target_physical_node;
  }
  if (initial_logical_to_physical_map.size() !=
      desired_logical_to_physical_map.size()) {
    std::stringstream ss;
    ss << "Initial and final logical->physical mappings have different sizes "
       << initial_logical_to_physical_map.size() << ", "
       << desired_logical_to_physical_map.size()
       << ". There are extra logical qubits in the final map missing from the "
          "initial map";
    throw std::runtime_error(ss.str());
  }
  if (node_mapping.size() != initial_logical_to_physical_map.size()) {
    std::stringstream ss;
    ss << "Converted " << initial_logical_to_physical_map.size()
       << " distinct logical qubits to " << node_mapping.size()
       << " distinct physical nodes; conversion error";
    throw std::runtime_error(ss.str());
  }
  const auto node_swaps = get_swaps(architecture, node_mapping);

  // Don't add unused nodes to the final circuit.
  std::set<Node> nodes_seen;
  for (const auto& swap : node_swaps) {
    nodes_seen.insert(swap.first);
    nodes_seen.insert(swap.second);
  }

  std::tuple<Circuit, unit_map_t, unit_map_t> result;

  // We rely on the algorithm to be correct,
  // i.e. it really has calculated the full desired mapping.
  //
  // NOTE: other nodes in the architecture might be involved in the swaps,
  // even if they were not mentioned in any of the input logical->physical maps.
  // But that's OK; if the caller wants to keep them fixed,
  // they should have put them into the input maps.
  std::get<1>(result) = initial_logical_to_physical_map;
  std::get<2>(result) = desired_logical_to_physical_map;

  for (const Node& node : nodes_seen) {
    std::get<0>(result).add_qubit(node);
  }
  // Now we can add the swaps.
  for (const auto& swap : node_swaps) {
    std::get<0>(result).add_op<Node>(OpType::SWAP, {swap.first, swap.second});
  }
  return result;
}

}  // namespace tket
