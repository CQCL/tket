#include "Mapping/LexiLabelling.hpp"

namespace tket {

bool LabellingRoutingMethod::check_method(
    const std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  std::set<Node> already_checked;
  // *it = {UnitID, {Vertex, Port}}
  for (auto it = mapping_frontier->quantum_boundary->get<TagKey>().begin();
       it != mapping_frontier->quantum_boundary->get<TagKey>().end(); ++it) {
    Edge e0 = mapping_frontier->circuit_.get_nth_out_edge(
        it->second.first, it->second.second);
    Vertex v0 = mapping_frontier->circuit_.target(e0);
    Node node = Node(it->first);
    // i.e. skip already checked vertices
    if (already_checked.find(node) == already_checked.end()) {
      already_checked.insert(node);
      // for coutner<n_edges check, to avoid iterating through many many qubits
      int n_edges =
          mapping_frontier->circuit_.n_in_edges_of_type(v0, EdgeType::Quantum);
      int counter = 1;  // 1 edge
      auto jt = it;
      ++jt;
      for (; jt != mapping_frontier->quantum_boundary->get<TagKey>().end() &&
             counter < n_edges;
           ++jt) {
        Edge e1 = mapping_frontier->circuit_.get_nth_out_edge(
            jt->second.first, jt->second.second);
        Vertex v1 = mapping_frontier->circuit_.target(e1);
        if (v0 == v1) {
          counter++;
          // confirms that there is at least one multi-qubit gate in the first
          // layer which is assigned to some Qubit not in the architecture
          if (!architecture->node_exists(Node(it->first)) ||
              !architecture->node_exists(Node(jt->first))) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

unit_map_t LabellingRoutingMethod::routing_method(
    std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  LexiRoute lr(architecture, mapping_frontier);
  lr.solve_labelling();
  return {};
}

nlohmann::json LabellingRoutingMethod::serialize() const {
  nlohmann::json j;
  j["name"] = "LabellingRoutingMethod";
  return j;
}

LabellingRoutingMethod LabellingRoutingMethod::deserialize(
    const nlohmann::json& /*j*/) {
  return LabellingRoutingMethod();
}

}  // namespace tket
