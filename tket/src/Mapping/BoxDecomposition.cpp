#include "Mapping/BoxDecomposition.hpp"

#include "Mapping/MappingFrontier.hpp"

namespace tket {

BoxDecomposition::BoxDecomposition(
    const ArchitecturePtr &_architecture,
    std::shared_ptr<MappingFrontier> &_mapping_frontier)
    : architecture_(_architecture), mapping_frontier_(_mapping_frontier) {}

void BoxDecomposition::solve() {
  // Box type vertices are later removed from DAG
  VertexList bin;

  std::shared_ptr<unit_frontier_t> frontier_edges =
      frontier_convert_vertport_to_edge(
          this->mapping_frontier_->circuit_,
          this->mapping_frontier_->quantum_boundary);
  CutFrontier next_cut =
      this->mapping_frontier_->circuit_.next_q_cut(frontier_edges);
  for (Vertex &vert : *next_cut.slice) {
    Op_ptr op = this->mapping_frontier_->circuit_.get_Op_ptr_from_Vertex(vert);
    bool conditional = op->get_type() == OpType::Conditional;
    if (conditional) {
      const Conditional &cond = static_cast<const Conditional &>(*op);
      op = cond.get_op();
    }
    if (!op->get_desc().is_box()) continue;
    const Box &b = static_cast<const Box &>(*op);
    Circuit replacement = *b.to_circuit();
    if (conditional) {
      this->mapping_frontier_->circuit_.substitute_conditional(
          replacement, vert, Circuit::VertexDeletion::No,
          Circuit::OpGroupTransfer::Merge);
    } else {
      this->mapping_frontier_->circuit_.substitute(
          replacement, vert, Circuit::VertexDeletion::No,
          Circuit::OpGroupTransfer::Merge);
    }
    bin.push_back(vert);
  }

  // Delete vertices
  this->mapping_frontier_->circuit_.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
}

BoxDecompositionRoutingMethod::BoxDecompositionRoutingMethod(){};

bool BoxDecompositionRoutingMethod::check_method(
    const std::shared_ptr<MappingFrontier> &mapping_frontier,
    const ArchitecturePtr &architecture) const {
  // Not currently in use
  (void)architecture;
  std::shared_ptr<unit_frontier_t> frontier_edges =
      frontier_convert_vertport_to_edge(
          mapping_frontier->circuit_, mapping_frontier->quantum_boundary);
  CutFrontier next_cut = mapping_frontier->circuit_.next_q_cut(frontier_edges);
  for (const Vertex &vert : *next_cut.slice) {
    Op_ptr op = mapping_frontier->circuit_.get_Op_ptr_from_Vertex(vert);
    if (op->get_desc().is_box() ||
        (op->get_type() == OpType::Conditional &&
         static_cast<const Conditional &>(*op).get_op()->get_desc().is_box()))
      return true;
  }
  return false;
}

unit_map_t BoxDecompositionRoutingMethod::routing_method(
    std::shared_ptr<MappingFrontier> &mapping_frontier,
    const ArchitecturePtr &architecture) const {
  BoxDecomposition bd(architecture, mapping_frontier);
  bd.solve();
  return {};
}

}  // namespace tket