#ifndef _TKET_BoxDecomposition_H_
#define _TKET_BoxDecomposition_H_

#include "Mapping/MappingFrontier.hpp"
#include "Mapping/RoutingMethod.hpp"

namespace tket {

class BoxDecomposition {
 public:
  /**
   * Class Constructor
   * @param _architecture Architecture object added operations must respect
   * @param _mapping_frontier Contains Circuit object to be modified
   */
  BoxDecomposition(
      const ArchitecturePtr& _architecture,
      std::shared_ptr<MappingFrontier>& _mapping_frontier);

  /**
   * Decompose any boxes in the next slice after the frontier
   */
  void solve();

 private:
  // Architecture all new physical operations must respect
  ArchitecturePtr architecture_;
  std::shared_ptr<MappingFrontier> mapping_frontier_;
};

class BoxDecompositionRoutingMethod : public RoutingMethod {
 public:
  /**
   * Decompose any boxes on the frontier
   */
  BoxDecompositionRoutingMethod();

  /**
   * @return true if method can route subcircuit, false if not
   */
  bool check_method(
      const std::shared_ptr<MappingFrontier>& mapping_frontier,
      const ArchitecturePtr& /*architecture*/) const override;

  /**
   * @param mapping_frontier Contains boundary of routed/unrouted circuit for
   * modifying
   * @param architecture Architecture providing physical constraints
   * @return Logical to Physical mapping at boundary due to modification.
   *
   */
  unit_map_t routing_method(
      std::shared_ptr<MappingFrontier>& mapping_frontier,
      const ArchitecturePtr& architecture) const override;

  nlohmann::json serialize() const override;

  static BoxDecompositionRoutingMethod deserialize(const nlohmann::json& /*j*/);
};

}  // namespace tket

#endif