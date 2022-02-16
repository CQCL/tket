#ifndef _TKET_AASRoute_H_
#define _TKET_AASRoute_H_

#include "ArchAwareSynth/SteinerForest.hpp"
#include "Mapping/LexicographicalComparison.hpp"
#include "Mapping/MappingFrontier.hpp"
#include "Mapping/RoutingMethod.hpp"

namespace tket {

class AASRouteError : public std::logic_error {
 public:
  explicit AASRouteError(const std::string& message)
      : std::logic_error(message) {}
};

// Child class of RoutingMethod, with overloaded methods for routing
// MappingFrontier objects
class AASRouteRoutingMethod : public RoutingMethod {
 public:
  /**
   * Checking and Routing methods redefined using LexiRoute. Only circuit depth,
   * corresponding to lookahead, is a required parameter.
   *
   * @param aaslookahead lookahead that should be used in the aas routing
   * @param cnotsynthtype type of cnot synthesis that should be used
   */
  AASRouteRoutingMethod(
      unsigned aaslookahead,
      aas::CNotSynthType cnotsynthtype = aas::CNotSynthType::Rec);

  /**
   * @return true if this method can route subcircuit, false if not
   * @param mapping_frontier mapping frontier contaning the circuit that should
   * be checked for routing
   * @param architecture architecture that is used for the routing check
   */
  bool check_method(
      const std::shared_ptr<MappingFrontier>& mapping_frontier,
      const ArchitecturePtr& architecture) const;

  /**
   * @param mapping_frontier Contains boundary of routed/unrouted circuit for
   * modifying
   * @param architecture Architecture providing physical constraints
   * @return Logical to Physical mapping at boundary due to modification.
   *
   */
  unit_map_t routing_method(
      std::shared_ptr<MappingFrontier>& mapping_frontier,
      const ArchitecturePtr& architecture) const;

  /**
   * @return cnot synth type of this routing method
   */
  aas::CNotSynthType get_cnotsynthtype() const;

  /**
   * @return aaslookahead of this routing method
   */
  unsigned get_aaslookahead() const;

  nlohmann::json serialize() const override;

  static AASRouteRoutingMethod deserialize(const nlohmann::json& j);

 private:
  // type of cnot synthesis that should be used
  aas::CNotSynthType cnotsynthtype_;
  // lookahead that should be used in the aas routing
  unsigned aaslookahead_;
};
}  // namespace tket

#endif