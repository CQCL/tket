#ifndef _TKET_LexiRoute_H_
#define _TKET_LexiRoute_H_

#include "Mapping/LexicographicalComparison.hpp"
#include "Mapping/MappingFrontier.hpp"
#include "Mapping/RoutingMethod.hpp"

namespace tket {

class LexiRouteError : public std::logic_error {
 public:
  explicit LexiRouteError(const std::string& message)
      : std::logic_error(message) {}
};

/**
 * A class for modifiying a Circuit held in a MappingFrontier object
 * with either an Architecture permitted single SWAP gate or BRIDGE gate.
 * Used in the LexiRouteRoutingMethod class which provides a subcircuit
 * modification method for MappingManager. Used in solution presented in "On the
 * qubit routing problem" -> arXiv:1902.08091
 */
class LexiRoute {
 public:
  /**
   * Class Constructor
   * @param _architecture Architecture object added operations must respect
   * @param _mapping_frontier Contains Circuit object to be modified
   */
  LexiRoute(
      const ArchitecturePtr& _architecture,
      std::shared_ptr<MappingFrontier>& _mapping_frontier);

  /**
   * When called, LexiRoute::solve will modify the Circuit held in
   * MappingFrontier object passed at class construction. Either a SWAP gate
   * will be inserted at the input boundary of the held Circuit or a CX gate
   * will be transformed into a BRIDGE gate. The added SWAP or BRIDGE gate will
   * be valid for the Architecture passed at class construction. Additionally,
   * an "unlabelled" Qubit in the Circuit may be relabelled to a Node in the
   * Architecture, or an "unlabelled" Qubit may have its path merged with an
   * ancilla qubit.
   * The decision making is based on the heuristic outlined in arXiv:1902.08091.
   *
   * @param lookahead Number of slices to lookahead at when determining best
   * SWAP or BRIDGE
   */
  void solve(unsigned lookahead);

 private:
  /**
   * this->interacting_uids_ attribute is a map where key is one UnitID
   * and value is the UnitID it needs to be adjacent to.
   * This map is implicitly updated whenever a logical SWAP is inserted.
   * set_interacting_uids determines this map for the first parallel set of
   * interacting UnitID in the Circuit held in this->mapping_frontier_
   * @param assigned_only If true, only include interactions where both UnitID
   * are in this->architecture_.
   */
  void set_interacting_uids(bool assigned_only = false);

  /**
   * Merges the qubit paths of "merge" and "ancilla" in mapping frontier circuit
   * such that the output of the final ancilla vertex leads into the input of
   * the first merge vertex.
   *
   * @param merge UnitID to which ancilla path is prepended
   * @param ancilla UnitID of ancilla opeartions
   */
  void merge_with_ancilla(const UnitID& merge, const UnitID& ancilla);

  /**
   * If there is some "free" Node in Architecture at distance "distances" on
   * the connectivity graph, assign (relable) UnitID assignee to it. "free"
   * => not in Circuit. If no unassigned node at distances from root, return
   * false.
   * @param assignee UnitID not in Architecture to relabel
   * @param root Node in Architecture
   * @param distances Distance at which to find free Node from root at
   * @return True if assigned, else False
   */
  bool assign_at_distance(
      const UnitID& assignee, const Node& root, unsigned distances);

  /**
   * If this->set_interacting_uids assigned_only bool is false then the
   * this->interacting_uids attribute may have key and value UnitID not in
   * this->architecture_.
   * update_labelling assigns these non-architecture UnitID to some Architecture
   * UnitID, updating the this->labelling_ attribute.
   * @return True if anything relabelled, else false
   */
  bool update_labelling();

  /**
   * Returns a set of pair of UnitID, each denoting a SWAP.
   * Returned SWAP have at least one UnitID in interacting_uids_.
   * This is such that enacting any of these SWAP will alter the distance
   * between some interacting UnitID.
   * @return std::pair<UnitID, UnitID> suitable for addition to Circuit
   */
  swap_set_t get_candidate_swaps();

  /**
   * Proposed swap will have two Node
   * Each of these Node may be in some interaction in the first layer of circuit
   * held in mapping_frontier. If either of these Node are in an interaction,
   * check whether said interaction is a CX interaction, and if the pair of Node
   * in the interaction are at distance 2. If true, compare lexicographical
   * distances between no swap and given swap assuming distance 2 interactions
   * are  complete. If no swap is better, update return object to reflect this.
   * @param swap Pair of Node comprising SWAP for checking
   * @param lookahead Number of steps of lookahead emplyed for comparison
   * @return Pair of bool, where true implies BRIDGE to be added
   */
  std::pair<bool, bool> check_bridge(
      const std::pair<Node, Node>& swap, unsigned lookahead);

  /**
   * Returns a pair of distances, where the distances are between n1 & p1, and
   * n2 & p2. Pair object is ordered such that the greatest distance is first.
   *
   * @param p0_first First Node in first interaction to find distance between
   * @param p0_second Second Node in first interaction to find distance between
   * @param p1_first First Node in second interaction to find distance between
   * @param p1_second Second Node in second interaction to find distance between
   * @return Pair of size_t, being distances on architecture graph
   */
  const std::pair<size_t, size_t> pair_distances(
      const Node& p0_first, const Node& p0_second, const Node& p1_first,
      const Node& p1_second) const;

  /**
   * It is always expected that at least one Node in a SWAP will be in some
   * interaction. This method checks that the given swap will strictly decrease
   * the distance for this interaction, and removes it from the swaps set if
   * not.
   *
   * @param swaps Potential swaps to remove from
   */
  void remove_swaps_decreasing(swap_set_t& swaps);

  // Architecture all new physical operations must respect
  ArchitecturePtr architecture_;
  //   Contains circuit for finding SWAP from and non-routed/routed boundary
  std::shared_ptr<MappingFrontier>& mapping_frontier_;
  //   Map between UnitID and UnitID they interact with at boundary
  unit_map_t interacting_uids_;
  //   Map between original circuit UnitID and new UnitID due to dynamic
  //   placement
  unit_map_t labelling_;
  //   Set tracking which Architecture Node are present in Circuit
  std::set<Node> assigned_nodes_;
};

// Child class of RoutingMethod, with overloaded methods for routing
// MappingFrontier objects
class LexiRouteRoutingMethod : public RoutingMethod {
 public:
  /**
   * Checking and Routing methods redefined using LexiRoute. Only circuit depth,
   * corresponding to lookahead, is a required parameter.
   *
   * @param _max_depth Number of layers of gates checked inr outed subcircuit.
   */
  LexiRouteRoutingMethod(unsigned _max_depth);

  /**
   * @return true if method can route subcircuit, false if not
   */
  bool check_method(
      const std::shared_ptr<MappingFrontier>& /*mapping_frontier*/,
      const ArchitecturePtr& /*architecture*/) const;

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

 private:
  unsigned max_depth_;
};


JSON_DECL(LexiRouteRoutingMethod)

}  // namespace tket

#endif