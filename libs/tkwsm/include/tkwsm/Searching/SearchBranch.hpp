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

#pragma once
#include <memory>

#include "tkwsm/GraphTheoretic/DomainInitialiser.hpp"
#include "tkwsm/Reducing/DerivedGraphsReducer.hpp"
#include "tkwsm/Reducing/DistancesReducer.hpp"
#include "tkwsm/Reducing/HallSetReduction.hpp"
#include "tkwsm/Reducing/ReducerWrapper.hpp"
#include "tkwsm/Searching/DomainsAccessor.hpp"
#include "tkwsm/Searching/NodeListTraversal.hpp"
#include "tkwsm/Searching/NodesRawData.hpp"
#include "tkwsm/Searching/WeightCalculator.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct ExtraStatistics;
class NearNeighboursData;
class NeighboursData;
class WeightChecker;

/** Represents a depth-first traversal.
 * However, it's really about reducing nodes etc.,
 * not the search logic. This contains the list of nodes so far,
 * which may be viewed as a stack, but the caller is responsible for
 * moving up and down the list.
 */
class SearchBranch {
 public:
  /** Constructs node[0], i.e. the initial node.
   * The search branch is also entirely responsible for ExtraStatistics.
   */
  SearchBranch(
      const DomainInitialiser::InitialDomains& initial_domains,
      const NeighboursData& pattern_ndata,
      NearNeighboursData& pattern_near_ndata,
      const NeighboursData& target_ndata, NearNeighboursData& target_near_ndata,
      unsigned max_distance_reduction_value, ExtraStatistics& extra_statistics);

  /** Extra parameters to configure reducing a single node. */
  struct ReductionParameters {
    WeightWSM max_weight;
  };

  /** If it returns false, the node is in an invalid state; but that's OK
   * because the caller should immediately discard it and move up.
   */
  bool reduce_current_node(const ReductionParameters& parameters);

  /** Keep decreasing the node index and reducing, until EITHER it's reduced
   * (returning true), OR we can't move up any more (and return false).
   */
  bool backtrack(const ReductionParameters& parameters);

  /** ASSUMING that the current node has been reduced,
   * make the given assignment [which must be valid],
   * and move down to a new node.
   * @param p_vertex A pattern vertex, pv.
   * @param t_vertex A target vertex, tv. We are making the new assignment
   * pv->tv.
   */
  void move_down(VertexWSM p_vertex, VertexWSM t_vertex);

  /** A relatively expensive operation. Simply returns ALL target vertices
   * which still exist in SOME domain, ANYWHERE in the search tree.
   * Thus, anything NOT in this set can be safely erased
   * (as far as THIS search is concerned, anyway).
   * Naturally, as the search progresses this set will reduce,
   * so for best results the caller should delay this as long as possible.
   * @return The set of all target vertices which might still occur during the
   * future search.
   */
  boost::dynamic_bitset<> get_used_target_vertices() const;

  void activate_weight_checker(WeightWSM total_p_edge_weights);

  const DomainsAccessor& get_domains_accessor() const;

  /** Should only rarely be used. */
  DomainsAccessor& get_domains_accessor_nonconst();

  std::size_t get_number_of_possible_assignments_from_initial_domains() const;

  std::size_t get_total_number_of_assignments_tried() const;

  /** Updates the statistics before returning them. */
  const ExtraStatistics& get_updated_extra_statistics();

 private:
  const NeighboursData& m_pattern_ndata;
  const NeighboursData& m_target_ndata;
  ExtraStatistics& m_extra_statistics;
  DerivedGraphsReducer m_derived_graphs_reducer;
  const WeightCalculator m_weight_calculator;
  HallSetReduction m_hall_set_reduction;

  NodesRawDataWrapper m_nodes_raw_data_wrapper;
  DomainsAccessor m_domains_accessor;
  NodeListTraversal m_node_list_traversal;

  // A bit crude, but to fit into the general framework we want each distance
  // to have a separate reducer object.
  // Since they all share references to the main data,
  // this is not really inefficient.
  std::vector<DistancesReducer> m_distance_reducers;
  std::vector<ReducerWrapper> m_reducer_wrappers;

  boost::dynamic_bitset<> m_work_set_for_reducers;

  /** Every pv->tv assignment in here has previously passed all
   * the simple checks for validity, so we need not repeat.
   * (Also, we don't need to store the impossible ones,
   * since the caller should erase them from all data as soon as they fail
   * the checks).
   */
  PossibleAssignments m_checked_assignments;

  std::unique_ptr<WeightChecker> m_weight_checker_ptr;

  // It's very tempting to use boost::dynamic_bitset<> here;
  // however, impossible_target_vertices are actually very rare
  // so it's not so important.
  std::set<VertexWSM> m_impossible_target_vertices;

  bool perform_single_assignment_checks_in_reduce_loop(
      std::size_t num_assignments_alldiff_processed);

  bool check_and_update_scalar_product_in_reduce_loop(
      const ReductionParameters& parameters,
      std::size_t num_assignments_alldiff_processed);

  bool perform_weight_nogood_check_in_reduce_loop(
      const ReductionParameters& parameters);

  ReductionResult perform_reducers_in_reduce_loop();

  bool perform_main_reduce_loop(const ReductionParameters& parameters);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
