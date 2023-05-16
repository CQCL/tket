// Copyright 2019-2023 Cambridge Quantum Computing
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
#include <boost/dynamic_bitset.hpp>

#include "tkwsm/Common/ReusableStorage.hpp"
#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class DomainsAccessor;

/** This is a kind of "Reducer" object, but not derived from
 * ReducerInterface because it doesn't care about individual assignments.
 *
 * This tries to find sets of pattern vertices {v(i)} such that
 *    |U| = |union_i Domain(v(i))| = |{v(i)}| > 1.
 * (This is termed a "Hall set", as in the 2015 Glasgow subgraph solver paper
 * "A Parallel, Backjumping Subgraph Isomorphism Algorithm using Supplemental
 * Graphs" by Ciaran McCreesh and Patrick Prosser. See also
 * Hall's marriage theorem).
 *
 * When this occurs, reduction is possible (since, Domain(v) must be
 * disjoint from U for all other v).
 * Moreover, it is a nogood if |U| < |{v(i)}|.
 *
 * We do not guarantee to detect Hall sets in all cases,
 * which would probably need exponential time; but it detects many
 * examples quite quickly.
 *
 * Note also that this is the only kind of reducer which can be useful
 * for unassigned PV.
 * All other reducers take action only when an assignment PV->TV is made.
 */
class HallSetReduction {
 public:
  HallSetReduction();

  /** Call at the start of reducing a new search node. */
  void clear();

  /** Reduce domains, if Hall sets are found. This can be called repeatedly
   * when reducing the same node.
   * @param accessor An object which provides read and write access to the
   * current node.
   * @return Information about the reduction.
   */
  ReductionResult reduce(DomainsAccessor& accessor);

 private:
  struct Data {
    VertexWSM pv;
    std::size_t domain_size;
  };

  // A partition is an object containing data about a set of pattern vertices
  // and their domains.
  // Each partition should have vertices and domains
  // disjoint from those of every other partition.
  //
  // (Thus, we can search for Hall sets, and reduce domains,
  // within each partition, independently of all others; they cannot
  // affect each other, as far as Hall set reduction is concerned).
  //
  // Note that, even if a Hall set is detected, it still may be worth
  // storing the domains and checking later for further reductions.
  // E.g., if initially we have Dom(x1) = Dom(x2) = Dom(x3) = {y1, y2, y3},
  // then y1,y2,y3 can be removed from all other domains.
  // But then we move x1,x2,x3 into their own new partition.
  // If in future some other reduction,
  // unrelated to Hall sets, reduces it to Dom(x1) = Dom(x2) = {y1,y2},
  // then Hall set reduction will detect this and force x3->y3 (or a nogood,
  // if y3 was removed from Dom(x3) by some other reduction).
  struct Partition {
    std::vector<Data> domains_data;

    // Use this to track consecutive occasions when
    // we are asked to recheck for a Hall set,
    // but the domain sizes are unchanged.
    unsigned attempt_count;

    enum class FillResult { UNCHANGED, CHANGED, ERASE_THIS_PARTITION, NOGOOD };

    FillResult initial_fill(const DomainsAccessor& accessor);

    // Does NOT update attempt_count; but checks existing
    // data to see if it changed.
    FillResult refill(const DomainsAccessor& accessor);

    enum class SearchResult { NOGOOD, NO_HALL_SET, HALL_SET };

    // Assume already sorted, and only called with nonempty data.
    SearchResult search_for_hall_set(
        const DomainsAccessor& accessor,
        boost::dynamic_bitset<>& union_of_domains_work_set) const;

    struct HallSetReductionData {
      ReductionResult result_to_return;
      bool changed;
      bool should_move_hall_pv_to_another_partition;
      bool should_erase_this_partition;
    };

    // If we know the Hall set size, we can get the PV;
    // they start at the back, and we know their number.
    // Reduce the domains of the OTHER PV in this partition (since, they are
    // by construction the ONLY PV which have domains intersecting the
    // given union_of_domains).
    HallSetReductionData reduce_with_hall_set(
        DomainsAccessor& accessor,
        const boost::dynamic_bitset<>& union_of_domains,
        const boost::dynamic_bitset<>& union_of_domains_complement,
        boost::dynamic_bitset<>& domain_mask_workset);

    // Moves the smallest domain sizes to the BACK (so we can pop them off;
    // this also removes size 1 domains, treated as a separate case,
    // rather than Hall sets).
    void sort_and_remove_singleton_domains();
  };

  const unsigned m_max_number_of_consecutive_failed_attempts;
  bool m_awaiting_initial_fill;
  boost::dynamic_bitset<> m_union_of_domains;

  // Ugly; should refactor to avoid these
  // (need a better domain intersection function).
  // We reuse these instead of reallocating new objects each time.
  // Note that UNLIKE std::set, copying, clearing, refilling
  // dynamic_bitsets is quite fast (because they're basically
  // just vectors).
  boost::dynamic_bitset<> m_union_of_domains_complement;
  boost::dynamic_bitset<> m_domain_mask_workset;

  std::set<ReusableStorageId> m_partition_ids;

  std::vector<ReusableStorageId> m_partition_ids_to_erase;
  std::vector<ReusableStorageId> m_new_partition_ids;

  ReusableStorage<Partition> m_partitions_storage;

  // Returns false if a nogood is found.
  bool fill_initial_partition(const DomainsAccessor& accessor);

  // A struct for dealing with control flow within the main reduce loop.
  struct WithinReduceLoopHandlePartitionRefill {
    bool nogood;
    bool should_skip;

    WithinReduceLoopHandlePartitionRefill(
        const DomainsAccessor& accessor, bool new_reduce_call,
        Partition& partition,
        std::vector<ReusableStorageId>& partition_ids_to_erase,
        ReusableStorageId partition_id,
        unsigned max_number_of_consecutive_failed_attempts);
  };

  // The "partition" reference might be invalidated by this function.
  // This assumes that we've already checked that there is a Hall set.
  ReductionResult within_reduce_loop_handle_hall_set_reduction(
      DomainsAccessor& accessor, Partition& partition,
      ReusableStorageId partition_id);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
