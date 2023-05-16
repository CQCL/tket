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

#include "tkwsm/Reducing/HallSetReduction.hpp"

#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/Searching/DomainsAccessor.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// NOTE:
// If a Hall set is detected in node[i],
// it might seem like a good idea to check node[i-1], node[i-2], ...
// also in the search branch, if the check is cheap.
// We need to experiment and see. Maybe it seems a little unlikely
// to be useful, since domain sizes are (mostly) monotonic,
// and it might seem that it would have been detected sooner.
// On the other hand, the Hall set detector might fail to detect
// many sets, purely down to a slight reordering of variables,
// so maybe the chance of finding the Hall set in previous nodes
// IS high enough to be worth trying.

HallSetReduction::HallSetReduction()
    : m_max_number_of_consecutive_failed_attempts(5),
      m_awaiting_initial_fill(true) {}

void HallSetReduction::clear() { m_awaiting_initial_fill = true; }

ReductionResult HallSetReduction::reduce(DomainsAccessor& accessor) {
  const bool new_reduce_call = m_awaiting_initial_fill;
  if (m_awaiting_initial_fill && !fill_initial_partition(accessor)) {
    return ReductionResult::NOGOOD;
  }
  if (m_partition_ids.empty()) {
    return ReductionResult::SUCCESS;
  }
  m_partition_ids_to_erase.clear();
  m_new_partition_ids.clear();

  ReductionResult result_to_return = ReductionResult::SUCCESS;

  for (auto id : m_partition_ids) {
    auto& partition = m_partitions_storage.get_nonconst_object(id);

    {
      const WithinReduceLoopHandlePartitionRefill refill_object(
          accessor, new_reduce_call, partition, m_partition_ids_to_erase, id,
          m_max_number_of_consecutive_failed_attempts);

      if (refill_object.nogood) {
        return ReductionResult::NOGOOD;
      }
      if (refill_object.should_skip) {
        continue;
      }
    }

    const auto search_result =
        partition.search_for_hall_set(accessor, m_union_of_domains);

    if (search_result == Partition::SearchResult::NOGOOD) {
      return ReductionResult::NOGOOD;
    }
    if (search_result == Partition::SearchResult::HALL_SET) {
      const auto hall_set_reduction_result =
          within_reduce_loop_handle_hall_set_reduction(accessor, partition, id);

      if (hall_set_reduction_result == ReductionResult::NOGOOD) {
        return ReductionResult::NOGOOD;
      }
      if (hall_set_reduction_result == ReductionResult::NEW_ASSIGNMENTS) {
        result_to_return = ReductionResult::NEW_ASSIGNMENTS;
      }
    }
  }
  for (auto id : m_partition_ids_to_erase) {
    TKET_ASSERT(m_partition_ids.erase(id) == 1);
    m_partitions_storage.release(id);
  }
  for (auto id : m_new_partition_ids) {
    TKET_ASSERT(m_partition_ids.insert(id).second);
  }
  return result_to_return;
}

HallSetReduction::WithinReduceLoopHandlePartitionRefill::
    WithinReduceLoopHandlePartitionRefill(
        const DomainsAccessor& accessor, bool new_reduce_call,
        Partition& partition,
        std::vector<ReusableStorageId>& partition_ids_to_erase,
        ReusableStorageId partition_id,
        unsigned max_number_of_consecutive_failed_attempts)
    : nogood(false), should_skip(false) {
  if (new_reduce_call) {
    return;
  }
  switch (partition.refill(accessor)) {
    case Partition::FillResult::NOGOOD:
      nogood = true;
      break;

    case Partition::FillResult::UNCHANGED:
      // No point in testing again, all data is unchanged
      // since the last check.
      should_skip = true;
      ++partition.attempt_count;
      if (partition.attempt_count >=
          max_number_of_consecutive_failed_attempts) {
        partition_ids_to_erase.push_back(partition_id);
      }
      break;

    case Partition::FillResult::CHANGED:
      partition.attempt_count = 0;
      break;

    case Partition::FillResult::ERASE_THIS_PARTITION:
      should_skip = true;
      partition_ids_to_erase.push_back(partition_id);
      break;
  }
}

ReductionResult HallSetReduction::within_reduce_loop_handle_hall_set_reduction(
    DomainsAccessor& accessor, Partition& partition,
    ReusableStorageId partition_id) {
  // First, erase the TV from all other domains in this partition.
  m_union_of_domains_complement = m_union_of_domains;
  m_union_of_domains_complement.flip();

  const auto reduction_data = partition.reduce_with_hall_set(
      accessor, m_union_of_domains, m_union_of_domains_complement,
      m_domain_mask_workset);

  if (reduction_data.result_to_return == ReductionResult::NOGOOD) {
    return ReductionResult::NOGOOD;
  }
  if (reduction_data.should_move_hall_pv_to_another_partition) {
    // Copy the Hall set data across to a new partition.
    // Of course, this may invalidate the "partition" reference!
    const auto new_id = m_partitions_storage.get_new_id();
    m_new_partition_ids.push_back(new_id);

    auto& new_partition = m_partitions_storage.get_nonconst_object(new_id);
    new_partition.attempt_count = 0;
    new_partition.domains_data.resize(m_union_of_domains.count());

    // CARE! The old partition reference might now be invalid!
    // So we have to get it again.
    const auto& old_partition_again =
        m_partitions_storage.get_object(partition_id);
    for (unsigned ii = 0; ii < new_partition.domains_data.size(); ++ii) {
      new_partition.domains_data[ii] = old_partition_again.domains_data.at(
          old_partition_again.domains_data.size() - ii - 1);
    }
  }
  if (reduction_data.should_erase_this_partition) {
    m_partition_ids_to_erase.push_back(partition_id);
  } else {
    // The old partition isn't actually being erased, so we're
    // retaining it; but we must delete the existing Hall set
    // (at the back), which is no longer relevant.
    auto& old_partition_again =
        m_partitions_storage.get_nonconst_object(partition_id);
    old_partition_again.domains_data.resize(
        old_partition_again.domains_data.size() - m_union_of_domains.count());
  }
  return reduction_data.result_to_return;
}

bool HallSetReduction::fill_initial_partition(const DomainsAccessor& accessor) {
  m_awaiting_initial_fill = false;
  for (auto id : m_partition_ids) {
    m_partitions_storage.release(id);
  }
  m_partition_ids.clear();
  const auto new_id = m_partitions_storage.get_new_id();

  switch (
      m_partitions_storage.get_nonconst_object(new_id).initial_fill(accessor)) {
    case Partition::FillResult::CHANGED:
      m_partition_ids.insert(new_id);
      break;
    case Partition::FillResult::ERASE_THIS_PARTITION:
      m_partitions_storage.release(new_id);
      break;
    case Partition::FillResult::NOGOOD:
      return false;
    case Partition::FillResult::UNCHANGED:
      // We just don't allow this to be returned for a new fill.
      TKET_ASSERT(false);
  }
  return true;
}

void HallSetReduction::Partition::sort_and_remove_singleton_domains() {
  std::sort(
      domains_data.begin(), domains_data.end(),
      [](const Data& lhs, const Data& rhs) {
        return lhs.domain_size > rhs.domain_size ||
               (lhs.domain_size == rhs.domain_size && lhs.pv > rhs.pv);
      });

  while (!domains_data.empty() && domains_data.back().domain_size < 2) {
    TKET_ASSERT(domains_data.back().domain_size == 1);
    domains_data.pop_back();
  }
}

HallSetReduction::Partition::FillResult
HallSetReduction::Partition::initial_fill(const DomainsAccessor& accessor) {
  domains_data.clear();
  for (VertexWSM pv : accessor.get_unassigned_pattern_vertices_superset()) {
    const auto domain_size = accessor.get_domain_size(pv);
    if (domain_size == 0) {
      return FillResult::NOGOOD;
    }
    if (domain_size == 1) {
      continue;
    }
    domains_data.emplace_back();
    domains_data.back().domain_size = domain_size;
    domains_data.back().pv = pv;
    attempt_count = 0;
  }
  if (domains_data.size() >= 3) {
    sort_and_remove_singleton_domains();
    if (domains_data.size() >= 3) {
      return FillResult::CHANGED;
    }
  }
  return FillResult::ERASE_THIS_PARTITION;
}

HallSetReduction::Partition::FillResult HallSetReduction::Partition::refill(
    const DomainsAccessor& accessor) {
  bool changed = false;
  for (auto& entry : domains_data) {
    const auto new_domain_size = accessor.get_domain_size(entry.pv);
    if (new_domain_size == 0) {
      return FillResult::NOGOOD;
    }
    if (new_domain_size != entry.domain_size) {
      changed = true;
      entry.domain_size = new_domain_size;
    }
  }
  if (!changed) {
    return FillResult::UNCHANGED;
  }
  if (domains_data.size() >= 3) {
    // You need at least 3 vertices for a useful Hall set.
    // (If there are vertices pv1, pv2 with
    // Dom(pv1) = Dom(pv2) = {a,b}, and no other domain including a or b,
    // then Hall set reduction can never reduce them further).
    // Note that singleton Hall sets Dom(pv) = {a} are treated
    // separately (i.e., by alldiff propagation).
    sort_and_remove_singleton_domains();
    if (domains_data.size() >= 3) {
      return FillResult::CHANGED;
    }
  }
  return FillResult::ERASE_THIS_PARTITION;
}

HallSetReduction::Partition::SearchResult
HallSetReduction::Partition::search_for_hall_set(
    const DomainsAccessor& accessor,
    boost::dynamic_bitset<>& union_of_domains) const {
  TKET_ASSERT(domains_data.size() > 2);

  // Note that only PV in this partition have domains which intersect
  // any possible combined union. Thus, a useful Hall set needs at least one
  // PV in this partition NOT to be included (for, they are the only PV which
  // could be reduced).

  if (domains_data.back().domain_size >= domains_data.size()) {
    // Even the smallest domain is too big to get a Hall set or nogood.
    return SearchResult::NO_HALL_SET;
  }

  union_of_domains = accessor.get_domain(domains_data.back().pv);

  for (std::size_t next_size = 2; next_size < domains_data.size();
       ++next_size) {
    // First, is it POSSIBLE that adding the next PV will create
    // a union of domains small enough? (The number of PVs must
    // hit or overtake the union size. We can get a guaranteed
    // lower bound on the union size. So, maybe we can prove
    // that ALL unions will be too big, for the order
    // in which we add them, and thus break out early).
    std::size_t union_size_lower_bound = union_of_domains.count();
    bool should_add_pv = false;

    // If we kept on adding PVs from this point, how big would the union
    // definitely get?
    for (auto hypothetical_size = next_size;
         hypothetical_size < domains_data.size(); ++hypothetical_size) {
      // Consider the hypothetical union with the next domain.
      union_size_lower_bound = std::max(
          union_size_lower_bound,
          domains_data[domains_data.size() - hypothetical_size].domain_size);

      if (union_size_lower_bound <= hypothetical_size) {
        // The union possibly could be small enough at this point.
        should_add_pv = true;
        break;
      }
    }
    if (!should_add_pv) {
      // At every stage, the union of domains was guaranteed to be too big.
      return SearchResult::NO_HALL_SET;
    }

    // Now, actually add the next PV.
    const VertexWSM& next_pv = domains_data[domains_data.size() - next_size].pv;
    union_of_domains |= accessor.get_domain(next_pv);
    const auto union_size = union_of_domains.count();
    if (union_size < next_size) {
      return SearchResult::NOGOOD;
    }
    if (union_size == next_size) {
      return SearchResult::HALL_SET;
    }
  }
  return SearchResult::NO_HALL_SET;
}

HallSetReduction::Partition::HallSetReductionData
HallSetReduction::Partition::reduce_with_hall_set(
    DomainsAccessor& accessor, const boost::dynamic_bitset<>& union_of_domains,
    const boost::dynamic_bitset<>& union_of_domains_complement,
    boost::dynamic_bitset<>& domain_mask_workset) {
  const unsigned number_of_pv = union_of_domains.count();
  TKET_ASSERT(number_of_pv > 1);
  TKET_ASSERT(number_of_pv < domains_data.size());

  HallSetReductionData reduction_data;
  reduction_data.changed = false;

  // As long as we have at least 3 vertices, there's the chance
  // of another Hall set.
  // Note that all OTHER vertices will have domains
  // disjoint from the Hall set at the end, so there's no point in keeping
  // the Hall set in THIS partition.
  reduction_data.should_move_hall_pv_to_another_partition = number_of_pv > 2;
  const unsigned number_of_vertices_to_reduce =
      domains_data.size() - number_of_pv;
  unsigned number_of_newly_assigned_vertices = 0;

  for (unsigned ii = 0; ii < number_of_vertices_to_reduce; ++ii) {
    // TODO: make an intersect function WITHOUT doing a swap,
    // to avoid this copy.
    domain_mask_workset = union_of_domains_complement;
    const auto intersect_result = accessor.intersect_domain_with_swap(
        domains_data[ii].pv, domain_mask_workset);

    if (intersect_result.reduction_result == ReductionResult::NOGOOD) {
      reduction_data.result_to_return = ReductionResult::NOGOOD;
      return reduction_data;
    }
    // Conveniently, the old domain size is already stored.
    if (intersect_result.changed) {
      reduction_data.changed = true;
      domains_data[ii].domain_size = intersect_result.new_domain_size;
    }
    // We continue to the end, just for this partition,
    // so at least everything is fully reduced before we go round the loop
    // again.
    if (intersect_result.reduction_result == ReductionResult::NEW_ASSIGNMENTS) {
      ++number_of_newly_assigned_vertices;
    }
  }

  if (number_of_newly_assigned_vertices > 0) {
    reduction_data.result_to_return = ReductionResult::NEW_ASSIGNMENTS;
  } else {
    reduction_data.result_to_return = ReductionResult::SUCCESS;
  }

  // How many vertices would be left in THIS partition
  // after removing the newly assigned ones (and the old Hall set, which will
  // be removed whatever happens)?
  reduction_data.should_erase_this_partition =
      domains_data.size() - number_of_pv - number_of_newly_assigned_vertices <
      3;

  return reduction_data;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
