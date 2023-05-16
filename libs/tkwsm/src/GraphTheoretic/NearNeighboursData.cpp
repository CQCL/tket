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

#include "tkwsm/GraphTheoretic/NearNeighboursData.hpp"

#include <algorithm>
#include <tkassert/Assert.hpp>

#include "tkwsm/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

NearNeighboursData::NearNeighboursData(const NeighboursData& ndata)
    : m_ndata(ndata) {
  m_data.resize(ndata.get_number_of_nonisolated_vertices());
  m_degree_counts_work_vector.reserve(m_data.size());
}

std::size_t NearNeighboursData::get_number_of_vertices() const {
  return m_ndata.get_number_of_nonisolated_vertices();
}

const boost::dynamic_bitset<>&
NearNeighboursData::get_vertices_at_exact_distance(
    VertexWSM v, unsigned distance) {
  TKET_ASSERT(distance > 0);
  auto& list = m_data.at(v).vertices_at_exact_distance;
  if (list.empty()) {
    // Initialise with element[0], i.e. neighbours.
    list.emplace_back();
    list.back().resize(m_data.size());
    for (const auto& vertex_weight_pair :
         m_ndata.get_neighbours_and_weights(v)) {
      // The vertex must not already exist in the bitset.
      TKET_ASSERT(!list.back().test_set(vertex_weight_pair.first));
    }
  }
  const unsigned index = distance - 1;
  while (index >= list.size()) {
    // The list is nonempty, but not long enough yet; so grow it
    // (if necessary; if we hit an empty bitset
    // then all later ones will also be empty).

    // NOTE: easier to use "auto" here, as we DON'T know
    // the position type.
    auto vertex = list.back().find_first();

    if (vertex >= m_data.size()) {
      // The last element of the list is already an empty set,
      // so no point filling up with more empty sets.
      return list.back();
    }
    // Make a new blank bitset of the correct size, and fill it.
    list.emplace_back();
    const auto& previous_bitset = list[list.size() - 2];
    auto& current_bitset = list.back();

    // Go through all the vertices and find their neighbours, and add those.
    for (current_bitset.resize(m_data.size()); vertex < m_data.size();
         vertex = previous_bitset.find_next(vertex)) {
      for (const auto& vertex_weight_pair :
           m_ndata.get_neighbours_and_weights(vertex)) {
        current_bitset.set(vertex_weight_pair.first);
      }
      // Now, current_bitset has vertices as distance <= d
      // from the ROOT vertex, but they are distance 1 from vertices
      // at distance d-1. Therefore, some will have distance d-1 or d-2
      // from the root vertex; just erase these (easy for bitsets!)
      current_bitset -= previous_bitset;
      if (list.size() >= 3) {
        current_bitset -= list[list.size() - 3];
      }
      current_bitset.set(v, false);
    }
  }
  return list[index];
}

// Get the raw vector of degrees of all the given vertices.
static void fill_degrees_vector(
    const boost::dynamic_bitset<>& vertices, const NeighboursData& ndata,
    std::vector<std::size_t>& degree_counts_work_vector) {
  degree_counts_work_vector.clear();
  for (auto vertex = vertices.find_first(); vertex < vertices.size();
       vertex = vertices.find_next(vertex)) {
    degree_counts_work_vector.emplace_back(ndata.get_degree(vertex));
  }
}

// Given the raw vector of degrees (unsorted), turn it into a sorted list
// of (degree, count) pairs, within "degree_counts"
static void fill_degree_counts_vector(
    std::vector<std::size_t>& degree_counts_work_vector,
    FilterUtils::DegreeCounts& degree_counts) {
  TKET_ASSERT(degree_counts.empty());
  std::sort(degree_counts_work_vector.begin(), degree_counts_work_vector.end());
  auto citer = degree_counts_work_vector.cbegin();

  while (citer != degree_counts_work_vector.cend()) {
    const std::size_t degree = *citer;
    const auto current_citer = citer;

    // When do we next hit a strictly greater degree?
    citer = std::upper_bound(citer, degree_counts_work_vector.cend(), degree);

    // This is valid even if citer equals cend.
    degree_counts.emplace_back(degree, std::distance(current_citer, citer));
  }
}

const FilterUtils::DegreeCounts&
NearNeighboursData::get_degree_counts_at_exact_distance(
    VertexWSM v, unsigned distance) {
  std::vector<FilterUtils::DegreeCounts>& list =
      m_data.at(v).degree_counts_for_exact_distance;
  TKET_ASSERT(distance > 0);
  const unsigned index = distance - 1;
  if (index < list.size()) {
    return list[index];
  }
  // The index is too big. But, has it already hit empty degree counts?
  // (In which case, all future ones will also be empty?)
  if (!list.empty() && list.back().empty()) {
    return list.back();
  }

  // The index is too big. So, start expanding the list.
  for (unsigned distance_for_calculation = list.size() + 1;
       distance_for_calculation <= distance; ++distance_for_calculation) {
    TKET_ASSERT(list.size() + 1 == distance_for_calculation);
    list.resize(distance_for_calculation);

    fill_degrees_vector(
        get_vertices_at_exact_distance(v, distance_for_calculation), m_ndata,
        m_degree_counts_work_vector);

    if (m_degree_counts_work_vector.empty()) {
      return list.back();
    }
    fill_degree_counts_vector(m_degree_counts_work_vector, list.back());
  }
  TKET_ASSERT(list.size() == distance);
  return list.back();
}

const boost::dynamic_bitset<>& NearNeighboursData::get_vertices_up_to_distance(
    VertexWSM v, unsigned distance) {
  auto& list = m_data.at(v).vertices_up_to_distance;
  TKET_ASSERT(distance > 0);
  if (list.empty()) {
    list.resize(1);
    list[0] = get_vertices_at_exact_distance(v, 1);
  }
  const unsigned index = distance - 1;
  if (index < list.size()) {
    return list[index];
  }
  // Here, the list is nonempty, but not yet big enough.
  while (index >= list.size()) {
    list.emplace_back();
    list.back() = list[list.size() - 2];
    const auto& new_domain = get_vertices_at_exact_distance(v, list.size());
    if (new_domain.none()) {
      // No point in growing the list; there are no more vertices
      // at greater distances, so it's stopped changing.
      return list.back();
    }
    list.back() |= new_domain;
  }
  TKET_ASSERT(list.size() == distance);
  return list.back();
}

const FilterUtils::DegreeCounts&
NearNeighboursData::get_degree_counts_up_to_distance(
    VertexWSM v, unsigned distance) {
  // It seems easier to calculate degree counts from vertices,
  // than try to combine degree counts at different distances.
  auto& list = m_data.at(v).degree_counts_up_to_max_distance;
  TKET_ASSERT(distance > 0);
  if (list.empty()) {
    list.emplace_back(get_degree_counts_at_exact_distance(v, 1));
  }
  const unsigned index = distance - 1;
  if (index < list.size()) {
    return list[index];
  }
  while (index >= list.size()) {
    // Currently, data is for distances up to and including list.size().
    if (get_vertices_at_exact_distance(v, list.size() + 1).none()) {
      // No new vertices to add, so no point in continuing.
      return list.back();
    }
    // Get results for the next distance.
    list.emplace_back();
    fill_degrees_vector(
        get_vertices_up_to_distance(v, list.size()), m_ndata,
        m_degree_counts_work_vector);

    fill_degree_counts_vector(m_degree_counts_work_vector, list.back());
  }
  TKET_ASSERT(list.size() == distance);
  return list.back();
}

std::size_t NearNeighboursData::get_n_vertices_up_to_distance(
    VertexWSM v, unsigned max_distance) {
  if (max_distance == 0) {
    return 0;
  }
  // PERFORMANCE NOTE: should we cache? We suspect that
  // this is fast enough that it's not worth caching.
  return get_vertices_up_to_distance(v, max_distance).count();
}

std::size_t NearNeighboursData::get_n_vertices_at_exact_distance(
    VertexWSM v, unsigned distance) {
  if (distance == 0) {
    return 0;
  }
  // PERFORMANCE NOTE: should we cache? Probably not worth it.
  return get_vertices_at_exact_distance(v, distance).count();
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
