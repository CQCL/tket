// Copyright 2019-2022 Cambridge Quantum Computing
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
#include "DerivedGraphsCalculator.hpp"
#include <forward_list>

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Store the data for derived graphs D2,D3, to be calculated
 * and built up lazily as required. */
class DerivedGraphsContainer {
public:

  /** All information related to a particular vertex. */
  struct VertexData {
    /** The edges in D2 connected to this vertex, sorted by vertex number. */
    DerivedGraphsCalculator::NeighboursAndCounts depth_2_neighbours;

    /** The edges in D3 connected to this vertex, sorted by vertex number. */
    DerivedGraphsCalculator::NeighboursAndCounts depth_3_neighbours;

    /** Number of triangles in the original graph containing this vertex. */
    std::size_t triangle_count;

    // The following initially are unfilled; let the caller fill them lazily.
    // The "counts" vectors give the second elements in "neighbours",
    // i.e. the weights in the derived graphs, but sorted.
    // Thus we can quickly check if PV->TV may be possible
    // (simimlar to degree sequence domination).
    // The degree sequences will be the sorted vertex degrees
    // in the derived graphs D2, D3, giving another useful filter.

    std::vector<std::size_t> depth_2_counts;
    std::vector<std::size_t> depth_3_counts;
    std::vector<std::size_t> depth_2_degree_sequence;
    std::vector<std::size_t> depth_3_degree_sequence;
  };

  /** Calculated lazily. Return a reference to the data for a given vertex
   * in the pattern graph, with D2,D3 edges filled in.
   * The key point is that the reference remains valid for the lifetime
   * of this object, even as other data is calculated and stored.
   * @param v A vertex in the pattern graph.
   * @param pattern_ndata Data for the original pattern graph.
   * @return The data for v, with D2,D3 edges filled in.
   */
  VertexData& get_pattern_v_data_permanent_reference(VertexWSM v, const NeighboursData& pattern_ndata);

  /** Return data about a vertex in the target graph, with D2,D3 edges
   * filled in. The reference remains valid even if other data is added.
   * @param v A vertex in the target graph.
   * @param target_ndata Data for the original target graph.
   * @return The data for v, with D2,D3 edges filled in.
   */
  VertexData& get_target_v_data_permanent_reference(VertexWSM v, const NeighboursData& target_ndata);

private:

  // To perform the calculation of D2,D3 edges when required.
  DerivedGraphsCalculator m_calculator;

  typedef std::forward_list<VertexData> List;
  typedef List::iterator Iter;

  // Because we store iterators to a linked list,
  // the references are never invalidated.
  // The key is the vertex; the value is an iterator to an element
  // in the List object containing all the data.
  std::map<VertexWSM, Iter> m_pattern_iters;
  std::map<VertexWSM, Iter> m_target_iters;

  // Contains all the data.
  List m_list_container;

  // Given the input data, check if there is already an entry in the map.
  // If not, calculate the data; then return it.
  VertexData& get_vertex_data_permanent_reference(VertexWSM v, const NeighboursData& ndata, std::map<VertexWSM, Iter>& map);
};



}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
