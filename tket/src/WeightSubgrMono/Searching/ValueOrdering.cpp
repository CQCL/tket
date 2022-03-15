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

#include "WeightSubgrMono/Searching/ValueOrdering.hpp"

#include <algorithm>

#include "Utils/Assert.hpp"
#include "WeightSubgrMono/Searching/FixedData.hpp"
#include "WeightSubgrMono/Searching/SharedData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/*
We are using "solution biased searching" as our heuristic. See the paper
"Sequential and Parallel Solution-Biased Search for Subgraph Algorithms".
The idea is, rather than always mapping into the target vertex with largest
possible degree, we should try lower degrees occasionally.
This is done simply by choosing a random vertex,
but letting the probability of choosing a vertex depend on the degree,
so that lower degrees are less likely.
*/

ValueOrdering::ValueOrdering() {
  m_data.resize(5);
  m_data.back().mass = 1;
  for(unsigned index = m_data.size()-1; index > 0; --index) {
    m_data[index-1].mass = 2*m_data[index].mass;
  }
}


void ValueOrdering::fill_data(const std::set<VertexWSM>& possible_values, SharedData& shared_data) {
  // Pass one, just to find the max degree.
  size_t max_degree = 0;
  for (auto tv : possible_values) {
    max_degree = std::max(max_degree, 
        shared_data.fixed_data.target_neighbours_data.get_degree(tv));
  }  
  // Pass two: record all those vertices with large enough degree.
  for(auto& entry : m_data) {
    entry.vertices.clear();
  }

  for (auto tv : possible_values) {
    const auto degree =
        shared_data.fixed_data.target_neighbours_data.get_degree(tv);
    if(degree + m_data.size() > max_degree) {
      const auto index = max_degree - degree;
      m_data[index].vertices.push_back(tv);
    }
  }
  TKET_ASSERT(!m_data[0].vertices.empty());
}


VertexWSM ValueOrdering::get_random_choice_from_data(SharedData& shared_data) const {
  // We need probability proportional to the mass; so get the total mass.
  std::size_t mass_sum = 0;
  for(const auto& entry : m_data) {
    mass_sum += entry.vertices.size() * entry.mass;
  }
  TKET_ASSERT(mass_sum > 0);
  const auto index = shared_data.rng.get_size_t(mass_sum-1);

  // Now, choose the vertex corresponding to this index.
  std::size_t preceding_mass = 0;
  for(const auto& entry : m_data) {
    std::size_t mass_after_these_vertices = preceding_mass + entry.vertices.size() * entry.mass;
    if(mass_after_these_vertices < index) {
      preceding_mass = mass_after_these_vertices;
      continue;
    }
    for(auto vertex : entry.vertices) {
      size_t mass_after_this_vertex = preceding_mass + entry.mass;
      if(mass_after_this_vertex >= index) {
        return vertex;
      }
      preceding_mass = mass_after_this_vertex;
    }
  }

  /*
  // It's an error if we reach here, although a pretty harmless one:
  // it just means our calculation of the solution biased heuristic is wrong.
  // This assert can just be removed in an "emergency" until the bug is sorted.
  std::cerr << "\n\nm_data: [";
  for(const auto& entry : m_data) {
    std::cerr << " m=" << entry.mass << ", v=" << entry.vertices.size() << "; ";
  }
  std::cerr << "]\n\nindex=" << index << ";  mass_sum=" << mass_sum << "\n";
  */
  TKET_ASSERT(false);
  return m_data.at(0).vertices.at(0);
}


VertexWSM ValueOrdering::get_target_value(
    const std::set<VertexWSM>& possible_values, SharedData& shared_data) {
  TKET_ASSERT(possible_values.size() >= 2);
  fill_data(possible_values, shared_data);
  return get_random_choice_from_data(shared_data);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
