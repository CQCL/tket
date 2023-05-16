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

#include "tkwsm/Searching/ValueOrdering.hpp"

#include <algorithm>
#include <tkassert/Assert.hpp>
#include <tkrng/RNG.hpp>

#include "tkwsm/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/*
We are using "solution biased searching" as our heuristic. See the paper
"Sequential and Parallel Solution-Biased Search for Subgraph Algorithms":

https://link.springer.com/chapter/10.1007/978-3-030-19212-9_2

A free version available from a co-author's website:

http://blairarchibald.co.uk/resources/publications/cpaior19.pdf

The idea is, rather than always mapping into the target vertex with largest
possible degree, we should try lower degrees occasionally.
This is done simply by choosing a random vertex,
but letting the probability of choosing a vertex depend on the degree,
so that lower degrees are less likely.
*/

ValueOrdering::ValueOrdering() {
  /* 2^{-5} = 1/32 ~ 0.03. A choice with probability <3% is probably
     not very significant for the overall search.

    The probability of picking a vertex (the mass) is proportional to
      2^degree.

    Therefore if a vertex has degree <= max degree - 5, then it has weighting
    <= 1/32 ~ 3% of any vertex of maximum degree.

    Thus, including these vertices also would usually only make a difference
    in a very small number of cases.

    As explained in the linked paper, the heuristic generally chooses
    highest degree vertices first, but occasionally chooses
    slightly lower degrees.

    This is based on a very rough estimate that the expected number of
    solutions in the subtree is roughly C . 2^degree.

    However, this is all very rough and experimental (as the paper says);
    a better estimate is more like A.B^degree for unknown B>1, but that's
    still very rough.

    As the paper says, it generally gives better results than non-exponential
    formulae but there's no deep underlying theory here.

    The selection "5" is also an arbitrary choice; we need a cutoff to avoid
    wasting time on vanishingly unlikely choices, and 3% seems small enough
    that it's a reasonable cutoff. (OK, we could have 1 vertex of degree 10
    and 1000 vertices of degree 5. But even in that case, I think we'd rather
    just try the degree 10 vertex first).
  */
  m_entries_for_high_degree_vertices.resize(5);
  m_entries_for_high_degree_vertices.back().mass = 1;
  for (unsigned index = m_entries_for_high_degree_vertices.size() - 1;
       index > 0; --index) {
    m_entries_for_high_degree_vertices[index - 1].mass =
        2 * m_entries_for_high_degree_vertices[index].mass;
  }
}

void ValueOrdering::fill_entries_for_high_degree_vertices(
    const boost::dynamic_bitset<>& possible_values,
    const NeighboursData& target_ndata) {
  // Multipass algorithm; simple and efficient enough.
  // Pass one: find the max degree.
  size_t max_degree = 0;
  for (auto tv = possible_values.find_first(); tv < possible_values.size();
       tv = possible_values.find_next(tv)) {
    max_degree = std::max(max_degree, target_ndata.get_degree(tv));
  }
  TKET_ASSERT(max_degree > 0);

  // Pass two: record all those vertices with large enough degree.
  for (HighDegreeVerticesData& entry : m_entries_for_high_degree_vertices) {
    entry.vertices.clear();
  }
  for (auto tv = possible_values.find_first(); tv < possible_values.size();
       tv = possible_values.find_next(tv)) {
    const auto degree = target_ndata.get_degree(tv);
    if (degree + m_entries_for_high_degree_vertices.size() > max_degree) {
      const auto index = max_degree - degree;
      m_entries_for_high_degree_vertices[index].vertices.push_back(tv);
    }
  }
  TKET_ASSERT(!m_entries_for_high_degree_vertices[0].vertices.empty());
}

VertexWSM ValueOrdering::get_random_choice_from_data(RNG& rng) const {
  // We need probability proportional to the mass; so get the total mass.
  // Consider each vertex as having k copies, k corresponding to the mass,
  // list them all, and choose one uniformly.

  std::size_t mass_sum = 0;
  for (const auto& entry : m_entries_for_high_degree_vertices) {
    mass_sum += entry.vertices.size() * entry.mass;
  }
  TKET_ASSERT(mass_sum > 0);
  const auto index = rng.get_size_t(mass_sum - 1);

  // Now, choose the vertex corresponding to this index.
  std::size_t preceding_mass = 0;
  for (const auto& entry : m_entries_for_high_degree_vertices) {
    std::size_t mass_after_these_vertices =
        preceding_mass + entry.vertices.size() * entry.mass;
    if (mass_after_these_vertices < index) {
      preceding_mass = mass_after_these_vertices;
      continue;
    }
    for (auto vertex : entry.vertices) {
      size_t mass_after_this_vertex = preceding_mass + entry.mass;
      if (mass_after_this_vertex >= index) {
        return vertex;
      }
      preceding_mass = mass_after_this_vertex;
    }
  }

  // It's an error if we reach here, although a "harmless" one:
  // it just means that our calculation of the solution biased heuristic
  // is wrong.
  // So, if ever this is reached, and it's an "emergency"
  // (no time to find the bug), just comment it out temporarily
  // until the bug is properly fixed.
  // The worst that could happen is reduced performance; it can't
  // lead to incorrect results.
  TKET_ASSERT(false);
  return m_entries_for_high_degree_vertices.at(0).vertices.at(0);
}

VertexWSM ValueOrdering::get_target_value(
    const boost::dynamic_bitset<>& possible_values,
    const NeighboursData& target_ndata, RNG& rng) {
  const BitsetInformation bitset_info(possible_values);
  // Not size 0 or 1.
  TKET_ASSERT(!bitset_info.empty && !bitset_info.single_element);
  fill_entries_for_high_degree_vertices(possible_values, target_ndata);
  return get_random_choice_from_data(rng);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
