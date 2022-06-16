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

#include <algorithm>
#include <catch2/catch_test_macros.hpp>

#include "../TestUtils/DebugFunctions.hpp"
#include "TokenSwapping/ExactMappingLookup.hpp"
#include "TokenSwapping/GeneralFunctions.hpp"
#include "TokenSwapping/VertexMappingFunctions.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

namespace {
struct ResultChecker {
  size_t failed_due_to_too_many_vertices = 0;
  size_t failed_due_to_table_missing_entry = 0;
  size_t success = 0;

  void check_failed_result(
      const ExactMappingLookup::Result& lookup_result,
      const VertexMapping& desired_mapping) {
    REQUIRE(!lookup_result.success);
    if (lookup_result.too_many_vertices) {
      CHECK(desired_mapping.size() >= 7);
      ++failed_due_to_too_many_vertices;
      return;
    }
    // There WERE enough edges. Why couldn't it find a solution?
    // The graph must have been too big.
    // The table should cover all 4-vertex mappings
    // (at least up to depth 12, and probably all).
    CHECK(desired_mapping.size() >= 5);
    ++failed_due_to_table_missing_entry;
  }

  void check_successful_result(
      const ExactMappingLookup::Result& lookup_result,
      const vector<Swap>& sorted_edges_vect, VertexMapping desired_mapping) {
    REQUIRE(lookup_result.success);
    ++success;
    // It succeeded. So, now we have to check it!
    CHECK(!lookup_result.too_many_vertices);

    // desired_mapping is a source->target mapping.
    // Interpret it to mean that mapping[i] = (current token on vertex i).
    // So initially, (token at i) = (target vertex).
    // Then, performing the swaps, all tokens should reach their home.
    for (const auto& swap : lookup_result.swaps) {
      REQUIRE(std::binary_search(
          sorted_edges_vect.cbegin(), sorted_edges_vect.cend(), swap));
      std::swap(desired_mapping[swap.first], desired_mapping[swap.second]);
    }
    CHECK(all_tokens_home(desired_mapping));
  }
};
}  // namespace

// We know that it succeeded and returned some swaps.
// Call it again with various max number of swaps limits.
static void recalculate_for_successful_problem_with_number_of_swaps_limits(
    const VertexMapping& desired_mapping, const vector<Swap>& edges_vect,
    const vector<Swap>& sorted_edges_vect, unsigned number_of_swaps,
    ExactMappingLookup& lookup, ResultChecker& checker) {
  for (unsigned max_number_of_swaps = 0; max_number_of_swaps < number_of_swaps;
       ++max_number_of_swaps) {
    const auto& lookup_result =
        lookup(desired_mapping, edges_vect, max_number_of_swaps);
    CHECK(!lookup_result.success);
  }
  for (unsigned max_number_of_swaps = number_of_swaps;
       max_number_of_swaps < number_of_swaps + 5; ++max_number_of_swaps) {
    const auto& lookup_result =
        lookup(desired_mapping, edges_vect, max_number_of_swaps);
    CHECK(lookup_result.success);
    CHECK(lookup_result.swaps.size() == number_of_swaps);
    checker.check_successful_result(
        lookup_result, sorted_edges_vect, desired_mapping);
  }
}

// A simple monotonic transformation, avoids contiguous vertices.
static unsigned get_vertex_number(unsigned ii) { return 10 * ii * (ii + 2); }

SCENARIO("Test exact mapping table lookup for wheel") {
  // A star is vertex 0, joined to 1,2,3,...,m.
  // A wheel also joins 1,2,...,m to make a cycle.
  VertexMapping desired_mapping;
  VertexMapping inverse_mapping;
  ExactMappingLookup lookup;

  // Maintain an unsorted vector, just in case sorting them makes a difference
  // (although it shouldn't).
  vector<Swap> all_edges;
  vector<Swap> all_edges_sorted;
  vector<size_t> vertices_used;
  ResultChecker checker;

  for (unsigned number_of_spokes = 3; number_of_spokes <= 6;
       ++number_of_spokes) {
    vertices_used.clear();
    vertices_used.push_back(0);
    all_edges.clear();
    for (unsigned ii = 1; ii <= number_of_spokes; ++ii) {
      const auto vv = get_vertex_number(ii);
      vertices_used.push_back(vv);
      all_edges.push_back(get_swap(0, vv));
    }
    // Complete the cycle on 1,2,...,m.
    all_edges.push_back(get_swap(vertices_used.back(), vertices_used[1]));
    for (unsigned ii = 1; ii < vertices_used.size(); ++ii) {
      all_edges.push_back(get_swap(vertices_used[ii - 1], vertices_used[ii]));
    }

    all_edges_sorted = all_edges;
    std::sort(all_edges_sorted.begin(), all_edges_sorted.end());
    desired_mapping.clear();

    // Set the SOURCE vertices.
    for (auto vv : vertices_used) {
      desired_mapping[vv];
    }
    for (int perm_counter = 0;;) {
      // Set the TARGET vertices.
      {
        unsigned ii = 0;
        for (auto& entry : desired_mapping) {
          entry.second = vertices_used[ii];
          ++ii;
        }
      }
      bool succeeded = false;
      unsigned number_of_swaps = 0;

      // We have a mapping. Try to look it up. Also, look up the inverse.
      inverse_mapping = get_reversed_map(desired_mapping);
      {
        // Care...because the result is stored internally,
        // another call to lookup will invalidate it!
        const auto& lookup_result = lookup(desired_mapping, all_edges);
        succeeded = lookup_result.success;
        if (lookup_result.success) {
          checker.check_successful_result(
              lookup_result, all_edges_sorted, desired_mapping);
          number_of_swaps = lookup_result.swaps.size();

          const auto& inverse_lookup_result =
              lookup(inverse_mapping, all_edges);
          CHECK(inverse_lookup_result.success);

          checker.check_successful_result(
              inverse_lookup_result, all_edges_sorted, inverse_mapping);
          CHECK(number_of_swaps == inverse_lookup_result.swaps.size());
        } else {
          // It failed. Why?
          checker.check_failed_result(lookup_result, desired_mapping);
          const auto& inverse_lookup_result =
              lookup(inverse_mapping, all_edges);
          checker.check_failed_result(inverse_lookup_result, inverse_mapping);
        }
      }

      if (succeeded) {
        recalculate_for_successful_problem_with_number_of_swaps_limits(
            desired_mapping, all_edges, all_edges_sorted, number_of_swaps,
            lookup, checker);
      }
      ++perm_counter;
      if (perm_counter > 10) {
        break;
      }
      if (!std::next_permutation(vertices_used.begin(), vertices_used.end())) {
        break;
      }
    }
  }

  CHECK(checker.failed_due_to_too_many_vertices == 22);
  CHECK(checker.failed_due_to_table_missing_entry == 0);
  CHECK(checker.success == 231);
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
