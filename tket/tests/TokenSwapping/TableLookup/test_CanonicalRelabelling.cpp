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
#include <numeric>
#include <set>

#include "PermutationTestUtils.hpp"
#include "TokenSwapping/CanonicalRelabelling.hpp"
#include "Utils/RNG.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

// Every element must represent the SAME mapping, up to an appropriate
// relabelling.
typedef vector<std::pair<VertexMapping, CanonicalRelabelling::Result>>
    EquivalentMappings;

// Everything in the OLD mapping does map to the expected vertex.
static void check_that_old_mapping_is_a_subset_of_expected(
    const VertexMapping& mapping,
    const CanonicalRelabelling::Result& relabelling,
    const std::array<unsigned, 6>& end_tokens) {
  for (const auto& orig_source_target_pair : mapping) {
    const auto& orig_source_v = orig_source_target_pair.first;
    const auto& orig_target_v = orig_source_target_pair.second;
    if (relabelling.old_to_new_vertices.count(orig_source_v) == 0) {
      // If this old vertex is unmentioned, it must be fixed.
      REQUIRE(orig_source_v == orig_target_v);
    } else {
      const auto new_source_v =
          relabelling.old_to_new_vertices.at(orig_source_v);
      const auto new_target_v =
          relabelling.old_to_new_vertices.at(orig_target_v);
      // end_tokens is the target->source mapping (the reverse of the usual).
      REQUIRE(end_tokens.at(new_target_v) == new_source_v);
    }
  }
}

// Everything in the expected new relabelled mapping agrees with the old
// mapping.
static void check_that_nonfixed_new_vertices_are_mentioned_in_old_mapping(
    const VertexMapping& mapping,
    const CanonicalRelabelling::Result& relabelling,
    const std::array<unsigned, 6>& end_tokens) {
  for (unsigned new_target_v = 0; new_target_v < end_tokens.size();
       ++new_target_v) {
    const auto new_source_v = end_tokens[new_target_v];
    if (new_source_v == new_target_v) {
      // Is it mentioned in the old mapping? If so, it must be fixed.
      if (new_source_v < relabelling.new_to_old_vertices.size()) {
        const auto old_source_v =
            relabelling.new_to_old_vertices.at(new_source_v);
        if (mapping.count(old_source_v) != 0) {
          // It IS mentioned, it MUST be fixed.
          REQUIRE(mapping.at(old_source_v) == old_source_v);
        }
      }
      continue;
    }
    // Different source, target, so the original mapping must mention this
    // (otherwise, the mapping would be incomplete).
    const auto old_source_v = relabelling.new_to_old_vertices.at(new_source_v);
    const auto old_target_v = relabelling.new_to_old_vertices.at(new_target_v);
    REQUIRE(mapping.at(old_source_v) == old_target_v);
  }
}

static void check_relabelling(const CanonicalRelabelling::Result& relabelling) {
  REQUIRE(
      relabelling.new_to_old_vertices.size() ==
      relabelling.old_to_new_vertices.size());
  REQUIRE(relabelling.new_to_old_vertices.size() >= 2);
  for (unsigned new_v = 0; new_v < relabelling.new_to_old_vertices.size();
       ++new_v) {
    const auto old_v = relabelling.new_to_old_vertices[new_v];
    REQUIRE(relabelling.old_to_new_vertices.at(old_v) == new_v);
  }
  for (const auto& old_new_pair : relabelling.old_to_new_vertices) {
    REQUIRE(
        relabelling.new_to_old_vertices.at(old_new_pair.second) ==
        old_new_pair.first);
  }
}

static void check_that_all_entries_have_the_same_permutation(
    unsigned permutation_hash, const EquivalentMappings& list) {
  REQUIRE(!list.empty());
  REQUIRE(permutation_hash >= 2);

  // end_tokens[i] tells us the SOURCE vertex of whatever token is now at vertex
  // i.
  const auto end_tokens =
      PermutationTestUtils::get_end_tokens_for_permutation(permutation_hash);

  for (const auto& entry : list) {
    const auto& mapping = entry.first;
    const auto& relabelling = entry.second;
    REQUIRE(relabelling.permutation_hash == permutation_hash);
    check_relabelling(relabelling);
    check_that_old_mapping_is_a_subset_of_expected(
        mapping, relabelling, end_tokens);
    check_that_nonfixed_new_vertices_are_mentioned_in_old_mapping(
        mapping, relabelling, end_tokens);
  }
}

// Create various random permutations on sets of size <= 6 of arbitrary labels,
// and see that the relabellings work.
SCENARIO("Relabelling test for random mappings") {
  const unsigned number_of_vertices = 5;
  vector<size_t> original_labels;

  // The generated mappings, together with the relabelling results.
  // The key is the permutation hash.
  std::map<unsigned, EquivalentMappings> entries;
  RNG rng;
  VertexMapping original_map;
  CanonicalRelabelling relabeller;

  for (unsigned nn = 0; nn < 200; ++nn) {
    original_map.clear();
    for (unsigned ii = 0; ii < number_of_vertices; ++ii) {
      original_map[rng.get_size_t(10000)];
    }
    original_labels.clear();
    for (const auto& entry : original_map) {
      original_labels.push_back(entry.first);
    }
    rng.do_shuffle(original_labels);
    {
      size_t ii = 0;
      for (auto& entry : original_map) {
        entry.second = original_labels[ii];
        ++ii;
      }
    }
    const auto result = relabeller(original_map);
    REQUIRE(!result.too_many_vertices);
    if (result.identity) {
      // Don't store identities.
      REQUIRE(all_tokens_home(original_map));
      REQUIRE(result.permutation_hash == 0);
      REQUIRE(result.old_to_new_vertices.empty());
      REQUIRE(result.new_to_old_vertices.empty());
    } else {
      REQUIRE(result.permutation_hash > 0);
      REQUIRE(result.old_to_new_vertices.size() == original_map.size());
      REQUIRE(result.new_to_old_vertices.size() == original_map.size());
      auto& list = entries[result.permutation_hash];
      list.push_back(std::make_pair(original_map, result));
    }
  }

  for (const auto& entry : entries) {
    check_that_all_entries_have_the_same_permutation(entry.first, entry.second);
  }
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
