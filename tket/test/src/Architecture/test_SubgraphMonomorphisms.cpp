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

#include <catch2/catch_test_macros.hpp>

#include "tket/Architecture/SubgraphMonomorphisms.hpp"

namespace tket {

using Connection = Architecture::Connection;

static void test(
    const Architecture& pattern, const Architecture& target,
    const std::vector<std::string>& expected_mappings) {
  const ArchitectureMapping pattern_arch_mapping(pattern);
  const ArchitectureMapping target_arch_mapping(target);

  SubgraphMonomorphisms::Parameters parameters;
  parameters.timeout_ms = 1000;
  parameters.max_number_of_mappings = 10000;

  const SubgraphMonomorphisms solutions(
      pattern_arch_mapping, target_arch_mapping, parameters);
  CHECK(solutions.mappings.size() == expected_mappings.size());

  CHECK(solutions.time_taken_ms < 1000);

  // We want a canonical string; we don't want to rely
  // on the order of the nodes and vertex labels.
  std::vector<unsigned> target_v(pattern_arch_mapping.number_of_vertices());
  std::vector<std::string> mapping_strings;
  std::map<Node, unsigned> node_map;
  for (unsigned ii = 0; ii <= 10; ++ii) {
    node_map[Node(ii)] = ii;
  }
  for (const auto& mapping : solutions.mappings) {
    REQUIRE(mapping.size() == target_v.size());
    for (unsigned ii = 0; ii < mapping.size(); ++ii) {
      const auto& pv_node = pattern_arch_mapping.get_node(ii);
      const auto& tv_node = target_arch_mapping.get_node(mapping[ii]);
      target_v.at(node_map.at(pv_node)) = node_map.at(tv_node);
    }
    std::stringstream ss;
    for (unsigned tv : target_v) {
      ss << tv;
    }
    mapping_strings.push_back(ss.str());
  }
  std::sort(mapping_strings.begin(), mapping_strings.end());
  CHECK(mapping_strings == expected_mappings);
}

SCENARIO("get all embeddings") {
  // Diamond with extra edge
  Architecture pattern(std::vector<Connection>{
      {Node(0), Node(1)},
      {Node(0), Node(2)},
      {Node(1), Node(2)},
      {Node(1), Node(3)},
      {Node(3), Node(4)},
      {Node(2), Node(3)}});

  // Four triangles - beginning of Sierpinski triangle!
  Architecture target(std::vector<Connection>{
      {Node(0), Node(1)},
      {Node(1), Node(2)},
      {Node(2), Node(3)},
      {Node(3), Node(4)},
      {Node(4), Node(5)},
      {Node(5), Node(0)},
      {Node(1), Node(3)},
      {Node(3), Node(5)},
      {Node(5), Node(1)},
  });

  test(
      pattern, target,
      std::vector<std::string>{
          "01532", "01534", "05132", "05134", "21350", "21354", "23150",
          "23154", "43510", "43512", "45310", "45312"});

  // Now add some isolated vertices.
  pattern.add_node(Node(5));
  test(
      pattern, target,
      std::vector<std::string>{
          "015324", "015342", "051324", "051342", "213504", "213540", "231504",
          "231540", "435102", "435120", "453102", "453120"});

  // Too many pattern vertices now!
  pattern.add_node(Node(6));
  test(pattern, target, std::vector<std::string>{});

  target.add_node(Node(6));
  test(
      pattern, target,
      std::vector<std::string>{
          "0153246", "0153426", "0513246", "0513426", "2135046", "2135406",
          "2315046", "2315406", "4351026", "4351206", "4531026", "4531206"});
}

}  // namespace tket
