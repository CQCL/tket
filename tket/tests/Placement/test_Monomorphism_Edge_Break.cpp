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

#include <catch2/catch_test_macros.hpp>

#include "Placement/Placement.hpp"

namespace tket {

SCENARIO("monomorphism edge break function using WSM") {
  Architecture triangle_with_leaf(std::vector<Architecture::Connection>{
      {Node(0), Node(1)},
      {Node(1), Node(2)},
      {Node(0), Node(2)},
      {Node(0), Node(3)}});

  const qubit_vector_t pattern_qubits{Qubit(0), Qubit(1), Qubit(2), Qubit(3)};
  QubitGraph q_graph(pattern_qubits);

  // Easier to read.
  std::map<std::string, unsigned> simple_map;
  for (unsigned ii = 0; ii <= 3; ++ii) {
    simple_map[std::string("q[") + std::to_string(ii) + "]"] = ii;
    simple_map[std::string("node[") + std::to_string(ii) + "]"] = ii;
  }

  const auto get_mappings_str = [&triangle_with_leaf, &q_graph, &simple_map](
                                    unsigned max_number) -> std::string {
    const auto mappings =
        monomorphism_edge_break(triangle_with_leaf, q_graph, max_number, 1000);
    CHECK(mappings.size() <= max_number);
    std::stringstream ss;
    unsigned ii = 0;
    for (const auto& bimap : mappings) {
      REQUIRE(bimap.size() == 4);
      ss << "Mapping[" << ii << "]:  { ";
      ++ii;
      unsigned jj = 0;
      for (const auto& pair : bimap.left) {
        REQUIRE(simple_map.at(pair.first.repr()) == jj);
        ++jj;
        ss << simple_map.at(pair.second.repr()) << " ";
      }
      ss << "}\n";
    }
    return ss.str();
  };
  CHECK(get_mappings_str(10) == "Mapping[0]:  { 0 1 2 3 }\n");

  q_graph.add_connection(Qubit(0), Qubit(1));
  CHECK(
      get_mappings_str(10) ==
      "Mapping[0]:  { 0 1 2 3 }\n"
      "Mapping[1]:  { 0 2 1 3 }\n"
      "Mapping[2]:  { 0 3 1 2 }\n"
      "Mapping[3]:  { 1 0 2 3 }\n"
      "Mapping[4]:  { 1 2 0 3 }\n"
      "Mapping[5]:  { 2 0 1 3 }\n"
      "Mapping[6]:  { 2 1 0 3 }\n"
      "Mapping[7]:  { 3 0 1 2 }\n");

  q_graph.add_connection(Qubit(1), Qubit(2));
  CHECK(
      get_mappings_str(20) ==
      "Mapping[0]:  { 0 1 2 3 }\n"
      "Mapping[1]:  { 0 2 1 3 }\n"
      "Mapping[2]:  { 1 0 2 3 }\n"
      "Mapping[3]:  { 1 0 3 2 }\n"
      "Mapping[4]:  { 1 2 0 3 }\n"
      "Mapping[5]:  { 2 0 1 3 }\n"
      "Mapping[6]:  { 2 0 3 1 }\n"
      "Mapping[7]:  { 2 1 0 3 }\n"
      "Mapping[8]:  { 3 0 1 2 }\n"
      "Mapping[9]:  { 3 0 2 1 }\n");

  q_graph.add_connection(Qubit(0), Qubit(2));
  CHECK(
      get_mappings_str(10) ==
      "Mapping[0]:  { 0 1 2 3 }\n"
      "Mapping[1]:  { 0 2 1 3 }\n"
      "Mapping[2]:  { 1 0 2 3 }\n"
      "Mapping[3]:  { 1 2 0 3 }\n"
      "Mapping[4]:  { 2 0 1 3 }\n"
      "Mapping[5]:  { 2 1 0 3 }\n");

  q_graph.add_connection(Qubit(0), Qubit(3));
  CHECK(
      get_mappings_str(10) ==
      "Mapping[0]:  { 0 1 2 3 }\n"
      "Mapping[1]:  { 0 2 1 3 }\n");

  q_graph.add_connection(Qubit(2), Qubit(3));
  // Even though the ORIGINAL problem has no solution,
  // the monomorphism edge break function erases pattern edges until it does.
  CHECK(
      get_mappings_str(10) ==
      "Mapping[0]:  { 1 3 0 2 }\n"
      "Mapping[1]:  { 2 3 0 1 }\n");
}

}  // namespace tket
