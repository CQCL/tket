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
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "Mapping/LexicographicalComparison.hpp"

namespace tket {

SCENARIO("Test LexicographicalComparison::LexicographicalComparison") {
  GIVEN("Five Node Architecture, interacting nodes all in architecture.") {
    std::vector<Node> nodes = {
        Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
        Node("test_node", 3), Node("test_node", 4)};
    // n0 -- n1 -- n2
    //       |
    //       n3
    //       |
    //       n4
    Architecture architecture(
        {{nodes[0], nodes[1]},
         {nodes[1], nodes[2]},
         {nodes[1], nodes[3]},
         {nodes[3], nodes[4]}});
    interacting_nodes_t interacting_nodes = {
        {nodes[0], nodes[3]},
        {nodes[3], nodes[0]},
        {nodes[2], nodes[4]},
        {nodes[4], nodes[2]}};

    ArchitecturePtr sc = std::make_shared<Architecture>(architecture);

    LexicographicalComparison lc_test(sc, interacting_nodes);

    lexicographical_distances_t distances =
        lc_test.get_lexicographical_distances();
    REQUIRE(distances.size() == 3);
    REQUIRE(distances[0] == 2);
    REQUIRE(distances[1] == 2);
    REQUIRE(distances[2] == 0);
  }
  GIVEN("Three Node architecture, some interacting node not in architecture.") {
    std::vector<Node> nodes = {
        Node("test_node", 0), Node("test_node", 1), Node("test_node", 2)};
    Architecture architecture({{nodes[0], nodes[1]}, {nodes[1], nodes[2]}});
    ArchitecturePtr sa = std::make_shared<Architecture>(architecture);
    interacting_nodes_t interacting_nodes = {
        {nodes[0], Node("bad_node", 4)}, {Node("test_node", 3), nodes[0]}};
    REQUIRE_THROWS_AS(
        LexicographicalComparison(sa, interacting_nodes),
        LexicographicalComparisonError);
  }
}

SCENARIO("Test LexicographicalComparison::increment_distances") {
  GIVEN("Three Node Architecture, varying standard increments.") {
    std::vector<Node> nodes = {Node(0), Node(1), Node(2)};
    Architecture architecture({{nodes[0], nodes[1]}, {nodes[1], nodes[2]}});
    interacting_nodes_t interactions = {
        {nodes[0], nodes[2]}, {nodes[2], nodes[0]}};
    ArchitecturePtr sa = std::make_shared<Architecture>(architecture);
    LexicographicalComparison lc_test(sa, interactions);

    lexicographical_distances_t distances =
        lc_test.get_lexicographical_distances();
    REQUIRE(distances[0] == 2);
    REQUIRE(distances[1] == 0);

    std::pair<Node, Node> interaction = {nodes[0], nodes[2]};
    lc_test.increment_distances(distances, interaction, -2);
    REQUIRE(distances[0] == 0);
    REQUIRE(distances[1] == 0);

    REQUIRE_THROWS_AS(
        lc_test.increment_distances(distances, interaction, -2),
        LexicographicalComparisonError);

    interaction = {nodes[1], nodes[0]};
    lc_test.increment_distances(distances, interaction, 2);
    REQUIRE(distances[0] == 0);
    REQUIRE(distances[1] == 2);
  }
}

SCENARIO(
    "Test LexicographicalComparison::get_updated_distances, five node "
    "architecture") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("test_node", 3), Node("test_node", 4)};
  // n0 -- n1 -- n2
  //       |
  //       n3
  //       |
  // n4
  Architecture architecture(
      {{nodes[0], nodes[1]},
       {nodes[1], nodes[2]},
       {nodes[1], nodes[3]},
       {nodes[3], nodes[4]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  interacting_nodes_t interacting_nodes = {
      {nodes[0], nodes[3]},
      {nodes[3], nodes[0]},
      {nodes[2], nodes[4]},
      {nodes[4], nodes[2]}};

  LexicographicalComparison lc_test(shared_arc, interacting_nodes);
  GIVEN("Two identical legal swap, one node in interaction.") {
    swap_t swap_12 = {nodes[1], nodes[2]};
    swap_t swap_21 = {nodes[1], nodes[2]};
    lexicographical_distances_t distances_12 =
        lc_test.get_updated_distances(swap_12);
    REQUIRE(distances_12.size() == 3);
    REQUIRE(distances_12[0] == 0);
    REQUIRE(distances_12[1] == 4);
    REQUIRE(distances_12[2] == 0);
    REQUIRE(distances_12 == lc_test.get_updated_distances(swap_21));
  }
  GIVEN("Two identical legal swap, both node in interaction.") {
    swap_t swap_34 = {nodes[3], nodes[4]};
    swap_t swap_43 = {nodes[4], nodes[3]};
    lexicographical_distances_t distances_34 =
        lc_test.get_updated_distances(swap_34);
    REQUIRE(distances_34.size() == 3);
    REQUIRE(distances_34[0] == 2);
    REQUIRE(distances_34[1] == 2);
    REQUIRE(distances_34[2] == 0);
    REQUIRE(distances_34 == lc_test.get_updated_distances(swap_43));
  }
  GIVEN("Illegal swap.") {
    // illegal swap -> as Node not in architecture will return unchanged
    swap_t swap_illegal = {Node("bad_node", 0), Node("bad_node", 9)};
    lexicographical_distances_t distances_illegal =
        lc_test.get_updated_distances(swap_illegal);
    REQUIRE(distances_illegal == lc_test.get_lexicographical_distances());
  }
  GIVEN("Swap between two qubits in already adjacent interaction.") {
    interacting_nodes_t interacting = {
        {nodes[0], nodes[1]}, {nodes[3], nodes[4]}};
    LexicographicalComparison lc_in(shared_arc, interacting);
    swap_t swap_01 = {nodes[0], nodes[1]};
    swap_t swap_10 = {nodes[1], nodes[0]};
    swap_t swap_34 = {nodes[3], nodes[4]};
    swap_t swap_43 = {nodes[4], nodes[3]};
    lexicographical_distances_t distances_01 =
        lc_in.get_updated_distances(swap_01);
    lexicographical_distances_t distances_10 =
        lc_in.get_updated_distances(swap_10);
    lexicographical_distances_t distances_34 =
        lc_in.get_updated_distances(swap_34);
    lexicographical_distances_t distances_43 =
        lc_in.get_updated_distances(swap_43);
    lexicographical_distances_t base_distances =
        lc_in.get_lexicographical_distances();
    lexicographical_distances_t comp = {0, 0, 4};
    REQUIRE(base_distances == comp);
    REQUIRE(distances_01 == base_distances);
    REQUIRE(distances_10 == base_distances);
    REQUIRE(distances_34 == base_distances);
    REQUIRE(distances_43 == base_distances);
  }
}

SCENARIO("Test LexicographicalComparison::remove_swaps_lexicographical") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("test_node", 3), Node("test_node", 4)};
  // n0 -- n1 -- n2
  //       |
  //       n3
  //       |
  //       n4
  Architecture architecture(
      {{nodes[0], nodes[1]},
       {nodes[1], nodes[2]},
       {nodes[1], nodes[3]},
       {nodes[3], nodes[4]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  interacting_nodes_t interacting_nodes = {
      {nodes[0], nodes[3]},
      {nodes[3], nodes[0]},
      {nodes[2], nodes[4]},
      {nodes[4], nodes[2]}};

  LexicographicalComparison lc_test(shared_arc, interacting_nodes);
  GIVEN("Single Swap.") {
    swap_t swap_01 = {nodes[0], nodes[1]};
    swap_set_t candidate_swaps = {swap_01};
    lc_test.remove_swaps_lexicographical(candidate_swaps);
    REQUIRE(candidate_swaps.size() == 1);
    REQUIRE(*candidate_swaps.begin() == swap_01);
  }
  GIVEN("Two Swap, both identical.") {
    swap_t swap_01 = {nodes[0], nodes[1]};
    swap_t swap_10 = {nodes[1], nodes[0]};
    swap_set_t candidate_swaps = {swap_01, swap_10};
    lc_test.remove_swaps_lexicographical(candidate_swaps);
    REQUIRE(candidate_swaps.size() == 2);
  }
  GIVEN("Swap on all edges.") {
    swap_t swap_01 = {nodes[0], nodes[1]};
    swap_t swap_12 = {nodes[1], nodes[2]};
    swap_t swap_13 = {nodes[1], nodes[3]};
    swap_t swap_34 = {nodes[3], nodes[4]};
    swap_set_t candidate_swaps = {swap_01, swap_12, swap_13, swap_34};
    lc_test.remove_swaps_lexicographical(candidate_swaps);
    REQUIRE(candidate_swaps.size() == 1);
  }
}
}  // namespace tket
