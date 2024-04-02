// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <boost/graph/adjacency_list.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <iostream>
#include <vector>

#include "tket/Architecture/Architecture.hpp"
#include "tket/Graphs/ArticulationPoints.hpp"

namespace tket {
namespace graphs {
namespace test_Architectures {

SCENARIO("Testing FullyConnected") {
  unsigned n_nodes = 10;
  FullyConnected arch(n_nodes);
  node_vector_t nodes_vec = arch.get_all_nodes_vec();
  node_set_t nodes(nodes_vec.begin(), nodes_vec.end());

  REQUIRE(arch.n_nodes() == nodes.size());
  for (const UnitID &uid : arch.nodes()) {
    REQUIRE(nodes.count(Node(uid)));
  }
  for (auto [n1, n2] : arch.get_all_edges_vec()) {
    REQUIRE(nodes.count(n1));
    REQUIRE(nodes.count(n2));
  }

  for (unsigned i = 0; i < n_nodes; i++) {
    for (unsigned j = 0; j < n_nodes; j++) {
      if (i != j) {
        Node n1("fcNode", i);
        Node n2("fcNode", j);
        REQUIRE(arch.edge_exists(n1, n2));
      }
    }
  }
  FullyConnected arch_named(2, "test_fc");
  REQUIRE(arch_named.get_all_nodes_vec()[0].reg_name() == "test_fc");
}

SCENARIO("Testing RingArch") {
  unsigned n_nodes = 10;
  RingArch arch(n_nodes);
  node_vector_t nodes_vec = arch.get_all_nodes_vec();
  node_set_t nodes(nodes_vec.begin(), nodes_vec.end());

  REQUIRE(arch.n_nodes() == nodes.size());
  for (const UnitID &uid : arch.nodes()) {
    REQUIRE(nodes.count(Node(uid)));
  }
  for (auto [n1, n2] : arch.get_all_edges_vec()) {
    REQUIRE(nodes.count(n1));
    REQUIRE(nodes.count(n2));
  }

  for (unsigned i = 0; i < n_nodes; i++) {
    Node n1("ringNode", i);
    Node n2("ringNode", (i + 1) % n_nodes);
    REQUIRE(arch.edge_exists(n1, n2));
  }

  RingArch arch_named(2, "test_ring");
  REQUIRE(arch_named.get_all_nodes_vec()[0].reg_name() == "test_ring");
}

SCENARIO("Testing SquareGrid") {
  unsigned ver = 5;
  unsigned hor = 5;
  unsigned layer = 2;
  SquareGrid arch(ver, hor, layer);
  node_vector_t nodes_vec = arch.get_all_nodes_vec();
  node_set_t nodes(nodes_vec.begin(), nodes_vec.end());

  REQUIRE(nodes.size() == arch.n_nodes());
  for (const UnitID &uid : arch.nodes()) {
    REQUIRE(nodes.count(Node(uid)));
  }
  for (auto [n1, n2] : arch.get_all_edges_vec()) {
    REQUIRE(nodes.count(n1));
    REQUIRE(nodes.count(n2));
  }

  for (const Node &n : nodes) {
    int row = n.index()[0], col = n.index()[1], l = n.index()[2];
    for (const Node &neigh : arch.get_neighbour_nodes(n)) {
      int row_neigh = neigh.index()[0], col_neigh = neigh.index()[1],
          l_neigh = neigh.index()[2];
      REQUIRE(
          std::abs(row - row_neigh) + std::abs(col - col_neigh) +
              std::abs(l - l_neigh) ==
          1);
    }
  }

  SquareGrid arch_named(2, 1, 1, "test_square_grid");
  REQUIRE(arch_named.get_all_nodes_vec()[0].reg_name() == "test_square_grid");
}

SCENARIO("Diameters") {
  GIVEN("an empty architecture") {
    Architecture arc;
    CHECK_THROWS(arc.get_diameter());
  }
  GIVEN("a singleton architecture") {
    Architecture arc;
    arc.add_node(Node(0));
    CHECK(arc.get_diameter() == 0);
  }
  GIVEN("a connected architecture") {
    Architecture arc(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(0)}});
    CHECK(arc.get_diameter() == 2);
  }
  GIVEN("a disconnected architecture") {
    // TKET-1425
    Architecture arc(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(0)},
         {Node(3), Node(4)}});
    CHECK_THROWS(arc.get_diameter());
  }
}

SCENARIO("connectivity") {
  GIVEN("simple architecture") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(2)},
         {Node(2), Node(3)}});
    MatrixXb connectivity(4, 4);
    connectivity << 0, 1, 1, 0,  // 0
        1, 0, 1, 0,              // 1
        1, 1, 0, 1,              // 2
        0, 0, 1, 0;              // 3

    REQUIRE(archi.get_connectivity() == connectivity);
  }
  GIVEN("connected architecture") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(2), Node(3)}});
    MatrixXb connectivity(4, 4);
    connectivity << 0, 1, 1, 1,  // 0
        1, 0, 1, 1,              // 1
        1, 1, 0, 1,              // 2
        1, 1, 1, 0;              // 3

    REQUIRE(archi.get_connectivity() == connectivity);
  }
}

SCENARIO("Test Architecture utility methods.") {
  GIVEN("Architecture::valid_operation, invalid and valid.") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}};
    Architecture architecture(edges);
    REQUIRE(!architecture.valid_operation({Node("test", 0)}));
    REQUIRE(architecture.valid_operation({Node(0)}));
    REQUIRE(architecture.valid_operation({Node(0), Node(1)}));
    REQUIRE(!architecture.valid_operation({Node(0), Node(2)}));
    REQUIRE(!architecture.valid_operation({Node(0), Node(1), Node(2)}));
  }
  GIVEN("Architecture::create_subarch") {
    std::vector<std::pair<unsigned, unsigned>> edges = {{0, 1}, {1, 2}};
    Architecture architecture(edges);
    Architecture subarc =
        architecture.create_subarch({Node(0), Node(1), Node(5)});
    REQUIRE(subarc.get_all_edges_vec().size() == 1);
  }
}
}  // namespace test_Architectures
}  // namespace graphs
}  // namespace tket
