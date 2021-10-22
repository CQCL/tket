// Copyright 2019-2021 Cambridge Quantum Computing
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
#include <catch2/catch.hpp>
#include <cmath>
#include <iostream>
#include <vector>

#include "Architecture/Architectures.hpp"
#include "Graphs/ArticulationPoints.hpp"
#include "Utils/UnitID.hpp"

namespace tket {
namespace graphs {
namespace test_Architectures {

SCENARIO("Testing FullyConnected") {
  GIVEN("A FullyConnected architecture constructed from an integer") {
    unsigned n_nodes = 10;
    node_vector_t nodes_vec =
        FullyConnected::get_nodes_canonical_order(n_nodes);
    node_set_t nodes(nodes_vec.begin(), nodes_vec.end());
    FullyConnected arch(n_nodes);

    REQUIRE(arch.n_uids() == nodes.size());
    for (const UnitID &uid : arch.get_all_uids()) {
      REQUIRE(nodes.count(Node(uid)));
    }
    for (auto [n1, n2] : arch.get_connections_vec()) {
      REQUIRE(nodes.count(n1));
      REQUIRE(nodes.count(n2));
    }

    for (unsigned i = 0; i < n_nodes; i++) {
      for (unsigned j = 0; j < n_nodes; j++) {
        if (i != j) {
          Node n1("fcNode", i);
          Node n2("fcNode", j);
          REQUIRE(arch.connection_exists(n1, n2));
        }
      }
    }
  }
  GIVEN("A FullyConnected architecture constructed from a node list") {
    Node n0{"a", 2};
    Node n1{"b", 0};
    Node n2{"c", 1};
    std::list<Node> nodes{n0, n1, n2};
    FullyConnected arch(nodes);
    REQUIRE(arch.n_uids() == nodes.size());
    REQUIRE(arch.get_degree(n1) == 2);
    node_set_t mindegs = arch.min_degree_uids();
    REQUIRE(mindegs.size() == 3);
    node_vector_t path1 = arch.get_path(n2, n0);
    REQUIRE(path1.size() == 2);
    node_vector_t path0 = arch.get_path(n2, n2);
    REQUIRE(path0.size() == 1);
    unsigned maxd = arch.get_max_depth(n1);
    REQUIRE(maxd == 1);
    unsigned outd = arch.get_out_degree(n0);
    REQUIRE(outd == 2);
    node_set_t neighbours = arch.get_neighbour_uids(n1);
    REQUIRE(neighbours.size() == 2);
    REQUIRE(!neighbours.contains(n1));
    std::set<std::pair<Node, Node>> edges = arch.get_connections_set();
    REQUIRE(edges.size() == 2 * 3);  // directed
    std::vector<std::size_t> dists = arch.get_distances(n2);
    REQUIRE(dists.size() == 2);
    REQUIRE(dists[0] == 1);
    REQUIRE(dists[1] == 1);
    unsigned d1 = arch.get_distance(n1, n2);
    REQUIRE(d1 == 1);
    unsigned d0 = arch.get_distance(n1, n1);
    REQUIRE(d0 == 0);
  }
}

SCENARIO("Testing RingArch") {
  using Arch = RingArch;
  unsigned n_nodes = 10;
  node_vector_t nodes_vec = Arch::get_nodes_canonical_order(n_nodes);
  node_set_t nodes(nodes_vec.begin(), nodes_vec.end());
  Arch arch(n_nodes);

  REQUIRE(arch.n_uids() == nodes.size());
  for (const UnitID &uid : arch.get_all_uids()) {
    REQUIRE(nodes.count(Node(uid)));
  }
  for (auto [n1, n2] : arch.get_connections_vec()) {
    REQUIRE(nodes.count(n1));
    REQUIRE(nodes.count(n2));
  }

  for (unsigned i = 0; i < n_nodes; i++) {
    Node n1("ringNode", i);
    Node n2("ringNode", (i + 1) % n_nodes);
    REQUIRE(arch.connection_exists(n1, n2));
  }
}

SCENARIO("Testing SquareGrid") {
  using Arch = SquareGrid;
  unsigned ver = 5;
  unsigned hor = 5;
  unsigned layer = 2;
  node_vector_t nodes_vec = Arch::get_nodes_canonical_order(ver, hor, layer);
  node_set_t nodes(nodes_vec.begin(), nodes_vec.end());
  Arch arch(ver, hor, layer);

  REQUIRE(nodes.size() == arch.n_uids());
  for (const UnitID &uid : arch.get_all_uids()) {
    REQUIRE(nodes.count(Node(uid)));
  }
  for (auto [n1, n2] : arch.get_connections_vec()) {
    REQUIRE(nodes.count(n1));
    REQUIRE(nodes.count(n2));
  }

  for (const Node &n : nodes) {
    int row = n.index()[0], col = n.index()[1], l = n.index()[2];
    for (const Node &neigh : arch.get_neighbour_uids(n)) {
      int row_neigh = neigh.index()[0], col_neigh = neigh.index()[1],
          l_neigh = neigh.index()[2];
      REQUIRE(
          abs(row - row_neigh) + abs(col - col_neigh) + abs(l - l_neigh) == 1);
    }
  }
}

SCENARIO("Diameters") {
  GIVEN("an empty architecture") {
    Architecture arc;
    CHECK_THROWS(arc.get_diameter());
  }
  GIVEN("a singleton architecture") {
    Architecture arc;
    arc.add_uid(Node(0));
    CHECK(arc.get_diameter() == 0);
  }
  GIVEN("a connected architecture") {
    Architecture arc({{0, 1}, {1, 2}, {2, 3}, {3, 0}});
    CHECK(arc.get_diameter() == 2);
  }
  GIVEN("a disconnected architecture") {
    // TKET-1425
    Architecture arc({{0, 1}, {1, 2}, {2, 0}, {3, 4}});
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

}  // namespace test_Architectures
}  // namespace graphs
}  // namespace tket
