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

#include <boost/graph/adjacency_list.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <vector>

#include "Graphs/DirectedGraph.hpp"
#include "Utils/UnitID.hpp"
#include "boost/range/iterator_range_core.hpp"

namespace tket {
namespace graphs {
namespace tests {
namespace test_DirectedGraph {

SCENARIO("Correct creation of DirectedGraph graphs") {
  GIVEN("Construct empty graph of nodes") {
    std::vector<Node> nodes{Node(3), Node(2), Node(5), Node(1)};
    DirectedGraph<Node> uidgraph(nodes);
    CHECK(uidgraph.n_nodes() == 4);
    CHECK(uidgraph.n_connected() == 0);
    CHECK(uidgraph.node_exists(Node(3)));
    CHECK(uidgraph.node_exists(Node(1)));
    CHECK(uidgraph.node_exists(Node(5)));
    CHECK(uidgraph.node_exists(Node(2)));
    CHECK(!uidgraph.node_exists(Node(4)));
    CHECK(!uidgraph.node_exists(Node(0)));
  }
  GIVEN("Construct Qubit graph from edges") {
    using Conn = DirectedGraph<Qubit>::Connection;

    std::vector<Conn> edges{
        {Qubit(0), Qubit(2)},
        {Qubit(3), Qubit(6)},
        {Qubit(6), Qubit(2)},
        {Qubit(2), Qubit(1)},
        {Qubit(1), Qubit(0)}};
    DirectedGraph<Qubit> uidgraph(edges);

    CHECK(uidgraph.n_nodes() == 5);
    CHECK(uidgraph.n_connected() == 5);

    CHECK(uidgraph.edge_exists(Qubit(0), Qubit(2)));
    CHECK(uidgraph.edge_exists(Qubit(3), Qubit(6)));
    CHECK(uidgraph.edge_exists(Qubit(6), Qubit(2)));
    CHECK(uidgraph.edge_exists(Qubit(2), Qubit(1)));
    CHECK(uidgraph.edge_exists(Qubit(1), Qubit(0)));
  }
  GIVEN("Construct graph using member functions") {
    std::vector<Node> uids = {Node(4), Node(1), Node(0), Node(1231)};
    DirectedGraph<Node> uidgraph;
    for (auto u : uids) uidgraph.add_node(u);

    uidgraph.add_connection(uids[0], uids[3], 3);
    uidgraph.add_connection(uids[2], uids[3], 0);

    CHECK(uidgraph.edge_exists(uids[0], uids[3]));
    CHECK(uidgraph.edge_exists(uids[2], uids[3]));
    CHECK(uidgraph.n_connections() == 2);
    CHECK(uidgraph.get_connection_weight(uids[0], uids[3]) == 3);
    CHECK(uidgraph.n_nodes() == 4);

    uidgraph.remove_connection(uids[0], uids[3]);
    uidgraph.remove_stray_nodes();

    CHECK(uidgraph.n_nodes() == 2);
    CHECK(uidgraph.n_connections() == 1);
  }
}

SCENARIO("Access underlying undirected connectivity") {
  GIVEN("some directed graph") {
    using Conn = DirectedGraph<Node>::Connection;
    std::vector<Conn> edges{{Node(0), Node(2)}, {Node(0), Node(4)},
                            {Node(3), Node(6)}, {Node(6), Node(3)},
                            {Node(6), Node(2)}, {Node(2), Node(1)},
                            {Node(1), Node(0)}};

    DirectedGraph<Node> uidgraph(edges);
    CHECK(uidgraph.n_connections() == edges.size());

    auto& g = uidgraph.get_undirected_connectivity();
    CHECK(boost::num_edges(g) == edges.size() - 1);
  }
}

SCENARIO("Disconnected graphs") {
  // TKET-1425
  GIVEN("a disconnected graph") {
    using Conn = DirectedGraph<Node>::Connection;
    std::vector<Conn> edges{{Node(0), Node(1)}, {Node(2), Node(3)}};
    DirectedGraph<Node> uidgraph(edges);
    CHECK(uidgraph.get_distance(Node(0), Node(0)) == 0);
    CHECK(uidgraph.get_distance(Node(2), Node(3)) == 1);
    CHECK_THROWS(uidgraph.get_distance(Node(0), Node(2)));
  }
}

}  // namespace test_DirectedGraph
}  // namespace tests
}  // namespace graphs
}  // namespace tket
