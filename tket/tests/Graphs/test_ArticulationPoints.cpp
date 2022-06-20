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
#include <boost/graph/adjacency_list.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <vector>

#include "Architecture/Architecture.hpp"
#include "Graphs/ArticulationPoints.hpp"

namespace tket {
namespace graphs {
namespace tests {
namespace test_ArticulationPoints {

using Edge = std::pair<unsigned, unsigned>;
using Vertex = unsigned;

static UndirectedConnGraph<Node> get_graph(std::vector<Edge> edges) {
  unsigned n_vertices = 0;
  for (auto [e1, e2] : edges) {
    if (e1 >= n_vertices) n_vertices = e1 + 1;
    if (e2 >= n_vertices) n_vertices = e2 + 1;
  }
  UndirectedConnGraph<Node> graph(n_vertices);
  for (unsigned i = 0; i < n_vertices; ++i) {
    graph[boost::vertex(i, graph)] = Node(i);
  }
  for (auto [e1, e2] : edges) {
    boost::add_edge(e1, e2, graph);
  }
  return graph;
}

SCENARIO("Build a BicomponentGraph") {
  UndirectedConnGraph<Node> graph;
  std::vector<std::vector<Vertex>> bicomps;
  std::vector<std::pair<unsigned, unsigned>> select_which;
  GIVEN("A simple graph") {
    graph = get_graph({{0, 1}, {1, 2}, {2, 3}, {1, 3}, {4, 0}});
    bicomps = {{0, 1}, {0, 4}, {1, 2, 3}};
    select_which = {{1, 2}, {4, 3}, {2, 3}};
  }
  GIVEN("A graph with 4 biconnected components") {
    bicomps = {{0, 1, 2}, {3, 4, 5, 6}, {7, 8, 9}, {10, 11}};
    std::vector<Vertex> links{3, 7, 4, 5};

    std::vector<Edge> edges;
    for (unsigned i = 0; i < links.size(); ++i) {
      auto &bicomp = bicomps[i];
      edges.push_back({bicomp[0], links[i]});
      for (unsigned j = 0; j < bicomp.size(); ++j) {
        for (unsigned k = j + 1; k < bicomp.size(); ++k) {
          edges.push_back({bicomp[j], bicomp[k]});
        }
      }
    }

    graph = get_graph(edges);
    select_which = {{0, 2}, {7, 4}, {10, 6}};
  }

  auto bg = detail::BicomponentGraph<Node>(graph);
  auto tester = detail::BicomponentGraphTester(&bg);
  auto &g = tester.get_graph();

  REQUIRE(boost::num_vertices(g) == tester.n_components());

  // Nodes in each comp must share a component vertex
  for (const std::vector<Vertex> &comp : bicomps) {
    std::set<unsigned> common_comps = tester.get_comps(Node(comp[0]));
    for (Vertex v : comp) {
      Node n(v);
      std::set<unsigned> new_set;
      std::set_intersection(
          common_comps.begin(), common_comps.end(), tester.get_comps(n).begin(),
          tester.get_comps(n).end(), std::inserter(new_set, new_set.begin()));
      common_comps = new_set;
    }
    REQUIRE(common_comps.size() == 1);
  }

  // select components and make sure selection is propagated
  for (auto [selecting, nb_selected] : select_which) {
    bg.select_comps(std::vector<Node>{Node(selecting)});
    bg.propagate_selected_comps();
    unsigned actual_nb_selected = 0;
    for (bool p : tester.get_selected_comps()) {
      actual_nb_selected += p;
    }
    REQUIRE(actual_nb_selected == nb_selected);
  }
}

SCENARIO("Run get_subgraph_aps") {
  UndirectedConnGraph<Node> graph = get_graph(
      {{0, 1},
       {1, 2},
       {0, 2},
       {2, 3},
       {3, 4},
       {2, 4},
       {4, 5},
       {5, 6},
       {6, 7},
       {7, 5},
       {4, 7},
       {7, 8},
       {8, 9},
       {7, 9}});
  UndirectedConnGraph<Node> subgraph =
      get_graph({{1, 2}, {0, 2}, {2, 3}, {6, 7}, {7, 5}});
  auto aps = get_subgraph_aps<Node>(graph, subgraph);
  auto expected_aps = std::set<Node>{Node(2), Node(4), Node(7)};

  REQUIRE(aps.size() == expected_aps.size());
  for (Node n : expected_aps) {
    REQUIRE(aps.count(n));
  }
}

SCENARIO("Find APs of disconnected nodes") {
  UndirectedConnGraph<Node> graph(2);
  for (unsigned i = 0; i < 2; ++i) {
    graph[boost::vertex(i, graph)] = Node(i);
  }

  UndirectedConnGraph<Node> subgraph(1);
  subgraph[boost::vertex(0, graph)] = Node(0);

  REQUIRE(get_subgraph_aps<Node>(graph, subgraph).size() == 0);
}

SCENARIO("Test APs in Architecture", "[architectures],[routing]") {
  Architecture arc({{0, 1}, {1, 2}, {2, 3}});
  node_set_t ap = arc.get_articulation_points();

  GIVEN("Removing node that preserves connected graph") {
    REQUIRE(ap.find(Node(0)) == ap.end());
    REQUIRE(ap.find(Node(3)) == ap.end());
  }
  GIVEN("Removing node doesn't preserve connected graph") {
    REQUIRE(ap.find(Node(1)) != ap.end());
    REQUIRE(ap.find(Node(2)) != ap.end());
  }

  Architecture arc2({{0, 1}, {0, 2}, {0, 3}, {2, 3}});
  ap = arc2.get_articulation_points();
  REQUIRE(ap.find(Node(0)) != ap.end());

  arc2.remove_node(Node(1));
  ap = arc2.get_articulation_points();
  REQUIRE(ap.find(Node(0)) == ap.end());
}

}  // namespace test_ArticulationPoints
}  // namespace tests
}  // namespace graphs
}  // namespace tket
