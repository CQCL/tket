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

#include <Graphs/Utils.hpp>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <catch2/catch_test_macros.hpp>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace tket {
namespace graphs {
namespace tests {
namespace test_GraphUtils {
using namespace tket::graphs::utils;

/* Random-ish graph used for testing
 *
 * VertexListT template type indicates the container type of the
 * graph to be created and returned (using boost container type tags)
 *   Options: boost::vecS, boost::listS, boost::setS, ...
 */
template <typename VertexListT>
static auto get_graph() {
  using Graph =
      boost::adjacency_list<boost::listS, VertexListT, boost::bidirectionalS>;
  constexpr unsigned num_vertices = 20;
  Graph g(num_vertices);

  std::vector<vertex<Graph>> VERTEX;
  for (unsigned i = 0; i < num_vertices; ++i) {
    VERTEX.push_back(boost::vertex(i, g));
  }

  // some edges for vertices 0...3
  for (unsigned m = 0; m < 4; ++m) {
    boost::add_edge(VERTEX[m], VERTEX[5], g);
    boost::add_edge(VERTEX[m], VERTEX[12], g);
    boost::add_edge(VERTEX[15], VERTEX[m], g);
  }

  // edges (m, 2) for 4 <= m <= 10
  for (unsigned m = 4; m <= 10; ++m) {
    boost::add_edge(VERTEX[m], VERTEX[2], g);
  }

  boost::add_edge(VERTEX[2], VERTEX[12], g);
  boost::add_edge(VERTEX[15], VERTEX[11], g);

  // cycle edges (m, m+1) for 12 <= m <= 18
  for (unsigned m = 12; m <= 18; ++m) {
    boost::add_edge(VERTEX[m], VERTEX[m + 1], g);
  }
  boost::add_edge(VERTEX[19], VERTEX[12], g);

  return g;
}

SCENARIO("The right class is instantiated") {
  GIVEN("vecS vertex container") {
    using Graph =
        boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>;

    Graph g(10);
    detail::graph_utils_impl<Graph> impl(g);

    REQUIRE(impl.to_index(1u) == 1u);
  }
  GIVEN("Explicit vertex index property") {
    using Graph =
        boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS>;
    using vertex_t = vertex<Graph>;

    Graph g(10);
    vertex_t third_vertex = boost::vertex(3, g);

    std::map<vertex_t, unsigned> pmap1;
    using PMap1 = decltype(pmap1);
    for (unsigned i = 0; i < 10; ++i) {
      auto v = boost::vertex(i, g);
      pmap1[v] = 9 - i;
    }
    detail::graph_utils_impl<Graph, PMap1> impl1(g, pmap1);

    REQUIRE(impl1.to_index(third_vertex) == 6u);
    impl1.remove_vertex(third_vertex);
    CHECK(impl1.to_index(boost::vertex(3, g)) == 5u);
    REQUIRE(impl1.to_index(boost::vertex(1, g)) == 7u);

    g = Graph(10);
    third_vertex = boost::vertex(3, g);
    std::map<vertex_t, std::string> pmap2;
    using PMap2 = decltype(pmap2);
    for (unsigned i = 0; i < 10; ++i) {
      auto v = boost::vertex(i, g);
      pmap2[v] = std::to_string(9 - i);
    }
    detail::graph_utils_impl<Graph, PMap2> impl2(g, pmap2);

    REQUIRE(pmap2[third_vertex] == "6");
    impl2.remove_vertex(boost::vertex(2, g));
    CHECK(pmap2[third_vertex] == "6");
    REQUIRE(pmap2[boost::vertex(1, g)] == "8");
  }
  GIVEN("Explicit vertex map with numeric indices") {
    using Graph =
        boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>;
    using vertex_t = vertex<Graph>;
    using Map = std::map<vertex_t, std::string>;

    Graph g(5);
    Map map;
    for (unsigned i = 0; i < 5; ++i) {
      map[i] = std::to_string(4 - i);
    }
    detail::graph_utils_impl_with_map<Graph, Map> impl(g, map);

    impl.remove_vertex(boost::vertex(2, g));
    CHECK(map[0] == "4");
    CHECK(map[1] == "3");
    CHECK(map[2] == "1");
    REQUIRE(map[3] == "0");
  }
  GIVEN("Explicit vertex map without numeric indices") {
    using Graph =
        boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS>;
    using vertex_t = vertex<Graph>;
    using Map = std::map<vertex_t, std::string>;

    Graph g(5);
    Map map;
    for (unsigned i = 0; i < 5; ++i) {
      auto v = boost::vertex(i, g);
      map[v] = std::to_string(4 - i);
    }
    detail::graph_utils_impl_with_map<Graph, Map> impl(g, map);

    vertex_t third_vertex = boost::vertex(3, g);
    CHECK(map[third_vertex] == "1");
    impl.remove_vertex(boost::vertex(2, g));
    REQUIRE(map[third_vertex] == "1");
  }
  GIVEN("Explicit vertex index property and explicit map") {
    using Graph =
        boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS>;
    using vertex_t = vertex<Graph>;
    using PMap = std::map<vertex_t, unsigned>;
    using Map = std::map<unsigned, std::string>;

    Graph g(5);
    Map map;
    PMap pmap;
    for (unsigned i = 0; i < 5; ++i) {
      auto v = boost::vertex(i, g);
      pmap[v] = i;
      map[i] = std::to_string(4 - i);
    }
    detail::graph_utils_impl_with_map<Graph, Map, PMap> impl(g, map, pmap);

    vertex_t third_vertex = boost::vertex(3, g);
    CHECK(pmap[third_vertex] == 3);
    REQUIRE(map[3] == "1");
    impl.remove_vertex(boost::vertex(2, g));
    CHECK(pmap[third_vertex] == 2);
    REQUIRE(map[2] == "1");
  }
}

SCENARIO("Using remove_vertex") {
  GIVEN("Removing a few vertices from a non-indexed graph") {
    auto g = get_graph<boost::listS>();

    std::set<vertex<decltype(g)>> removed;

    unsigned num_vertices = boost::num_vertices(g);
    // remove some vertices
    for (unsigned m : {16, 13, 9, 6, 2, 0}) {
      auto v = boost::vertex(m, g);
      boost::clear_vertex(v, g);
      remove_vertex(v, g);
      num_vertices--;
      CHECK(boost::num_vertices(g) == num_vertices);
    }
  }
  GIVEN(
      "Removing a few vertices from a graph with non-indexed property "
      "map") {
    auto g = get_graph<boost::listS>();
    std::map<vertex<decltype(g)>, std::string> pmap;
    for (unsigned m = 0; m < boost::num_vertices(g); ++m) {
      pmap[boost::vertex(m, g)] = std::to_string(m);
    }

    auto old_pmap = pmap;
    std::set<vertex<decltype(g)>> removed;

    // remove some vertices
    for (unsigned m : {16, 13, 9, 6, 2, 0}) {
      auto v = boost::vertex(m, g);
      boost::clear_vertex(v, g);
      remove_vertex(v, g, pmap);
      removed.insert(v);
    }

    for (auto [k, v] : old_pmap) {
      // either the vertex is in the new property map
      // or it was deleted
      CHECK(pmap.count(k) != removed.count(k));
    }
  }
  GIVEN("Removing a few vertices from an indexed graph") {
    auto g_vec = get_graph<boost::vecS>();
    auto g_list = get_graph<boost::listS>();
    std::map<vertex<decltype(g_list)>, unsigned> g_list_ind;
    for (unsigned m = 0; m < boost::num_vertices(g_list); ++m) {
      g_list_ind[boost::vertex(m, g_list)] = m;
    }

    // remove some vertices
    for (unsigned m : {16, 9, 6, 2}) {
      auto v_vec = boost::vertex(m, g_vec);
      boost::clear_vertex(v_vec, g_vec);
      remove_vertex(v_vec, g_vec);

      auto v_list = boost::vertex(m, g_list);
      boost::clear_vertex(v_list, g_list);
      remove_vertex(v_list, g_list, g_list_ind);
    }

    // check that both graph are identical
    for (unsigned m = 0; m < boost::num_vertices(g_list); ++m) {
      auto v_list = boost::vertex(m, g_list);
      auto v_vec = m;
      auto get_nth_adj = [](auto v, auto g, unsigned n) {
        auto [from, to] = boost::adjacent_vertices(v, g);
        return std::vector<decltype(*from)>(from, to)[n];
      };
      for (unsigned n = 0; n < boost::out_degree(v_list, g_list); ++n) {
        unsigned adj_vec = get_nth_adj(v_vec, g_vec, n);
        unsigned adj_list = g_list_ind[get_nth_adj(v_list, g_list, n)];
        CHECK(adj_vec == adj_list);
      }
    }
  }
}

SCENARIO("Using remove_edge") {
  GIVEN("Removing a few edges with remove_stray_vertices = false") {
    auto g = get_graph<boost::vecS>();
    auto n_edges = boost::num_edges(g);
    auto n_vertices = boost::num_vertices(g);
    remove_edge(boost::edge(4, 2, g).first, g);
    remove_edge(boost::edge(18, 19, g).first, g);
    n_edges -= 2;
    CHECK(boost::num_edges(g) == n_edges);
    REQUIRE(boost::num_vertices(g) == n_vertices);

    remove_edge(boost::edge(17, 18, g).first, g);
    --n_edges;
    CHECK(boost::num_edges(g) == n_edges);
    REQUIRE(boost::num_vertices(g) == n_vertices);
  }
  GIVEN("Removing a few edges with remove_stray_vertices = true") {
    auto g = get_graph<boost::vecS>();
    auto n_edges = boost::num_edges(g);
    auto n_vertices = boost::num_vertices(g);

    remove_edge(boost::edge(2, 12, g).first, g, true);
    remove_edge(boost::edge(17, 18, g).first, g, true);
    n_edges -= 2;
    CHECK(boost::num_edges(g) == n_edges);
    REQUIRE(boost::num_vertices(g) == n_vertices);

    remove_edge(boost::edge(16, 17, g).first, g, true);
    --n_edges;
    --n_vertices;
    CHECK(boost::num_edges(g) == n_edges);
    REQUIRE(boost::num_vertices(g) == n_vertices);
  }
}

SCENARIO("Using symmetrise") {
  GIVEN("A directed vecS graph") {
    using SymGraph =
        boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>;
    auto g = get_graph<boost::vecS>();
    auto sym = symmetrise<SymGraph>(g);

    for (auto e : boost::make_iterator_range(boost::edges(g))) {
      auto v = boost::source(e, g);
      auto w = boost::target(e, g);

      bool exists1 = boost::edge(v, w, sym).second;
      bool exists2 = boost::edge(w, v, sym).second;
      CHECK((exists1 && exists2));
    }
    for (auto e : boost::make_iterator_range(boost::edges(sym))) {
      auto v = boost::source(e, sym);
      auto w = boost::target(e, sym);

      bool exists1 = boost::edge(v, w, g).second;
      bool exists2 = boost::edge(w, v, g).second;
      CHECK((exists1 || exists2));
    }
  }
}

SCENARIO("Degree helper functions") {
  auto g = get_graph<boost::vecS>();
  REQUIRE(min_degree(g) == 1);
  REQUIRE(max_degree(g) == 11);

  auto max_deg_set = max_degree_nodes(g);
  CHECK(max_deg_set.count(2));
  REQUIRE(max_deg_set.size() == 1);

  auto min_deg_set = min_degree_nodes(g);
  CHECK(min_deg_set.count(4));
  CHECK(min_deg_set.count(6));
  CHECK(min_deg_set.count(7));
  CHECK(min_deg_set.count(8));
  CHECK(min_deg_set.count(9));
  CHECK(min_deg_set.count(10));
  CHECK(min_deg_set.count(11));
  REQUIRE(min_deg_set.size() == 7);
}

SCENARIO("Degree helper functions 2") {
  using Graph =
      boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>;
  using Edge = std::pair<int, int>;

  enum { A, B, C, D, E, N };
  const char* name = "ABCDE";
  Edge edge_array[] = {{A, B}, {A, D}, {C, A}, {D, C},
                       {C, E}, {B, D}, {D, E}, {A, E}};
  const int num_vertices = N;
  const int num_edges = sizeof(edge_array) / sizeof(edge_array[0]);

  Graph g(num_vertices);
  for (int i = 0; i < num_edges; ++i) {
    boost::add_edge(edge_array[i].first, edge_array[i].second, g);
  }

  GIVEN("max_degree test") {
    CHECK(max_degree(g) == 4);
    auto nodes = max_degree_nodes(g);
    CHECK(nodes.size() == 2);
    CHECK(nodes.find(D) != nodes.end());
    CHECK(nodes.find(A) != nodes.end());
  }

  GIVEN("min_degree test") {
    CHECK(min_degree(g) == 2);
    auto nodes = min_degree_nodes(g);
    CHECK(nodes.size() == 1);
    CHECK(nodes.find(B) != nodes.end());
  }
}

SCENARIO("Test symmetrise") {
  using Graph =
      boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>;
  using UndirGraph =
      boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
  using Edge = std::pair<int, int>;

  enum { A, B, C, D, E, N };
  const char* name = "ABCDE";
  Edge edge_array[] = {{A, B}, {A, D}, {C, A}, {D, C},
                       {C, E}, {B, D}, {D, E}, {A, E}};
  const int num_vertices = N;
  const int num_edges = sizeof(edge_array) / sizeof(edge_array[0]);

  Graph g(num_vertices);
  for (int i = 0; i < num_edges; ++i) {
    boost::add_edge(edge_array[i].first, edge_array[i].second, g);
  }

  auto g_sym = symmetrise<UndirGraph>(g);
  auto [e1, e1_exists] = boost::edge(B, A, g_sym);
  CHECK(e1_exists);
  auto [e2, e2_exists] = boost::edge(D, B, g_sym);
  CHECK(e2_exists);
  auto [e3, e3_exists] = boost::edge(B, D, g_sym);
  CHECK(e3_exists);
  auto [e4, e4_exists] = boost::edge(E, B, g_sym);
  CHECK(!e4_exists);
}

}  // namespace test_GraphUtils
}  // namespace tests
}  // namespace graphs
}  // namespace tket
