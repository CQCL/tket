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
#include <boost/property_map/property_map.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <catch2/catch_test_macros.hpp>
#include <map>
#include <vector>

#include "Graphs/TreeSearch.hpp"
#include "Graphs/Utils.hpp"

namespace tket {
namespace graphs {
namespace tests {
namespace test_TreeSearch {

template <typename VertexListT>
static auto get_graph() {
  using Graph =
      boost::adjacency_list<boost::listS, VertexListT, boost::bidirectionalS>;

  constexpr unsigned num_vertices = 20;
  Graph g(num_vertices);

  std::vector<utils::vertex<Graph>> VERTEX(num_vertices);
  for (unsigned i = 0; i < num_vertices; ++i) {
    VERTEX[i] = boost::vertex(i, g);
  }

  auto make_cycle = [&g, &VERTEX](unsigned from, unsigned to) {
    const unsigned delta = to - from;
    for (unsigned i = 0; i < delta; ++i) {
      unsigned i_plus_one = (i + 1) % delta;
      unsigned v1 = i + from;
      unsigned v2 = i_plus_one + from;
      boost::add_edge(VERTEX[v1], VERTEX[v2], g);
    }
  };

  // cycle for vertices 0...9
  make_cycle(0, 10);

  // cycle for vertices 10...14
  make_cycle(10, 15);

  // edge between two cycles
  boost::add_edge(VERTEX[9], VERTEX[10], g);
  boost::add_edge(VERTEX[12], VERTEX[8], g);

  // disconnected cycle 15...19
  make_cycle(15, 20);

  return g;
}

SCENARIO("BFS get_parents") {
  GIVEN("Using implicit vecS indices") {
    auto g = get_graph<boost::vecS>();

    auto bfs = run_bfs(0, g);
    auto parents = bfs.get_parents();

    for (unsigned i = 1; i < 10; ++i) {
      CHECK(parents[i] == (i - 1) % 10);
    }
    CHECK(parents[10] == 9);
    for (unsigned i = 11; i < 15; ++i) {
      CHECK(parents[i] == i - 1);
    }
    for (unsigned i = 15; i < 20; ++i) {
      CHECK(parents[i] == i);
    }
  }
  GIVEN("Using explicit indices") {
    auto g = get_graph<boost::listS>();
    unsigned num_vertices = boost::num_vertices(g);
    std::map<utils::vertex<decltype(g)>, unsigned> VERTEX;

    for (unsigned i = 0; i < num_vertices; ++i) {
      VERTEX[boost::vertex(i, g)] = num_vertices - 1 - i;
    }

    auto bfs = run_bfs(boost::vertex(0, g), g, VERTEX);
    auto parents = bfs.get_parents();

    auto is_parent = [&VERTEX, &parents, num_vertices](
                         unsigned node, unsigned child) {
      return VERTEX[parents[num_vertices - child - 1]] ==
             num_vertices - node - 1;
    };
    for (unsigned i = 1; i < 10; ++i) {
      CHECK(is_parent((i - 1) % 10, i));
      // CHECK(parents[num_vertices - i] == num_vertices - (i-1)%10);
    }
    CHECK(is_parent(9, 10));
    for (unsigned i = 11; i < 15; ++i) {
      CHECK(is_parent(i - 1, i));
      // CHECK(parents[num_vertices - i] == num_vertices - (i-1));
    }
    for (unsigned i = 15; i < 20; ++i) {
      CHECK(is_parent(i, i));
      // CHECK(parents[num_vertices - i] == num_vertices - i);
    }
  }
  GIVEN("Using temporary object") {
    auto g = get_graph<boost::vecS>();
    auto parents = run_bfs(0, g).get_parents();
    auto g2 = get_graph<boost::vecS>();
    auto bfs = run_bfs(0, g);
    auto parents2 = bfs.get_parents();

    const unsigned num_vertices = boost::num_vertices(g);

    for (unsigned i = 0; i < num_vertices; ++i) {
      CHECK(parents[i] == parents2[i]);
    }
  }
}

SCENARIO("BFS get_dists") {
  GIVEN("Using implicit vecS indices") {
    auto g = get_graph<boost::vecS>();

    auto dists = run_bfs(0, g).get_dists();
    CHECK(dists[0] == 0);
    CHECK(dists[3] == 3);
    CHECK(dists[12] == 12);

    auto bfs = run_bfs(4, g);
    const unsigned num_vertices = boost::num_vertices(g);
    for (unsigned i = 0; i < num_vertices; ++i) {
      CHECK(bfs.get_dists()[i] == bfs.get_dist(i));
    }
  }
}

SCENARIO("BFS path to root") {
  auto g = get_graph<boost::vecS>();

  auto path = run_bfs(11, g).path_to_root(4);
  decltype(path) res{4, 3, 2, 1, 0, 9, 8, 12, 11};

  for (unsigned i = 0; i < res.size(); ++i) {
    CHECK(res[i] == path[i]);
  }
  auto path2 = run_bfs(11, g).path_from_root(4);
  decltype(path) res2{11, 12, 8, 9, 0, 1, 2, 3, 4};

  for (unsigned i = 0; i < res2.size(); ++i) {
    CHECK(res2[i] == path2[i]);
  }
}

SCENARIO("BFS depth") {
  auto g = get_graph<boost::vecS>();

  auto bfs = run_bfs(4, g);

  CHECK(bfs.max_depth() == 10);
  CHECK(bfs.max_depth_vertex() == 14);
}

SCENARIO("DFS get_parents") {
  GIVEN("Using implicit vecS indices") {
    auto g = get_graph<boost::vecS>();

    auto dfs = run_dfs(0, g);
    auto parents = dfs.get_parents();

    for (unsigned i = 1; i < 10; ++i) {
      CHECK(parents[i] == (i - 1) % 10);
    }
    CHECK(parents[10] == 9);
    for (unsigned i = 11; i < 15; ++i) {
      CHECK(parents[i] == i - 1);
    }
    CHECK(parents[15] == 15);
    for (unsigned i = 16; i < 20; ++i) {
      CHECK(parents[i] == i - 1);
    }
  }
  GIVEN("Using explicit indices") {
    auto g = get_graph<boost::listS>();
    unsigned num_vertices = boost::num_vertices(g);
    std::map<utils::vertex<decltype(g)>, unsigned> VERTEX;

    for (unsigned i = 0; i < num_vertices; ++i) {
      VERTEX[boost::vertex(i, g)] = num_vertices - 1 - i;
    }

    auto dfs = run_dfs(boost::vertex(0, g), g, VERTEX);
    auto parents = dfs.get_parents();

    auto is_parent = [&VERTEX, &parents, num_vertices](
                         unsigned node, unsigned child) {
      return VERTEX[parents[num_vertices - child - 1]] ==
             num_vertices - node - 1;
    };
    for (unsigned i = 1; i < 10; ++i) {
      CHECK(is_parent((i - 1) % 10, i));
    }
    CHECK(is_parent(9, 10));
    for (unsigned i = 11; i < 15; ++i) {
      CHECK(is_parent(i - 1, i));
    }
    CHECK(is_parent(15, 15));
    for (unsigned i = 16; i < 20; ++i) {
      CHECK(is_parent(i - 1, i));
    }
  }
  GIVEN("Using temporary object") {
    auto g = get_graph<boost::vecS>();
    auto parents = run_dfs(0, g).get_parents();
    auto g2 = get_graph<boost::vecS>();
    auto dfs = run_dfs(0, g);
    auto parents2 = dfs.get_parents();

    const unsigned num_vertices = boost::num_vertices(g);

    for (unsigned i = 0; i < num_vertices; ++i) {
      CHECK(parents[i] == parents2[i]);
    }
  }
}

SCENARIO("DFS get_dists") {
  GIVEN("Using implicit vecS indices") {
    auto g = get_graph<boost::vecS>();

    auto dists = run_dfs(0, g).get_dists();
    CHECK(dists[0] == 0);
    CHECK(dists[3] == 3);
    CHECK(dists[12] == 12);

    auto dfs = run_dfs(4, g);
    const unsigned num_vertices = boost::num_vertices(g);
    for (unsigned i = 0; i < num_vertices; ++i) {
      CHECK(dfs.get_dists()[i] == dfs.get_dist(i));
    }
  }
}

SCENARIO("DFS path to root") {
  auto g = get_graph<boost::vecS>();

  auto path = run_dfs(11, g).path_to_root(4);
  decltype(path) res{4, 3, 2, 1, 0, 9, 8, 12, 11};

  for (unsigned i = 0; i < res.size(); ++i) {
    CHECK(res[i] == path[i]);
  }
  auto path2 = run_dfs(11, g).path_from_root(4);
  decltype(path) res2{11, 12, 8, 9, 0, 1, 2, 3, 4};

  for (unsigned i = 0; i < res2.size(); ++i) {
    CHECK(res2[i] == path2[i]);
  }
}

SCENARIO("DFS depth") {
  auto g = get_graph<boost::vecS>();

  auto dfs = run_dfs(4, g);

  CHECK(dfs.max_depth() == 10);
  CHECK(dfs.max_depth_vertex() == 14);
}

SCENARIO("longest simplest path") {
  auto g = get_graph<boost::vecS>();

  auto res = longest_simple_path(g);
  unsigned i = 0;
  for (unsigned v : res) {
    CHECK(v == i);
    ++i;
  }
  CHECK(res.size() == 15);
  CHECK(res[0] == 0);
}

}  // namespace test_TreeSearch
}  // namespace tests
}  // namespace graphs
}  // namespace tket
