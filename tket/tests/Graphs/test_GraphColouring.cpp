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

#include "EdgeSequence.hpp"
#include "EdgeSequenceColouringParameters.hpp"
#include "GraphTestingRoutines.hpp"
#include "Graphs/AdjacencyData.hpp"
#include "Graphs/GraphColouring.hpp"
#include "RandomGraphGeneration.hpp"
#include "RandomPlanarGraphs.hpp"
#include "Utils/RNG.hpp"

using std::map;
using std::vector;

namespace tket {
namespace graphs {
namespace tests {

SCENARIO("Test many colourings: random trees") {
  RNG rng;
  RandomTreeParameters params;
  AdjacencyData adjacency_data;
  EdgeSequenceColouringParameters colouring_parameters;
  EdgeSequence edge_sequence(adjacency_data, rng);

  for (size_t number_of_vertices = 5; number_of_vertices < 100;
       number_of_vertices += 20) {
    for (params.approx_number_of_children_per_node = 2;
         params.approx_number_of_children_per_node <= 4;
         ++params.approx_number_of_children_per_node) {
      for (params.approx_number_of_spawns = 5;
           params.approx_number_of_spawns <= 10;
           params.approx_number_of_spawns += 3) {
        for (int nn = 0; nn < 10; ++nn) {
          adjacency_data.clear(number_of_vertices);
          colouring_parameters.test_colourings(params, edge_sequence);
        }
      }
    }
  }
  CHECK(colouring_parameters.total_number_of_colourings == 4840);
}

SCENARIO("Test many colourings: random dense graphs") {
  RNG rng;
  RandomDenseGraphParameters params;
  AdjacencyData adjacency_data;
  EdgeSequenceColouringParameters colouring_parameters;
  EdgeSequence edge_sequence(adjacency_data, rng);

  for (size_t number_of_vertices = 2; number_of_vertices < 15;
       number_of_vertices += 5) {
    for (int nn = 0; nn < 10; ++nn) {
      adjacency_data.clear(number_of_vertices);
      colouring_parameters.test_colourings(params, edge_sequence);
    }
  }
  CHECK(colouring_parameters.total_number_of_colourings == 700);
}

SCENARIO("Test many colourings: random dense graphs with known colours") {
  RNG rng;
  RandomColouredDenseGraphParameters params;
  AdjacencyData adjacency_data;
  EdgeSequenceColouringParameters colouring_parameters;
  EdgeSequence edge_sequence(adjacency_data, rng);

  for (size_t number_of_vertices = 2; number_of_vertices < 15;
       number_of_vertices += 5) {
    for (params.max_number_of_colours_to_use = 1;
         params.max_number_of_colours_to_use < number_of_vertices;
         ++params.max_number_of_colours_to_use) {
      for (int nn = 0; nn < 10; ++nn) {
        adjacency_data.clear(number_of_vertices);
        colouring_parameters.test_colourings(params, edge_sequence);
      }
    }
  }
  CHECK(colouring_parameters.total_number_of_colourings == 4912);
}

SCENARIO("Test many colourings: random k-partite graphs") {
  RNG rng;
  RandomColouredKPartiteGraphParameters params;
  AdjacencyData adjacency_data;
  EdgeSequenceColouringParameters colouring_parameters;
  EdgeSequence edge_sequence(adjacency_data, rng);
  size_t total_number_of_edges = 0;
  size_t total_number_of_colourings = 0;

  for (params.number_of_vertex_sets = 1; params.number_of_vertex_sets < 5;
       ++params.number_of_vertex_sets) {
    for (params.number_of_vertices_in_each_set = 1;
         params.number_of_vertices_in_each_set < 5;
         ++params.number_of_vertices_in_each_set) {
      for (params.percentage_of_added_edges = 10;
           params.percentage_of_added_edges < 100;
           params.percentage_of_added_edges += 20) {
        params.add_edges(edge_sequence);

        total_number_of_edges += edge_sequence.edges.size();
        total_number_of_colourings +=
            colouring_parameters.test_colourings(params, edge_sequence);
      }
    }
  }
  CHECK(total_number_of_edges == 755);
  CHECK(total_number_of_colourings == 744);
}

SCENARIO("Test many colourings: random fibrous graphs") {
  RNG rng;
  RandomFibrousGraphParameters params;
  AdjacencyData adjacency_data;
  EdgeSequenceColouringParameters colouring_parameters;
  EdgeSequence edge_sequence(adjacency_data, rng);

  for (size_t number_of_vertices = 5; number_of_vertices < 50;
       number_of_vertices += 10) {
    for (params.number_of_strands = 1; params.number_of_strands < 10;
         ++params.number_of_strands) {
      for (int nn = 0; nn < 20; ++nn) {
        adjacency_data.clear(number_of_vertices);
        colouring_parameters.test_colourings(params, edge_sequence);
      }
    }
  }
  CHECK(colouring_parameters.total_number_of_colourings == 3413);
}

SCENARIO("Test many colourings: trivial graphs") {
  // Not actually used, no randomness
  RNG rng;
  AdjacencyData adjacency_data;
  EdgeSequenceColouringParameters colouring_parameters;
  EdgeSequence edge_sequence(adjacency_data, rng);

  INFO("Edgeless graphs");
  {
    EdgelessGraph params;

    // Very cheap to colour, so do lots
    for (size_t number_of_vertices = 1; number_of_vertices < 1000;
         ++number_of_vertices) {
      adjacency_data.clear(number_of_vertices);
      colouring_parameters.test_colourings(params, edge_sequence);
    }
  }
  INFO("Complete graphs (all edges filled)");
  {
    CompleteGraph params;
    for (size_t number_of_vertices = 5; number_of_vertices < 10;
         number_of_vertices += 4) {
      adjacency_data.clear(number_of_vertices);
      colouring_parameters.test_colourings(params, edge_sequence);
    }
  }
  CHECK(colouring_parameters.total_number_of_colourings == 46);
}

SCENARIO("Test fixed graph") {
  // Whenever a particular graph colouring fails somehow,
  // you can copy and paste it into here to add to the tests
  {
    const map<size_t, vector<size_t>> data = {
        // clang-format off
            { 0, { 2, 4, 5, 6, 7, 9, 12, } },
            { 1, { 6, 9, 10, 12, } },
            { 2, { 3, 4, 5, 7, 8, 9, 10, 11, } },
            { 3, { 5, 8, 9, 12, } },
            { 4, { 6, 7, 8, 10, 12, } },
            { 5, { 6, 8, 9, 10, 11, 12, } },
            { 6, { 7, 10, 11, 12, } },
            { 7, { 9, 11, } },
            { 8, { 10, 11, 12, } },
            { 9, { 10, 11, 12, } },
            { 10, { 12, } },
        // clang-format on
    };

    const vector<size_t> known_colouring{
        0, 2, 1, 0, 2, 2, 1, 3, 4, 4, 0, 0, 3,
    };

    GraphTestingRoutines::calculate_and_check_optimal_colouring(
        GraphColouringResult(known_colouring), AdjacencyData(data));
  }
  {
    const map<size_t, vector<size_t>> data = {
        // clang-format off
            { 0, { 1, 3, } },
            { 1, { 2, 3, 4, } },
            { 2, { 4, } },
            { 3, { 4, } },
        // clang-format on
    };
    const vector<size_t> known_colouring{0, 1, 2, 2, 0};

    GraphTestingRoutines::calculate_and_check_optimal_colouring(
        GraphColouringResult(known_colouring), AdjacencyData(data));
  }
  {
    const map<size_t, vector<size_t>> data = {
        // clang-format off
            { 1, { 3, 6, 8, 11, 13, 17, } },
            { 4, { 8, } },
            { 24, { }, }
        // clang-format on
    };

    const vector<size_t> known_colouring{
        0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0,
        1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    };

    GraphTestingRoutines::calculate_and_check_optimal_colouring(
        GraphColouringResult(known_colouring), AdjacencyData(data));
  }
  {
    map<size_t, vector<size_t>> data = {
        // clang-format off
            { 2, { 9, 10, } },
            { 10, { 13, } },
            { 14, {} }
        // clang-format on
    };
    const vector<size_t> known_colouring{
        0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0,
    };

    GraphTestingRoutines::calculate_and_check_optimal_colouring(
        GraphColouringResult(known_colouring), AdjacencyData(data));
  }
}

SCENARIO("Test random planar graphs") {
  const size_t grid_width = 20;
  const size_t max_number_of_regions = 50;
  RNG rng;
  RandomPlanarGraphs planar_graphs;
  size_t number_of_colourings = 0;

  for (int nn = 0; nn < 50; ++nn) {
    planar_graphs.reset(grid_width);
    size_t current_number_of_regions = 10 * grid_width * grid_width;

    for (size_t paranoid_loop_counter = 10 * grid_width * grid_width;
         paranoid_loop_counter > 0; --paranoid_loop_counter) {
      const auto new_number_of_regions = planar_graphs.merge_squares(rng);
      REQUIRE(new_number_of_regions <= current_number_of_regions);
      if (new_number_of_regions == current_number_of_regions) {
        continue;
      }
      current_number_of_regions = new_number_of_regions;
      if (current_number_of_regions < 5) {
        break;
      }
      if (current_number_of_regions > max_number_of_regions) {
        continue;
      }
      // NOW, at this stage, we have a new planar graph; so colour it.
      // Actually, this will miss colouring a few distinct graphs
      // with the same number of regions, but never mind.
      const auto raw_data = planar_graphs.get_region_data();
      const AdjacencyData adjacency_data(raw_data);
      const auto colouring =
          GraphColouringRoutines::get_colouring(adjacency_data);
      ++number_of_colourings;
      GraphTestingRoutines::require_valid_suboptimal_colouring(
          colouring, adjacency_data);

      INFO("number_of_colourings=" << number_of_colourings);
      // The Four Colour Theorem guarantees this (for an optimal colouring).
      REQUIRE(colouring.number_of_colours <= 4);
    }
  }
  CHECK(number_of_colourings == 2300);
}

// The Mycielski graph construction is a way to generate
// triangle-free graphs of high chromatic number,
// and hence good for testing colouring. (The max clique size is just 2).
// See https://en.wikipedia.org/wiki/Mycielskian
static AdjacencyData get_Mycielski_graph(const AdjacencyData& graph) {
  const size_t original_vertices = graph.get_number_of_vertices();
  const size_t vertex_w = 2 * original_vertices;

  AdjacencyData new_graph(2 * original_vertices + 1);
  for (size_t ii = 0; ii < original_vertices; ++ii) {
    REQUIRE(new_graph.add_edge(ii, vertex_w));
    const auto& original_neighbours = graph.get_neighbours(ii);

    for (size_t jj : original_neighbours) {
      // Add the original edges
      new_graph.add_edge(ii, jj);

      // Add u(i) -> v(j)
      new_graph.add_edge(ii, original_vertices + jj);
    }
  }
  return new_graph;
}

// Pass in an initial seed graph with known chromatic number
static void test_Mycielski_graph_sequence(
    AdjacencyData graph, size_t chromatic_number, size_t number_of_graphs) {
  const auto initial_number_of_vertices = graph.get_number_of_vertices();
  const auto initial_number_of_edges = graph.get_number_of_edges();
  const auto initial_chromatic_number = chromatic_number;

  for (size_t counter = 0; counter < number_of_graphs; ++counter) {
    if (counter != 0) {
      graph = get_Mycielski_graph(graph);
      ++chromatic_number;
    }
    const auto colouring = GraphColouringRoutines::get_colouring(graph);
    GraphTestingRoutines::require_valid_suboptimal_colouring(colouring, graph);

    INFO(
        "Counter=" << counter << ", v=" << graph.get_number_of_vertices()
                   << ", e=" << graph.get_number_of_edges()
                   << ". Initial graph: v=" << initial_number_of_vertices
                   << ", e=" << initial_number_of_edges << ", chromatic number "
                   << initial_chromatic_number);

    REQUIRE(colouring.number_of_colours == chromatic_number);
  }
}

SCENARIO("Test Mycielski graphs") {
  AdjacencyData graph(2);

  // Simple edge.
  graph.add_edge(0, 1);

  // Even though the graphs are large
  // (the last has 767 vertices and 22196 edges!),
  // our algorithm still colours them in a fraction of a second.
  // Some graphs have many vertices, but are not very dense.
  test_Mycielski_graph_sequence(graph, 2, 9);

  // V shape, still triangle-free.
  graph.clear(3);
  graph.add_edge(0, 1);
  graph.add_edge(1, 2);
  // Goes up to v=1023, e=35062 in a fraction of a second!
  test_Mycielski_graph_sequence(graph, 2, 9);
}

}  // namespace tests
}  // namespace graphs
}  // namespace tket
