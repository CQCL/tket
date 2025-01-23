// Copyright Quantinuum
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
#include <array>
#include <catch2/catch_test_macros.hpp>

#include "SquareGridGeneration.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// Also test get_graph_edge_weights
static WeightWSM get_total_weights(const SquareGrid& grid) {
  WeightWSM total = 0;
  for (auto ww : grid.horiz_weights) {
    REQUIRE(ww > 0);
    total += ww;
  }
  for (auto ww : grid.vert_weights) {
    REQUIRE(ww > 0);
    total += ww;
  }
  const auto gmap = grid.get_graph_edge_weights();
  REQUIRE(gmap.size() == grid.horiz_weights.size() + grid.vert_weights.size());
  WeightWSM other_total = 0;
  for (const auto& entry : gmap) {
    other_total += entry.second;
  }
  REQUIRE(total == other_total);
  return total;
}

SCENARIO("square grid rotate 4 times equals identity") {
  RNG rng;
  for (int ii = 0; ii < 5; ++ii) {
    SquareGrid grid;
    grid.width = 4;
    grid.height = 7;
    grid.fill_weights(rng);
    const auto original_copy = grid;
    const auto total_weight = get_total_weights(grid);
    for (int nn = 0; nn < 4; ++nn) {
      grid = grid.get_rotated_grid();
      REQUIRE(total_weight == get_total_weights(grid));
    }
    REQUIRE(grid.horiz_weights == original_copy.horiz_weights);
    REQUIRE(grid.vert_weights == original_copy.vert_weights);
  }
}

SCENARIO("square grid refl/rotate twice equals identity") {
  RNG rng;
  for (int ii = 0; ii < 5; ++ii) {
    SquareGrid grid;
    grid.width = 4;
    grid.height = 4;
    grid.fill_weights(rng);
    const auto original_copy = grid;
    const auto total_weight = get_total_weights(grid);
    for (int nn = 0; nn < 2; ++nn) {
      grid = grid.get_reflected_grid();
      REQUIRE(total_weight == get_total_weights(grid));
      grid = grid.get_rotated_grid();
      REQUIRE(total_weight == get_total_weights(grid));
    }
    REQUIRE(grid.horiz_weights == original_copy.horiz_weights);
    REQUIRE(grid.vert_weights == original_copy.vert_weights);
  }
}

// Width 3, height 2 vertical edge indices:
//
//     +---+---+---+
//     1   3   5   7
//     +---+---+---+
//     0   2   4   6
//     +---+---+---+
//
// Horizontal indices:
//
//     +-6-+-7-+-8-+
//     |   |   |   |
//     +-3-+-4-+-5-+
//     |   |   |   |
//     +-0-+-1-+-2-+

SCENARIO("square grid reflection") {
  RNG rng;
  for (int ii = 0; ii < 5; ++ii) {
    SquareGrid grid;
    grid.width = 3;
    grid.height = 2;
    grid.fill_weights(rng);
    REQUIRE(grid.horiz_weights.size() == 9);
    REQUIRE(grid.vert_weights.size() == 8);

    const auto other_grid = grid.get_reflected_grid();
    REQUIRE(other_grid.width == grid.width);
    REQUIRE(other_grid.height == grid.height);
    REQUIRE(other_grid.horiz_weights.size() == grid.horiz_weights.size());
    REQUIRE(other_grid.vert_weights.size() == grid.vert_weights.size());
    REQUIRE(get_total_weights(grid) == get_total_weights(other_grid));

    // Which horiz edges are mapped into each other?
    const std::vector<std::pair<unsigned, unsigned>> horiz_equal_pairs{
        {0, 2}, {1, 1}, {3, 5}, {4, 4}, {6, 8}, {7, 7}};
    for (const auto& entry : horiz_equal_pairs) {
      REQUIRE(
          other_grid.horiz_weights.at(entry.first) ==
          grid.horiz_weights.at(entry.second));
      REQUIRE(
          other_grid.horiz_weights.at(entry.second) ==
          grid.horiz_weights.at(entry.first));
    }
    const std::vector<std::pair<unsigned, unsigned>> vert_equal_pairs{
        {0, 6}, {2, 4}, {1, 7}, {3, 5}};
    for (const auto& entry : vert_equal_pairs) {
      REQUIRE(
          other_grid.vert_weights.at(entry.first) ==
          grid.vert_weights.at(entry.second));
      REQUIRE(
          other_grid.vert_weights.at(entry.second) ==
          grid.vert_weights.at(entry.first));
    }
  }
}

// Width 2, height 4 vert/horiz indices:
//
//  +8+9+
//  3 7 11
//  +6+7+
//  2 6 10
//  +4+5+
//  1 5 9
//  +2+3+
//  0 4 8
//  +0+1+
//
// Width 4, height 2 indices:
//
//  + 8 + 9 + 10+ 11+
//  1   3   5   7   9
//  + 4 + 5 + 6 + 7 +
//  0   2   4   6   8
//  + 0 + 1 + 2 + 3 +
//
// ...and the original, rotated:
//
//  + 11+ 10+ 9 + 8 +
//  9   7   5   3   1
//  + 7 + 6 + 5 + 4 +
//  8   6   4   2   0
//  + 3 + 2 + 1 + 0 +
//

SCENARIO("square grid rotation") {
  RNG rng;
  for (int ii = 0; ii < 5; ++ii) {
    SquareGrid grid;
    grid.width = 2;
    grid.height = 4;
    grid.fill_weights(rng);
    REQUIRE(grid.horiz_weights.size() == 10);
    REQUIRE(grid.vert_weights.size() == 12);

    const auto other_grid = grid.get_rotated_grid();
    REQUIRE(other_grid.width == grid.height);
    REQUIRE(other_grid.height == grid.width);
    REQUIRE(other_grid.horiz_weights.size() == grid.vert_weights.size());
    REQUIRE(other_grid.vert_weights.size() == grid.horiz_weights.size());
    REQUIRE(get_total_weights(grid) == get_total_weights(other_grid));

    //     i:      an original horiz edge index
    // Element[i]: the new vert edge index it becomes
    const std::vector<unsigned> horiz_to_vert_data{8, 9, 6, 7, 4,
                                                   5, 2, 3, 0, 1};
    for (unsigned ii = 0; ii < horiz_to_vert_data.size(); ++ii) {
      REQUIRE(
          other_grid.vert_weights.at(horiz_to_vert_data[ii]) ==
          grid.horiz_weights.at(ii));
    }
    // element[i] is the new horiz edge index, for original vert edge index i
    const std::vector<unsigned> vert_to_horiz_data{3, 2, 1,  0,  7, 6,
                                                   5, 4, 11, 10, 9, 8};

    for (unsigned ii = 0; ii < vert_to_horiz_data.size(); ++ii) {
      REQUIRE(
          other_grid.horiz_weights.at(vert_to_horiz_data[ii]) ==
          grid.vert_weights.at(ii));
    }
  }
}

//  1x1 square; next to reflected; rotated
//
//  # b #   # b #   # B #
//  A   B   B   A   b   a
//  # a #   # a #   # A #

SCENARIO("square grid reflection, rotation on 1x1 square") {
  RNG rng;
  std::vector<unsigned> numbers;
  const auto add_numbers = [&numbers](const SquareGrid& gg) {
    numbers.push_back(gg.width);
    numbers.push_back(gg.height);
    numbers.push_back(gg.horiz_weights.size());
    numbers.push_back(gg.vert_weights.size());
    for (auto ww : gg.horiz_weights) {
      numbers.push_back(ww);
    }
    for (auto ww : gg.vert_weights) {
      numbers.push_back(ww);
    }
    numbers.push_back(0);
  };

  for (int ii = 0; ii < 5; ++ii) {
    SquareGrid grid;
    grid.width = 1;
    grid.height = 1;
    grid.fill_weights(rng);
    add_numbers(grid);
    const auto reflected_grid = grid.get_reflected_grid();
    add_numbers(reflected_grid);
    const auto rotated_grid = grid.get_rotated_grid();
    add_numbers(rotated_grid);
  }
  REQUIRE(
      numbers == std::vector<unsigned>{
                     1, 1, 2, 2, 8, 3, 7, 9, 0, 1, 1, 2, 2, 8, 3, 9, 7, 0, 1, 1,
                     2, 2, 7, 9, 3, 8, 0, 1, 1, 2, 2, 1, 4, 3, 1, 0, 1, 1, 2, 2,
                     1, 4, 1, 3, 0, 1, 1, 2, 2, 3, 1, 4, 1, 0, 1, 1, 2, 2, 5, 4,
                     3, 6, 0, 1, 1, 2, 2, 5, 4, 6, 3, 0, 1, 1, 2, 2, 3, 6, 4, 5,
                     0, 1, 1, 2, 2, 2, 5, 5, 8, 0, 1, 1, 2, 2, 2, 5, 8, 5, 0, 1,
                     1, 2, 2, 5, 8, 5, 2, 0, 1, 1, 2, 2, 5, 4, 7, 3, 0, 1, 1, 2,
                     2, 5, 4, 3, 7, 0, 1, 1, 2, 2, 7, 3, 4, 5, 0});
}

SCENARIO("square grid: check gdata conversion") {
  RNG rng;
  SquareGrid grid;
  grid.width = 2;
  grid.height = 3;
  std::vector<WeightWSM> sorted_weights;
  std::vector<WeightWSM> sorted_weights_again;

  for (grid.width = 1; grid.width < 10; ++grid.width) {
    for (grid.height = 1; grid.height < 10; ++grid.height) {
      grid.fill_weights(rng);
      sorted_weights.clear();
      for (auto ww : grid.horiz_weights) {
        sorted_weights.push_back(ww);
      }
      for (auto ww : grid.vert_weights) {
        sorted_weights.push_back(ww);
      }
      std::sort(sorted_weights.begin(), sorted_weights.end());
      const auto gdata = grid.get_graph_edge_weights();
      REQUIRE(gdata.size() == sorted_weights.size());
      sorted_weights_again.clear();
      for (const auto& entry : gdata) {
        sorted_weights_again.push_back(entry.second);
      }
      std::sort(sorted_weights_again.begin(), sorted_weights_again.end());
      REQUIRE(sorted_weights == sorted_weights_again);
    }
  }
}

SCENARIO(
    "square grid subgraph_isomorphism_min_scalar_product picks out reversed "
    "cycles") {
  SquareGrid cycle;
  cycle.width = 1;
  cycle.height = 1;
  // "random" increasing weights.
  const std::array<unsigned, 4> cycle_weights{1, 3, 7, 20};
  cycle.horiz_weights = {cycle_weights[0], cycle_weights[2]};
  cycle.vert_weights = {cycle_weights[3], cycle_weights[1]};
  // To get the minimum scalar product, the orders must be opposite.
  unsigned min_sc_prod = 0;
  for (unsigned ii = 0; ii < cycle_weights.size(); ++ii) {
    min_sc_prod +=
        cycle_weights[ii] * cycle_weights[cycle_weights.size() - 1 - ii];
  }

  SquareGrid big_grid;
  for (big_grid.width = 1; big_grid.width < 5; ++big_grid.width) {
    for (big_grid.height = 1; big_grid.height < 5; ++big_grid.height) {
      big_grid.resize_weight_vectors();
      for (unsigned dx = 0; dx < big_grid.width; ++dx) {
        for (unsigned dy = 0; dy < big_grid.height; ++dy) {
          // Make a copy of 1,2,3,4 within the edges,
          // starting at point (dx,dy).
          unsigned value = 9999;
          for (auto& ww : big_grid.horiz_weights) {
            ww = value--;
          }
          for (auto& ww : big_grid.vert_weights) {
            ww = value--;
          }
          const unsigned horiz_start = dy * big_grid.width + dx;
          big_grid.horiz_weights[horiz_start] = cycle_weights[0];
          big_grid.horiz_weights[horiz_start + big_grid.width] =
              cycle_weights[2];

          const unsigned vert_start = dx * big_grid.height + dy;
          big_grid.vert_weights[vert_start] = cycle_weights[3];
          big_grid.vert_weights[vert_start + big_grid.height] =
              cycle_weights[1];

          const auto sc_prod =
              cycle.get_subgraph_isomorphism_min_scalar_product(big_grid);
          REQUIRE(sc_prod == min_sc_prod);
        }
      }
    }
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
