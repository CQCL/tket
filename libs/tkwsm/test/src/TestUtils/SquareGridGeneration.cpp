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

#include "SquareGridGeneration.hpp"

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include <tkwsm/Common/GeneralUtils.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

void SquareGrid::resize_weight_vectors() {
  REQUIRE(width > 0);
  REQUIRE(height > 0);
  horiz_weights.resize(width * (height + 1));
  vert_weights.resize((width + 1) * height);
}

void SquareGrid::fill_weights(RNG& rng) {
  resize_weight_vectors();
  for (auto& ww : horiz_weights) {
    ww = rng.get_size_t(1, 9);
  }
  for (auto& ww : vert_weights) {
    ww = rng.get_size_t(1, 9);
  }
}

// The vertical edges are labelled as follows:
//
//     +---+---+
//    1|  3|   |5
//     +---+---+
//    0|  2|   |4
//     +---+---+
//
// ...i.e., vertical edge[0] is between points (0,0) and (0,1) here.
// vertical edge[5] is between points (2,1) and (2,2) here.

// The horizontal edges are labelled as follows:
//
//      3   4   5
//    +---+---+---+
//    |   |   |   |
//    +---+---+---+
//      0   1   2
//
// ...Thus, horiz edge[0] is between points (0,0) and (1,0) here,
// horiz edge[4] is between points (1,1) and (2,1) here.

SquareGrid SquareGrid::get_reflected_grid() const {
  SquareGrid other(*this);
  for (unsigned xx = 0; xx <= width; ++xx) {
    const unsigned vert_edge_start = xx * height;
    const unsigned other_vert_edge_start = (width - xx) * height;
    for (unsigned yy = 0; yy < height; ++yy) {
      other.vert_weights[other_vert_edge_start + yy] =
          vert_weights[vert_edge_start + yy];
    }
  }
  // Each horizontal slice is left in place, except that the indices reverse.
  for (unsigned yy = 0; yy <= height; ++yy) {
    const unsigned horiz_edge_start = yy * width;
    for (unsigned xx = 0; xx < width; ++xx) {
      other.horiz_weights[horiz_edge_start + xx] =
          horiz_weights[horiz_edge_start + (width - 1) - xx];
    }
  }
  return other;
}

// Width 3, height 2 vertical/horiz edge indices:
//
//     +-6-+-7-+-8-+
//     1   3   5   7
//     +-3-+-4-+-5-+
//     0   2   4   6
//     +-0-+-1-+-2-+
//
// Width 2, height 3; next to original indices, rotated:
//
//     +-6-+-7-+    +-7-+-6-+
//     2   5   8    8   5   2
//     +-4-+-5-+    +-5-+-4-+
//     1   4   7    7   4   1
//     +-2-+-3-+    +-3-+-2-+
//     0   3   6    6   3   0
//     +-0-+-1-+    +-1-+-0-+
//
// ...so we just interchange horiz/vert weights, then reflect

SquareGrid SquareGrid::get_rotated_grid() const {
  SquareGrid other;
  other.width = height;
  other.height = width;
  other.horiz_weights = vert_weights;
  other.vert_weights = horiz_weights;
  return other.get_reflected_grid();
}

GraphEdgeWeights SquareGrid::get_graph_edge_weights() const {
  GraphEdgeWeights map;
  for (unsigned xx = 0; xx <= width; ++xx) {
    for (unsigned yy = 0; yy <= height; ++yy) {
      // Start at point (x,y).
      const VertexWSM vv = yy * (width + 1) + xx;
      if (xx < width) {
        map[get_edge(vv, vv + 1)] = horiz_weights.at(yy * width + xx);
      }
      if (yy < height) {
        map[get_edge(vv, vv + width + 1)] = vert_weights.at(xx * height + yy);
      }
    }
  }
  return map;
}

WeightWSM SquareGrid::get_scalar_product_translated_into_other(
    const SquareGrid& other, unsigned dx, unsigned dy) const {
  if (width + dx > other.width || height + dy > other.height) {
    return 0;
  }
  WeightWSM total = 0;
  // Horiz weights.
  for (unsigned yy = 0; yy <= height; ++yy) {
    const unsigned horiz_edges_start = yy * width;
    const unsigned other_horiz_edges_start = (yy + dy) * other.width + dx;

    for (unsigned xx = 0; xx < width; ++xx) {
      const auto aa = horiz_weights.at(horiz_edges_start + xx);
      const auto bb = other.horiz_weights.at(other_horiz_edges_start + xx);
      total += aa * bb;
    }
  }
  // Vert weights.
  for (unsigned xx = 0; xx <= width; ++xx) {
    const unsigned vert_edges_start = xx * height;
    const unsigned other_vert_edges_start = (xx + dx) * other.height + dy;
    for (unsigned yy = 0; yy < height; ++yy) {
      const auto aa = vert_weights.at(vert_edges_start + yy);
      const auto bb = other.vert_weights.at(other_vert_edges_start + yy);
      total += aa * bb;
    }
  }
  return total;
}

WeightWSM SquareGrid::get_min_scalar_product_translated_into_other(
    const SquareGrid& other) const {
  if (width > other.width || height > other.height) {
    return 0;
  }
  WeightWSM total;
  set_maximum(total);
  const auto max_total = total;
  for (unsigned dx = 0; dx + width <= other.width; ++dx) {
    for (unsigned dy = 0; dy + height <= other.height; ++dy) {
      total = std::min(
          total, get_scalar_product_translated_into_other(other, dx, dy));
    }
  }
  REQUIRE(total < max_total);
  return total;
}

std::string SquareGrid::str() const {
  std::stringstream ss;
  ss << "\nWidth " << width << ", height " << height << "\n"
     << horiz_weights.size() << " horiz. edges, weights: [";
  for (auto ww : horiz_weights) {
    ss << ww << " ";
  }
  ss << "]\n" << vert_weights.size() << " vert. edges, weights: [";
  for (auto ww : vert_weights) {
    ss << ww << " ";
  }
  ss << "]\n";
  return ss.str();
}

WeightWSM SquareGrid::get_subgraph_isomorphism_min_scalar_product(
    const SquareGrid& other) const {
  WeightWSM total;
  set_maximum(total);
  const auto max_total = total;

  const auto add_scalar_product_with_rotations = [&other,
                                                  &total](SquareGrid grid) {
    for (int nn = 0; nn < 4; ++nn) {
      if (nn > 0) {
        grid = grid.get_rotated_grid();
      }
      const auto product =
          grid.get_min_scalar_product_translated_into_other(other);
      if (product > 0) {
        total = std::min(total, product);
      }
    }
  };
  add_scalar_product_with_rotations(*this);
  const auto reflected_grid = get_reflected_grid();
  add_scalar_product_with_rotations(reflected_grid);
  if (max_total == total) {
    return 0;
  }
  REQUIRE(total > 0);
  return total;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
