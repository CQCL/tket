// Copyright 2019-2023 Cambridge Quantum Computing
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

#pragma once
#include <string>
#include <tkrng/RNG.hpp>

#include "GraphGeneration.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// Consider the 2D integer lattice Z^2, with rectangles
// R = { (x,y) in Z^2  :  a <= x <= b,  c <= y <=d },
// and horizontal/vertical edges.
// It's easy to prove, by considering distances,
// that rectangles can only embed into another in the obvious ways,
// i.e. square-to-square, i.e. a translation by (dx,dy),
// followed possibly by a rotation and/or reflection.
struct SquareGrid {
  // Width N means the possible x co-ordinates are {0,1, ..., N},
  // i.e. there are N+1 vertical lines going across.
  unsigned width;
  unsigned height;
  std::vector<WeightWSM> horiz_weights;
  std::vector<WeightWSM> vert_weights;

  // Assumes that width, height have been filled already.
  // Merely resizes horiz_weights, vert_weights to the correct size.
  void resize_weight_vectors();

  std::string str() const;

  // Assumes that width, height have been filled already;
  // then resizes the other vectors appropriately and fills with random weights.
  void fill_weights(RNG& rng);

  GraphEdgeWeights get_graph_edge_weights() const;

  // Flip the grid and weights horizontally (i.e., the mirror axis is vertical).
  SquareGrid get_reflected_grid() const;

  // Rotate the grid 90 degrees.
  SquareGrid get_rotated_grid() const;

  // If translated by (dx,dy) and embedded into the other,
  // gets the scalar product.
  // Returns 0 if embedding is impossible.
  // NOTE: if we were really fancy we'd use Fast Fourier Transforms.
  WeightWSM get_scalar_product_translated_into_other(
      const SquareGrid& other, unsigned x_offset, unsigned y_offset) const;

  // Go through all possible translations of this grid,
  // and find the scalar product with the other grid,
  // then return the minimum.
  WeightWSM get_min_scalar_product_translated_into_other(
      const SquareGrid& other) const;

  // Allowing ALL translations, rotations, reflections, return the minimum
  // possible scalar product of an embedding of this grid into the other (or 0
  // if impossible).
  WeightWSM get_subgraph_isomorphism_min_scalar_product(
      const SquareGrid& other) const;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
