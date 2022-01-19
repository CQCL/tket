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

#pragma once

#include <limits>

#include "Graphs/GraphColouring.hpp"

namespace tket {
namespace graphs {
namespace tests {

struct EdgeSequence;

/**
 * Base class for any set of parameters
 * to generate a random graph of particular type,
 * by repeatedly adding edges in a sequence.
 * This class generates the edges, but doesn't store them itself.
 * For testing, we want to generate many different kinds of graphs.
 */
struct RandomGraphParameters {
  /**
   * At the moment when "add_edges" is called, it may be possible
   * to create a valid but possibly suboptimal colouring easily
   * (for the full graph with all edges).
   * If empty, then no colouring has been set.
   * If nonempty, it is a valid but possibly suboptimal colouring.
   */
  GraphColouringResult known_colouring;

  /**
   * There are some types of graph, e.g. trees, where we don't immediately
   * know an actual colouring, but we can prove that it's always possible
   * with a certain number of colours.
   * The graph is guaranteed to be colourable
   * with at most this number of colours.
   */
  std::size_t max_chromatic_number = std::numeric_limits<std::size_t>::max();

  /**
   * The caller must ensure that "edge_sequence" is correctly initialised
   * with the number of vertices, etc.
   * For things like, e.g. trees, the class assumes that
   * the graph is empty at the start. However, the caller
   * might choose to ignore this, in which case you can get messy graphs
   * with trees, etc. etc. stuck on top of each other, which of course
   * might no longer be trees (and it will mess up the known_colouring
   * and max_chromatic_number, if any).
   * @param edge_sequence The storage for the graph data, and giving access to
   * an RNG.
   * @return false if we should not call again.
   */
  virtual bool add_edges(EdgeSequence& edge_sequence) = 0;
};

/**
 * I just made this name up, "fibrous" is not standard.
 * Grow random single "fibres" (or "strands"), which can join up at the end
 * to create cycles. Extra crossings with itself and other fibres
 * may occur.
 */
struct RandomFibrousGraphParameters : public RandomGraphParameters {
  std::size_t number_of_strands = 10;

  /**
   * Keep adding edges one-by-one to the strand until a coin flip
   * with this percentage of success fails.
   */
  std::size_t percentage_for_each_strand_to_grow = 80;

  /**
   * At the end, decide whether to join the last vertex in this strand
   * to the first, to ensure a cycle. (Of course, many shorter cycles
   * may already exist, due to chance overlappings).
   */
  std::size_t percentage_for_strand_to_become_a_cycle = 50;

  virtual bool add_edges(EdgeSequence& edge_sequence) override;
};

/**
 * Generates a tree, so ensure that it doesn't ever create a cycle,
 * as well as being connected. (Actually, it will be a single tree,
 * plus extra isolated vertices).
 */
struct RandomTreeParameters : public RandomGraphParameters {
  std::size_t approx_number_of_children_per_node = 2;
  std::size_t approx_number_of_spawns = 10;

  virtual bool add_edges(EdgeSequence& edge_sequence) override;
};

/**
 * Simply add random i-j edges one-by-one,
 * giving up when it can't easily add any more
 * (and thus, has probably reached a fairly dense graph).
 */
struct RandomDenseGraphParameters : public RandomGraphParameters {
  /**
   * This gives a natural approx upper limit on the graph density.
   * We keep trying to add a random edge, and once we've failed this many times
   * consecutively, the graph density is probably quite high: probably
   * quite a few edges are in existence, making new edges harder to find.
   * e.g., if this equals 20, then if it ends due to this being exceeded,
   * it means that the chance of a random edge not being present
   * is ~ 1/20, so that the density is (very approximately) ~100%-5% = 95%.
   * (Of course, this is all very approximate; a fun exercise to compute
   * the exact probabilities and expectations...)
   */
  std::size_t max_number_of_consecutive_add_edge_attempts = 10;

  virtual bool add_edges(EdgeSequence& edge_sequence) override;
};

/**
 * First colour the vertices randomly, THEN add the edges.
 * Hence, we always have a known colouring.
 */
struct RandomColouredDenseGraphParameters : public RandomGraphParameters {
  std::size_t max_number_of_consecutive_add_edge_attempts = 10;
  std::size_t max_number_of_colours_to_use = 5;

  virtual bool add_edges(EdgeSequence& edge_sequence) override;
};

/**
 * The class of graphs generated may well be very similar to
 * RandomColouredDenseGraphParameters, with quite different probabilities,
 * but the implementation is a bit different.
 * Generate m sets of k vertices, all with the same colour.
 */
struct RandomColouredKPartiteGraphParameters : public RandomGraphParameters {
  /**
   * Equals the number of colours assigned: every vertex in a single set
   * has the same colour.
   */
  std::size_t number_of_vertex_sets = 1;

  std::size_t number_of_vertices_in_each_set = 2;

  /** Each possible edge will be added independently
   * with a certain probability.
   */
  std::size_t percentage_of_added_edges = 20;

  /**
   * Clears out all existing data and sets the number of vertices
   * appropriately, since existing edges would mess up the known colouring.
   */
  virtual bool add_edges(EdgeSequence& edge_sequence) override;
};

/**
 * A trivial graph: no edges on n vertices!
 * Not random, but it's handy to reuse the interface.
 */
struct EdgelessGraph : public RandomGraphParameters {
  virtual bool add_edges(EdgeSequence& edge_sequence) override;
};

/**
 * A trivial graph: every vertex connected to every other by an edge.
 * Not random, but handy to reuse the interface.
 */
struct CompleteGraph : public RandomGraphParameters {
  virtual bool add_edges(EdgeSequence& edge_sequence) override;
};

}  // namespace tests
}  // namespace graphs
}  // namespace tket
