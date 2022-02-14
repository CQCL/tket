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

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace tket {
namespace graphs {

// NOTE: GRAPH DATA STRUCTURES/FORMATS:
// in the long term, ideally we should use boost graph data formats
// (i.e., adjacency_list) or some other type for ALL graphs in tket.
//
// Assuming that we do like boost::adjacency_list, which has a lot
// of extra template parameters, it seems there are 4 options.
// In increasing order of work/code complexity:
//
//   (i): do nothing (the short term option!)
//
//  (ii): create boost::... -> {our data structures} conversion routines,
//        and continue to use our routines unchanged.
//        (e.g., brute force colouring is the bottleneck, not conversion).
//
// (iii): create {general boost::...} -> {special restricted boost::...}
//        conversion routines, and change our algorithms to use the
//        restricted boost objects. (With all template parameters fixed,
//        and all irrelevant extra vertex data stripped, so at least
//        we can be sure of efficiency in implementation).
//
//  (iv): change all routines to use general boost::... directly,
//        with templates everywhere, and POSSIBLY unavoidable inefficiencies,
//        because the boost graph stuff is very general.
//
// Another option is to go Java-style and use interfaces/virtual functions
// everywhere instead of templates.

/**
 * Data for an undirected graph. Once stored inside here, the data has
 * automatically been checked and cleaned: the vertices in neighbour lists
 * actually exist, there are no missing vertices, no duplicate edge data, etc.
 * The vertices are {0,1,2,...,v-1}.
 * The number of vertices must be known at the start
 * (or, it can be reset later, but only by clearing all data),
 * so it is not completely dynamic. The constructors throw upon invalid data.
 */
class AdjacencyData {
 public:
  /**
   * Initialise with N vertices, no edges.
   */
  explicit AdjacencyData(std::size_t number_of_vertices = 0);

  /**
   * Construct from a known graph in simple raw format. A sparse format:
   * we don't explicitly need to list vertices if they have no edges.
   * @param raw_data A mapping (vertex i) -> {list of neighbours of i}.
   *      If i->j is an edge, j->i will be automatically deduced also.
   * @param number_of_vertices Optional: if not specified, the number of
   * vertices will be deduced.
   */
  AdjacencyData(
      const std::map<std::size_t, std::vector<std::size_t>>& raw_data,
      std::size_t number_of_vertices = 0);

  /**
   * The vertices are {0,1,2,...,n}.
   * @param raw_data element[i] lists vertices j such that there is an edge
   * i->j. There is no need to list j->i also, it will be automatically deduced.
   * @param allow_loops Optional: default is false. Loops i->i make no sense for
   * colouring.
   */
  AdjacencyData(
      const std::vector<std::vector<std::size_t>>& raw_data,
      // For colouring, loops i -> i are obviously not allowed
      bool allow_loops = false);

  /**
   * Changes the number of vertices and clears all data.
   * @param number_of_vertices The number of vertices, v.
   *   The actual vertices will be {0,1,2,...,v-1}.
   */
  void clear(std::size_t number_of_vertices);

  /** For a given vertex v, return all vertices j such that j-v is an edge.
   */
  const std::set<std::size_t>& get_neighbours(std::size_t vertex) const;

  /** Returns the total number of vertices in the graph.
   */
  std::size_t get_number_of_vertices() const;

  /** Returns the total number of edges in the graph (i->j and j->i counting as
   * one edge).
   */
  std::size_t get_number_of_edges() const;

  /** Returns true if and only if the edge i-j exists. */
  bool edge_exists(std::size_t i, std::size_t j) const;

  /** You must set the number of vertices BEFORE calling this.
   * If the edge i-j does not already exist, adds it, and return true.
   * If it already existed, do nothing and return false.
   * If (i,j) is invalid, throws.
   */
  bool add_edge(std::size_t i, std::size_t j);

  /** to_string is useful for debugging. You can copy the graph data
   * and easily paste it back into C++ code.
   */
  std::string to_string() const;

 private:
  // Element i gives all neighbours j for vertex i, including j<i for speed.
  std::vector<std::set<std::size_t>> m_cleaned_data;
};

}  // namespace graphs
}  // namespace tket
