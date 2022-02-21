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

#include "AdjacencyData.hpp"

#include <algorithm>
#include <set>
#include <sstream>
#include <stdexcept>

#include "Utils/Assert.hpp"

using std::exception;
using std::map;
using std::runtime_error;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

namespace tket {
namespace graphs {

string AdjacencyData::to_string() const {
  map<std::size_t, set<std::size_t>> data_to_display;
  for (std::size_t i = 0; i < m_cleaned_data.size(); ++i) {
    const auto& neighbours = m_cleaned_data[i];
    auto& neighbours_to_display = data_to_display[i];
    for (std::size_t v : neighbours) {
      if (i <= v) {
        neighbours_to_display.insert(v);
      }
    }
    if (neighbours_to_display.empty()) {
      // Don't display it after all
      data_to_display.erase(i);
    }
  }

  stringstream ss;
  ss << "\nThere are " << m_cleaned_data.size()
     << " vertices in total.\nVertex neighbours:\n{";

  for (const auto& entry_to_display : data_to_display) {
    ss << "\n    { " << entry_to_display.first << ", { ";
    for (auto v : entry_to_display.second) {
      ss << v << ", ";
    }
    ss << "} },";
  }
  ss << "\n}\n";
  return ss.str();
}

const set<std::size_t>& AdjacencyData::get_neighbours(
    std::size_t vertex) const {
  // GCOVR_EXCL_START
  TKET_ASSERT(
      vertex < m_cleaned_data.size() ||
      AssertMessage()
          << "AdjacencyData: get_neighbours called with invalid vertex "
          << vertex << "; there are only " << m_cleaned_data.size()
          << " vertices");
  // GCOVR_EXCL_STOP
  return m_cleaned_data[vertex];
}

std::size_t AdjacencyData::get_number_of_vertices() const {
  return m_cleaned_data.size();
}

std::size_t AdjacencyData::get_number_of_edges() const {
  // Each edge i-j is counted twice, i->j and j->i, except for loops.
  std::size_t edges = 0;
  std::size_t loops = 0;

  for (std::size_t ii = 0; ii < m_cleaned_data.size(); ++ii) {
    edges += m_cleaned_data[ii].size();
    loops += m_cleaned_data[ii].count(ii);
  }
  return loops + (edges - loops) / 2;
}

bool AdjacencyData::add_edge(std::size_t i, std::size_t j) {
  const bool exists = edge_exists(i, j);
  if (exists) {
    return false;
  }
  m_cleaned_data[i].insert(j);
  m_cleaned_data[j].insert(i);
  return true;
}

bool AdjacencyData::edge_exists(std::size_t i, std::size_t j) const {
  // GCOVR_EXCL_START
  TKET_ASSERT(
      (i < m_cleaned_data.size() && j < m_cleaned_data.size()) ||
      AssertMessage() << "edge_exists called with vertices " << i << ", " << j
                      << ", but there are only " << m_cleaned_data.size()
                      << " vertices");
  // GCOVR_EXCL_STOP
  return m_cleaned_data[i].count(j) != 0;
}

void AdjacencyData::clear(std::size_t number_of_vertices) {
  m_cleaned_data.resize(number_of_vertices);
  for (auto& entry : m_cleaned_data) {
    entry.clear();
  }
}

AdjacencyData::AdjacencyData(std::size_t number_of_vertices) {
  m_cleaned_data.resize(number_of_vertices);
}

AdjacencyData::AdjacencyData(
    const map<std::size_t, vector<std::size_t>>& raw_data,
    std::size_t number_of_vertices) {
  for (const auto& entry : raw_data) {
    number_of_vertices = std::max(number_of_vertices, entry.first + 1);
    for (std::size_t neighbour : entry.second) {
      number_of_vertices = std::max(number_of_vertices, neighbour + 1);
    }
  }
  m_cleaned_data.resize(number_of_vertices);
  for (const auto& entry : raw_data) {
    for (std::size_t neighbour : entry.second) {
      add_edge(entry.first, neighbour);
    }
  }
}

AdjacencyData::AdjacencyData(
    const vector<vector<std::size_t>>& raw_data, bool allow_loops) {
  m_cleaned_data.resize(raw_data.size());

  for (std::size_t i = 0; i < m_cleaned_data.size(); ++i) {
    for (std::size_t j : raw_data[i]) {
      // GCOVR_EXCL_START
      TKET_ASSERT(
          i != j || allow_loops ||
          AssertMessage() << "Vertex " << i << " out of "
                          << m_cleaned_data.size() << " has a loop.");
      TKET_ASSERT(
          j < m_cleaned_data.size() ||
          AssertMessage() << "Vertex " << i << " has illegal neighbour vertex "
                          << j << ", the size is " << m_cleaned_data.size());
      // GCOVR_EXCL_STOP
      m_cleaned_data[i].insert(j);
      m_cleaned_data[j].insert(i);
    }
  }
}

}  // namespace graphs
}  // namespace tket
