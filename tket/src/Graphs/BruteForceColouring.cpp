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

#include "BruteForceColouring.hpp"

#include <set>
#include <sstream>
#include <stdexcept>

#include "ColouringPriority.hpp"
#include "Utils/Assert.hpp"

using std::map;
using std::set;
using std::size_t;
using std::vector;

namespace tket {
namespace graphs {

struct BruteForceColouring::Impl {
  struct NodeColouringData {
    vector<size_t> allowed_colours;

    // The index in "allowed_colours", not the actual colour
    size_t current_colour_index = 0;

    bool is_valid_colour() const {
      return current_colour_index < allowed_colours.size();
    }

    size_t get_colour() const { return allowed_colours[current_colour_index]; }
  };

  // This exactly mirrors the nodes in the ColouringPriority object;
  // it is just extra data specifically related to colouring.
  vector<NodeColouringData> colouring_data;

  // KEY: is the vertex, VALUE: the colour
  map<size_t, size_t> colours;

  // Fills in the colour possibilities,
  // possibly increasing suggested_number_of_colours.
  // Returns false if it failed.
  bool initial_colouring_setup(
      const ColouringPriority& priority, size_t& suggested_number_of_colours) {
    const auto& initial_clique = priority.get_initial_clique();

    suggested_number_of_colours =
        std::max(suggested_number_of_colours, initial_clique.size());

    const auto& nodes = priority.get_nodes();

    if (suggested_number_of_colours > nodes.size()) {
      return false;
    }
    colouring_data.resize(nodes.size());

    for (size_t i = 0; i < initial_clique.size(); ++i) {
      colouring_data[i].allowed_colours.assign(1, i);
      colours[nodes[i].vertex] = i;
    }

    // Now, recursively keep increasing the number of colours until we at least
    // have nonempty colour possibility lists. (It may still be impossible, of
    // course, due to the detailed graph structure).

    set<size_t> forbidden_colours;

    for (size_t i = initial_clique.size(); i < nodes.size(); ++i) {
      forbidden_colours.clear();

      for (size_t node_index : nodes[i].earlier_neighbour_node_indices) {
        const auto& earlier_colours =
            colouring_data[node_index].allowed_colours;

        // It is only the initial CLIQUE vertices which have fixed colours;
        // as the number of possible colours increases, EVERY other vertex
        // has more colour possibilities; so we must NOT think that,
        // just because CURRENTLY a vertex has only one colour,
        // that it will ALWAYS be that way!
        if (initial_clique.count(nodes[node_index].vertex) != 0) {
          TKET_ASSERT(earlier_colours.size() == 1);
          forbidden_colours.insert(earlier_colours[0]);
        }
      }
      auto& possible_colours = colouring_data[i].allowed_colours;
      possible_colours.clear();
      for (size_t col = 0; col < suggested_number_of_colours; ++col) {
        if (forbidden_colours.count(col) == 0) {
          possible_colours.push_back(col);
        }
      }
      if (possible_colours.empty()) {
        ++suggested_number_of_colours;
        return initial_colouring_setup(priority, suggested_number_of_colours);
      }
    }
    return true;
  }

  void fill_colour_map(const ColouringPriority& priority) {
    const auto& nodes = priority.get_nodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
      colours[nodes[i].vertex] = colouring_data[i].get_colour();
    }
  }

  bool attempt_brute_force_colouring(const ColouringPriority& priority) {
    for (auto& data : colouring_data) {
      data.current_colour_index = 0;
    }
    const auto& nodes = priority.get_nodes();
    const size_t number_of_nodes = nodes.size();

    for (size_t current_node_index = 0;;) {
      const auto& current_node = nodes[current_node_index];
      auto& current_colouring_node = colouring_data[current_node_index];

      if (current_colouring_node.is_valid_colour()) {
        // We have a candidate colour, now test it for consistency.
        const auto current_col = current_colouring_node.get_colour();
        bool colour_is_impossible = false;

        // TODO: would a set of colours be faster here?
        for (size_t earlier_node_index :
             current_node.earlier_neighbour_node_indices) {
          if (colouring_data[earlier_node_index].get_colour() == current_col) {
            colour_is_impossible = true;
            break;
          }
        }
        if (!colour_is_impossible) {
          // Advance to the next node to colour.
          ++current_node_index;
          if (current_node_index < number_of_nodes) {
            colouring_data[current_node_index].current_colour_index = 0;
            continue;
          }
          // We've hit the end! We are finished.
          return true;
        }
        ++current_colouring_node.current_colour_index;
        continue;
      }

      // We must backtrack.
      if (current_node_index == 0) {
        return false;
      }
      --current_node_index;

      // Advance the colour.
      ++colouring_data[current_node_index].current_colour_index;
    }
    // We cannot actually reach here, the outer loop only returns, never breaks.
  }
};

BruteForceColouring::~BruteForceColouring() {}

BruteForceColouring::BruteForceColouring(
    const ColouringPriority& priority, size_t suggested_number_of_colours)
    : m_pimpl(std::make_unique<BruteForceColouring::Impl>()) {
  const auto number_of_nodes = priority.get_nodes().size();
  if (suggested_number_of_colours >= number_of_nodes) {
    // We've been given permission to use many colours;
    // so just use them all!
    for (size_t i = 0; i < number_of_nodes; ++i) {
      m_pimpl->colours[priority.get_nodes()[i].vertex] = i;
    }
    return;
  }

  const auto initial_suggested_number_of_colours = suggested_number_of_colours;

  try {
    if (!m_pimpl->initial_colouring_setup(
            priority, suggested_number_of_colours)) {
      throw std::runtime_error("initial_colouring_setup failed");
    }

    // From now on, every time we fail to colour with the specified number,
    // we simply have to add the extra colour onto
    // each list of possible colours;
    // no need to redo "initial_colouring_setup".

    for (; suggested_number_of_colours <= number_of_nodes;
         ++suggested_number_of_colours) {
      if (m_pimpl->attempt_brute_force_colouring(priority)) {
        // We've succeeded!
        m_pimpl->fill_colour_map(priority);
        return;
      }
      // It's impossible with this number of colours,
      // so try again with one more.
      // If we were really fancy we might consider
      // a binary search on the number of colours.
      for (size_t i = priority.get_initial_clique().size(); i < number_of_nodes;
           ++i) {
        m_pimpl->colouring_data[i].allowed_colours.push_back(
            suggested_number_of_colours);
      }
    }
    throw std::runtime_error("suggested_number_of_colours hit number_of_nodes");
  } catch (const std::exception& e) {
    // GCOVR_EXCL_START
    TKET_ASSERT(
        AssertMessage() << "initial_suggested_number_of_colours = "
                        << initial_suggested_number_of_colours
                        << ", reached suggested_number_of_colours = "
                        << suggested_number_of_colours << ", had "
                        << number_of_nodes << " nodes. Error: " << e.what()
                        << priority.print_raw_data());
    // GCOVR_EXCL_STOP
  }
}

const map<size_t, size_t>& BruteForceColouring::get_colours() const {
  return m_pimpl->colours;
}

}  // namespace graphs
}  // namespace tket
