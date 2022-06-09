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
#include <list>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "TokenSwapping/VectorListHybridSkeleton.hpp"
#include "Utils/RNG.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

// A slower implementation of VectorListHybridSkeleton
// using linked lists
struct VLHS_tester_reimplementation {
  // Each node will contain the index it was given.
  mutable std::list<size_t> data;

  void clear() { data.clear(); }
  size_t size() const { return data.size(); }
  size_t front_index() const { return data.front(); }
  size_t back_index() const { return data.back(); }

  std::list<size_t>::iterator find(size_t index) {
    for (auto iter = data.begin(); iter != data.end(); ++iter) {
      if (*iter == index) {
        return iter;
      }
    }
    throw std::runtime_error(
        std::string("index ") + std::to_string(index) + " not found");
  }

  std::list<size_t>::const_iterator find(size_t index) const {
    for (auto citer = data.cbegin(); citer != data.cend(); ++citer) {
      if (*citer == index) {
        return citer;
      }
    }
    throw std::runtime_error(
        std::string("index ") + std::to_string(index) + " not found");
  }

  std::optional<size_t> next(size_t index) const {
    auto citer = find(index);
    ++citer;
    if (citer == data.cend()) {
      return {};
    }
    return *citer;
  }

  std::optional<size_t> previous(size_t index) const {
    auto citer = find(index);
    --citer;
    if (citer != data.cend()) {
      return {};
    }
    return *citer;
  }

  void erase(size_t index) {
    auto iter = find(index);
    data.erase(iter);
  }

  void insert_for_empty_list(size_t new_index) {
    REQUIRE(data.empty());
    data.push_front(new_index);
  }

  void insert_after(size_t index, size_t new_index) {
    auto iter = find(index);
    // We can only insert BEFORE an iter with STL
    ++iter;
    if (iter == data.end()) {
      // We were at the back.
      data.push_back(new_index);
      return;
    }
    // We're now after the node, we insert before
    data.insert(iter, new_index);
  }

  void insert_before(size_t index, size_t new_index) {
    auto iter = find(index);
    data.insert(iter, new_index);
  }
};

// Keep track of which indices have currently not yet been erased
struct ValidIndices {
  std::set<size_t> indices;

  bool contains(size_t index) const { return indices.count(index) != 0; }
  void check_and_insert_new_index(size_t index) {
    REQUIRE(index != VectorListHybridSkeleton::get_invalid_index());
    REQUIRE(indices.count(index) == 0);
    indices.insert(index);
  }

  void check_and_erase_index(size_t index) {
    REQUIRE(indices.count(index) != 0);
    indices.erase(index);
  }

  size_t get_index(RNG& rng) const {
    REQUIRE(!indices.empty());
    auto citer = indices.cbegin();
    for (size_t ii = rng.get_size_t(indices.size() - 1); ii != 0; --ii) {
      ++citer;
    }
    return *citer;
  }
};

void require_equal_indices(
    size_t index, const std::optional<size_t>& index_opt) {
  if (index == VectorListHybridSkeleton::get_invalid_index()) {
    REQUIRE(!index_opt);
    return;
  }
  REQUIRE(index_opt);
  REQUIRE(index_opt.value() == index);
}

bool are_equal(
    const VectorListHybridSkeleton& vlhs,
    const VLHS_tester_reimplementation& tester,
    const ValidIndices& valid_indices) {
  if (vlhs.size() != tester.size()) {
    return false;
  }
  if (vlhs.size() == 0) {
    return true;
  }
  auto citer = tester.data.cbegin();
  for (auto index = vlhs.front_index();
       index != VectorListHybridSkeleton::get_invalid_index();
       index = vlhs.next(index)) {
    if (*citer != index) {
      return false;
    }
    REQUIRE(valid_indices.contains(index));
    ++citer;
  }
  REQUIRE(citer == tester.data.cend());
  REQUIRE(*tester.data.cbegin() == vlhs.front_index());
  REQUIRE(*tester.data.crbegin() == vlhs.back_index());
  return true;
}

SCENARIO("Random operations preserve VLHS") {
  RNG rng;
  VLHS_tester_reimplementation tester;
  VectorListHybridSkeleton vlhs;
  ValidIndices valid_indices;
  REQUIRE(are_equal(vlhs, tester, valid_indices));

  for (int op_counter = 0; op_counter < 10000; ++op_counter) {
    INFO("counter=" << op_counter);
    if (op_counter + 1 % 100 == 0) {
      vlhs.clear();
      tester.clear();
      valid_indices.indices.clear();
    }
    bool should_insert = rng.check_percentage(50);
    if (valid_indices.indices.empty()) {
      should_insert = true;
    }
    if (valid_indices.indices.size() > 10) {
      should_insert = false;
    }
    if (should_insert) {
      if (valid_indices.indices.empty()) {
        vlhs.insert_for_empty_list();
        const auto new_index = vlhs.front_index();
        REQUIRE(new_index == vlhs.back_index());
        tester.insert_for_empty_list(new_index);
        valid_indices.check_and_insert_new_index(new_index);
      } else {
        const auto index = valid_indices.get_index(rng);
        const bool insert_after = rng.check_percentage(50);

        if (insert_after) {
          vlhs.insert_after(index);
          const auto new_index = vlhs.next(index);
          tester.insert_after(index, new_index);
          valid_indices.check_and_insert_new_index(new_index);
        } else {
          vlhs.insert_before(index);
          const auto new_index = vlhs.previous(index);
          tester.insert_before(index, new_index);
          valid_indices.check_and_insert_new_index(new_index);
        }
      }
    } else {
      // We erase instead.
      const auto index = valid_indices.get_index(rng);
      vlhs.erase(index);
      tester.erase(index);
      valid_indices.check_and_erase_index(index);
    }
    REQUIRE(are_equal(vlhs, tester, valid_indices));
  }
}

static std::string get_fixed_ops_str(bool do_fast_clear) {
  std::stringstream ss;
  VectorListHybridSkeleton vlhs;
  ss << vlhs.debug_str();
  vlhs.insert_for_empty_list();
  ss << "\nInsert: " << vlhs.debug_str();
  vlhs.insert_after(vlhs.front_index());
  ss << "\nInsert after front: " << vlhs.debug_str();
  const auto id = vlhs.front_index();
  vlhs.insert_before(id);
  ss << "\nInsert before front: " << vlhs.debug_str();
  vlhs.insert_after(id);
  ss << "\nInsert after " << id << ": " << vlhs.debug_str();
  vlhs.erase(3);
  ss << "\nErase 3: " << vlhs.debug_str();
  if (do_fast_clear) {
    vlhs.fast_clear();
    ss << "\nFast clear: " << vlhs.debug_str();
  } else {
    vlhs.clear();
    ss << "\nClear: " << vlhs.debug_str();
  }
  vlhs.insert_for_empty_list();
  ss << "\nInsert: " << vlhs.debug_str();
  return ss.str();
}

SCENARIO("Some fixed ops") {
  // The only difference should be in the internal link values.
  const std::string common_prefix{
      "VLHS: size 0, front NULL back NULL, del.front NULL\n"
      "Active links: forward []\n"
      "Backward ()\n"
      "Del.links: {}\n"
      "Insert: VLHS: size 1, front 0 back 0, del.front NULL\n"
      "Active links: forward [0->]\n"
      "Backward (0->)\n"
      "Del.links: {}\n"
      "Insert after front: VLHS: size 2, front 0 back 1, del.front NULL\n"
      "Active links: forward [0->1->]\n"
      "Backward (1->0->)\n"
      "Del.links: {}\n"
      "Insert before front: VLHS: size 3, front 2 back 1, del.front NULL\n"
      "Active links: forward [2->0->1->]\n"
      "Backward (1->0->2->)\n"
      "Del.links: {}\n"
      "Insert after 0: VLHS: size 4, front 2 back 1, del.front NULL\n"
      "Active links: forward [2->0->3->1->]\n"
      "Backward (1->3->0->2->)\n"
      "Del.links: {}\n"
      "Erase 3: VLHS: size 3, front 2 back 1, del.front 3\n"
      "Active links: forward [2->0->1->]\n"
      "Backward (1->0->2->)\n"
      "Del.links: {3->}\n"};
  const std::string fast_clear_suffix{
      "Fast clear: VLHS: size 0, front NULL back NULL, del.front 2\n"
      "Active links: forward []\n"
      "Backward ()\n"
      "Del.links: {2->0->1->3->}\n"
      "Insert: VLHS: size 1, front 2 back 2, del.front 0\n"
      "Active links: forward [2->]\n"
      "Backward (2->)\n"
      "Del.links: {0->1->3->}"};
  const std::string clear_suffix{
      "Clear: VLHS: size 0, front NULL back NULL, del.front 0\n"
      "Active links: forward []\n"
      "Backward ()\n"
      "Del.links: {0->1->2->3->}\n"
      "Insert: VLHS: size 1, front 0 back 0, del.front 1\n"
      "Active links: forward [0->]\n"
      "Backward (0->)\n"
      "Del.links: {1->2->3->}"};
  const auto fast_clear_str = get_fixed_ops_str(true);
  CHECK(fast_clear_str == common_prefix + fast_clear_suffix);

  const auto clear_str = get_fixed_ops_str(false);
  CHECK(clear_str == common_prefix + clear_suffix);
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
