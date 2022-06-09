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

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <sstream>

#include "TokenSwapping/VectorListHybrid.hpp"
#include "Utils/RNG.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

typedef VectorListHybrid<unsigned> List;
typedef List::ID ID;

SCENARIO("Reversing a list") {
  RNG rng;
  List list;
  auto copied_elements = list.to_vector();
  vector<unsigned> copied_elements_again;
  REQUIRE(copied_elements.empty());
  for (int count = 0; count < 1000; ++count) {
    const unsigned x = rng.get_size_t(1000);
    switch (x % 7) {
      // Should we delete?
      case 0:
        if (list.size() != 0) {
          const auto id = list.front_id().value();
          list.erase(id);
        }
        break;
      case 1:
        if (list.size() != 0) {
          const auto id = list.back_id().value();
          list.erase(id);
        }
        break;
      case 2:
        list.clear();
        break;
      default:
        break;
    }
    if (x % 2 == 0) {
      list.push_front(x);
    } else {
      list.push_back(x);
    }
    copied_elements = list.to_vector();
    list.reverse();
    copied_elements_again = list.to_vector();
    std::reverse(copied_elements.begin(), copied_elements.end());
    REQUIRE(copied_elements == copied_elements_again);
  }
}

// Write the contents to a string for testing, possibly including IDs.
static std::string repr(const List& list, bool include_ids) {
  std::stringstream ss;
  ss << "[size " << list.size() << ": ";
  for (auto id_opt = list.front_id(); id_opt;) {
    const auto id = id_opt.value();
    id_opt = list.next(id);
    ss << list.at(id) << " ";
  }
  if (include_ids) {
    ss << "; ids: ";
    for (auto id_opt = list.front_id(); id_opt;) {
      const auto id = id_opt.value();
      id_opt = list.next(id);
      ss << id << " ";
    }
  }
  ss << "]";
  return ss.str();
}

// In "operations", a positive number p means go to position p % size() in the
// list, and insert a number there. A negative number n means do the same thing
// with abs(n) % size(), but erase instead of insert. Returns a string
// representing the elements which were erased/inserted, again using negative
// numbers to denote erasure. Does NOT give the IDs.
static std::string perform_operation(
    const vector<int>& operations, List& list, unsigned& next_element) {
  std::stringstream ss;
  ss << "[";
  for (int position_code : operations) {
    REQUIRE(position_code != 0);
    const auto size = list.size();
    if (size == 0) {
      if (position_code > 0) {
        list.push_back(next_element);
        ss << "new: " << next_element << " ";
        next_element += 100;
        continue;
      }
      // Cannot erase from an empty list!
      ss << "; ";
      continue;
    }
    // It's nonempty.
    unsigned position = std::abs(position_code);
    position %= size;
    ID id = list.front_id().value();
    for (unsigned nn = 0; nn < position; ++nn) {
      const auto next_id = list.next(id);
      REQUIRE(next_id);
      id = next_id.value();
    }
    ss << "at " << position << ": ";
    if (position_code > 0) {
      ss << next_element << " ";
      const ID new_id = list.insert_after(id);
      list.at(new_id) = next_element;
      next_element += 100;
      continue;
    }
    ss << "-" << list.at(id) << " ";
    list.erase(id);
  }
  ss << "]";
  return ss.str();
}

namespace {
struct Result {
  std::string initial_op_str;
  std::string list_str_after_one_op;
  std::string list_str_after_one_op_without_ids;
  std::string op_str_after_two_ops;
  std::string list_str_after_two_ops;
  std::string list_str_after_two_ops_without_ids;

  Result(const vector<int>& operations, List& list, unsigned& next_element)
      : initial_op_str(perform_operation(operations, list, next_element)),
        list_str_after_one_op(repr(list, true)),
        list_str_after_one_op_without_ids(repr(list, false)),
        op_str_after_two_ops(perform_operation(operations, list, next_element)),
        list_str_after_two_ops(repr(list, true)),
        list_str_after_two_ops_without_ids(repr(list, false)) {}

  void check_equal_contents_without_ids(const Result& other) const {
    CHECK(initial_op_str == other.initial_op_str);
    CHECK(
        list_str_after_one_op_without_ids ==
        other.list_str_after_one_op_without_ids);
    CHECK(op_str_after_two_ops == other.op_str_after_two_ops);
    CHECK(
        list_str_after_two_ops_without_ids ==
        other.list_str_after_two_ops_without_ids);
  }

  void check_equal_id_data(const Result& other) const {
    CHECK(list_str_after_one_op == other.list_str_after_one_op);
    CHECK(list_str_after_two_ops == other.list_str_after_two_ops);
  }

  void check_different_id_data(const Result& other) const {
    CHECK(list_str_after_one_op != other.list_str_after_one_op);
    CHECK(list_str_after_two_ops != other.list_str_after_two_ops);
  }
};
}  // namespace

// We want to test that lists have equal or different contents,
// with/without clear/fast_clear, etc.
// The same sequences of logical operations
// (erase, insert, etc.) applied to a new list or a fast_cleared list might NOT
// preserve IDs, but should preserve the contents. With clear(), it should ALSO
// preserve IDs.
SCENARIO("Inserting, erasing, clearing tests") {
  // These are just some random numbers.
  const vector<int> operations{-10, -4, 1,  3,  -8, 2,  -2, -3, -5, -9,
                               -6,  -2, -7, 2,  5,  -8, 6,  -4, 10, 7,
                               -10, -1, 5,  6,  9,  1,  4,  -7, -1, 4,
                               8,   -9, 8,  -3, -5, -6, 9,  3,  7,  10};

  List list;
  unsigned next_element = 999;
  const Result result_with_new_object(operations, list, next_element);

  // Also test clearing empty objects.
  {
    // bits 00 mean do nothing, 01 means clear, 11 means fast clear.
    const vector<unsigned> clear_options{
        0,    // nothing,
        0x5,  // clear, clear,
        0x7,  // fast clear, clear,
        0xD,  // clear, fast clear,
        0xF,  // fast clear, fast clear
        0x15  // clear, clear, clear
    };
    for (unsigned option : clear_options) {
      List empty_list;
      unsigned copy = option;
      while (copy != 0) {
        const unsigned code = copy & 0x3;
        copy >>= 2;
        switch (code) {
          case 1:
            empty_list.clear();
            break;
          case 3:
            empty_list.fast_clear();
            break;
          default: {
            REQUIRE(false);
          }
        }
      }
      next_element = 999;
      const Result result_with_empty_list(operations, empty_list, next_element);
      result_with_empty_list.check_equal_contents_without_ids(
          result_with_new_object);
      result_with_empty_list.check_equal_id_data(result_with_new_object);
    }
  }
  // Now repeat the operations.
  list.clear();
  {
    INFO("second time, cleared list");
    next_element = 999;
    const Result result_with_cleared_object(operations, list, next_element);
    result_with_cleared_object.check_equal_contents_without_ids(
        result_with_new_object);
    result_with_cleared_object.check_equal_id_data(result_with_new_object);
  }
  list.fast_clear();
  {
    INFO("third time, fast cleared list");
    next_element = 999;
    const Result result_with_cleared_object(operations, list, next_element);
    result_with_cleared_object.check_equal_contents_without_ids(
        result_with_new_object);
    result_with_cleared_object.check_different_id_data(result_with_new_object);
  }
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
