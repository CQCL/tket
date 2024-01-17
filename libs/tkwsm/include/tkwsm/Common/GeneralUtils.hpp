// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <algorithm>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "SpecialExceptions.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// GCOVR_EXCL_START
// A bug in test coverage? The following functions DEFINITELY are
// used many times, but not according to test coverage...

/** Sets the numeric variable to its maximum possible value.
 * This helps avoid mixing types accidentally, e.g. writing
 * x = std::numeric_limits<long>::max() where x is actually an int.
 * @param val the numeric variable, to be set to max()
 */
template <class T>
void set_maximum(T& val) {
  val = std::numeric_limits<T>::max();
}

/** Simply check if the variable does have its maximum possible value.
 * @param val the numeric variable
 * @return True if it equals its maximum possible value
 */
template <class T>
bool is_maximum(const T& val) {
  return val == std::numeric_limits<T>::max();
}
// GCOVR_EXCL_STOP

/** Handy for testing; a string represention of a std container.
 * @param elems A container of elements.
 * @param max_elems_to_print The maximum number to print, before terminating
 * early.
 * @return a string representation of the elements.
 */
template <class Container>
std::string str(const Container& elems, std::size_t max_elems_to_print = 10) {
  std::stringstream ss;
  if (elems.size() > 3) {
    ss << elems.size() << " elems: ";
  }
  ss << "[ ";
  std::size_t number_printed = 0;
  for (const auto& elem : elems) {
    ss << elem << " ";
    ++number_printed;
    if (number_printed == elems.size()) {
      break;
    }
    if (number_printed >= max_elems_to_print) {
      ss << "...";
      break;
    }
  }
  ss << "]";
  return ss.str();
}

// GCOVR_EXCL_START
// A bug in test coverage?
/** Handy for testing.
 * @param elems An ordinary std::vector of comparable elements.
 * @return True if they are in increasing order, with all values distinct.
 */
template <class T>
bool is_sorted_and_unique(const std::vector<T>& elems) {
  return std::is_sorted(elems.cbegin(), elems.cend()) &&
         (std::adjacent_find(elems.cbegin(), elems.cend()) == elems.cend());
}
// GCOVR_EXCL_STOP

/** Checks if the map has this key.
 * @param map A std::map
 * @param key A key to check
 * @return The value in the map corresponding to the key if it exists, or null
 * optional object if the key does not exist.
 */
template <class Key, class Value>
std::optional<Value> get_optional_value(
    const std::map<Key, Value>& map, const Key& key) {
  const auto citer = map.find(key);
  if (citer == map.cend()) {
    return {};
  }
  return citer->second;
}

// GCOVR_EXCL_START
// There MUST be a bug in test coverage, surely?
// Below, "get_sum_or_throw" definitely calls "get_checked_sum" in all cases.
// Yet the SAME coverage report (for the same file!)
// reported the first as covered (including the actual lines calling
// the second),
// but the second as uncovered (every line!)

/** For an unsigned integer type, returns x+y if the value is correct,
 * or null if overflow would occur (so the actual value of x+y does not fit in
 * the type).
 * @param x First number
 * @param y Second number
 * @return x+y, if the value fits in the unsigned integer type; otherwise null.
 */
template <typename UINT>
std::optional<UINT> get_checked_sum(UINT x, UINT y) {
  if (x == 0) {
    return y;
  }
  if (y == 0) {
    return x;
  }
  const UINT sum = x + y;
  // x,y >= 1, and if   n = x add y  (the C add operation)
  // then  n = x+y (mod M)  (as ordinary integers),
  // with M being the max possible value.
  // If overflow didn't occur, then  n=x+y is the correct value,
  // and so n > x (of course!)
  //
  // If overflow occurs, then
  //  M <= x+y <= 2M  (as ordinary integers).
  // So, if   x+y = 2M   then x=y=M, and then
  //   n = 0 < x.
  // Otherwise, x+y < 2M, so that  0 <= x+y-M < M,
  // and thus   n = x+y-M = x - (M-y)  (as ordinary integers).
  // Clearly, now  n <= x (and also n <= y).
  if (sum <= x) {
    return {};
  }
  return sum;
}

/** For an unsigned integer type, returns x*y
 * if the value is small enough to fit inside a UINT.
 * Otherwise, returns null.
 * @param x First number
 * @param y Second number
 * @return x*y, if the value fits in the unsigned integer type; otherwise null.
 */
template <typename UINT>
std::optional<UINT> get_checked_product(UINT x, UINT y) {
  switch (x) {
    case 0:
      return 0;
    case 1:
      return y;
    default:
      break;
  }
  switch (y) {
    case 0:
      return 0;
    case 1:
      return x;
    default:
      break;
  }
  if (y > std::numeric_limits<UINT>::max() / x) {
    // By definition of unsigned integer type division for x > 0,
    // M div x = n
    // (where div means the C++ integer division operation)
    // if and only if   n <= M/x < n+1  (as real numbers).
    // Therefore nx <= M,  (n+1)x > M  (as real numbers),
    // so  yx <= M  (as real numbers)  <==>  y <= n.
    //
    // (OK, int division can be a lot slower than multiplication.
    // But I think it's unavoidable; we cannot just inspect
    // x MULT y directly, as for addition.
    // E.g., if x ~ sqrt(M), y = 1.5x,  x even, then
    //   xy = 1.5x^2 ~ 1.5 M will be converted to ~0.5 M,
    // which is far bigger than both x and y).
    return {};
  }
  return x * y;
}
// GCOVR_EXCL_STOP

/** Throws the special IntegerOverflow exception if the values are too big.
 * @param x First number
 * @param y Second number
 * @return x+y, if the value fits in the unsigned integer type; otherwise throws
 * an IntegerOverflow exception.
 */
template <typename UINT>
UINT get_sum_or_throw(UINT x, UINT y) {
  const auto sum_opt = get_checked_sum(x, y);
  if (sum_opt) {
    return sum_opt.value();
  }
  std::stringstream ss;
  ss << "(" << x << " + " << y << ")";
  throw IntegerOverflow(ss.str());
}

/** Throws the special IntegerOverflow exception if the values are too big.
 * @param x First number
 * @param y Second number
 * @return x*y, if the value fits in the unsigned integer type; otherwise throws
 * an IntegerOverflow exception.
 */
template <typename UINT>
UINT get_product_or_throw(UINT x, UINT y) {
  const auto product_opt = get_checked_product(x, y);
  if (product_opt) {
    return product_opt.value();
  }
  std::stringstream ss;
  ss << "(" << x << " * " << y << ")";
  throw IntegerOverflow(ss.str());
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
