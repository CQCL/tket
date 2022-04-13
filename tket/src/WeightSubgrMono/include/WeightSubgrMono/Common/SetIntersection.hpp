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
#include <algorithm>
#include <limits>
#include <set>
#include <utility>
#include <vector>

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** This contains functions for intersecting sets, and related things.
 * By interleaving iterators and using lower_bound,
 * the intersection of sets and similar tasks can often be quite a bit faster
 * than the obvious method of just going through
 * the smaller container one-by-one and checking against the larger container.
 *
 * These are usually quite good if the intersection is much smaller
 * than both sets.
 *
 * In the worst case, these algorithms are all O([min size].log [max size]),
 * and so asymptotically no worse than the naive
 * methods; but in practice they should be faster.
 *
 * Also, it is slightly faster to intersect a std::set and a sorted
 * std::vector than two std::sets.
 */

template <class T>
bool disjoint(const std::set<T>& set1, const std::set<T>& set2);

template <class T>
void fill_intersection(
    const std::set<T>& set, const std::vector<T>& sorted_vect,
    std::set<T>& result);

/** Assume that Int is an integer type (signed or unsigned), and that the vector
 * is sorted according to first elements (i.e., the second elements are simply
 * ignored). This is equivalent to lexicograph sorting, of course, since we
 * assume that the T values are distinct.
 */
template <class T, class Int>
void fill_intersection_ignoring_second_elements(
    const std::set<T>& set, const std::vector<std::pair<T, Int>>& sorted_vect,
    std::set<T>& result);

/** Assume that "ExtraData" is some kind of object that contains a T value,
 * plus some kind of extra data, irrelevant (for the purposes of
 * this intersection, at least).
 * Assume that we can convert back-and-forth between
 * T values and ExtraData objects (at least for the purposes of
 * this intersection). Thus, pass in T -> ExtraData and ExtraData -> T
 * functions.
 */
template <
    class T, class ExtraData, class GetTFromExtraData, class GetExtraDataFromT>
void fill_intersection(
    const std::set<T>& set, const std::vector<ExtraData>& sorted_vect,
    std::set<T>& result, GetTFromExtraData get_t,
    GetExtraDataFromT get_extra_data) {
  result.clear();
  if (set.empty() || sorted_vect.empty()) {
    return;
  }
  auto vect_citer = sorted_vect.cbegin();

  for (;;) {
    const T vect_t_value = get_t(*vect_citer);
    auto set_citer = set.lower_bound(vect_t_value);
    if (set_citer == set.cend()) {
      break;
    }
    // We have s >= v, where s in S, v in V.
    if (vect_t_value == *set_citer) {
      // s=v is a common element.
      result.insert(vect_t_value);
      // Advance in S.
      ++set_citer;
      if (set_citer == set.cend()) {
        break;
      }
      // Now, t>s with t in S.
    }
    // If we had s=v above, then now t>v.
    // Otherwise, s>v.
    // In both cases, set_citer now refers to z>v.
    // Get y >= z, with y in V.
    // We can advance v [vector iterator advancement is very fast O(1),
    // unlike O(log N) for std::set].
    ++vect_citer;

    // Note: we're searching over a vector subinterval,
    // so slightly faster than the whole vector.
    vect_citer = std::lower_bound(
        vect_citer, sorted_vect.cend(), get_extra_data(*set_citer));
    if (vect_citer == sorted_vect.cend()) {
      break;
    }
  }
  return;
}

// Further implementations:

// Slightly streamlined, because we return as soon as
// a common element is found.
template <class T>
bool disjoint(const std::set<T>& set1, const std::set<T>& set2) {
  if (set1.empty() || set2.empty()) {
    return true;
  }

  // We always have x1 in S1.
  auto citer1 = set1.cbegin();

  // Break when we run out of elements.
  for (;;) {
    // Get x2 >= x1, where x(j) in S(j).
    auto citer2 = set2.lower_bound(*citer1);
    if (citer2 == set2.cend()) {
      break;
    }
    if (*citer1 == *citer2) {
      // x1=x2; a common element.
      return false;
    }
    // We now know that x2 > x1.
    // Get x1' >= x2.
    citer1 = set1.lower_bound(*citer2);
    if (citer1 == set1.cend()) {
      break;
    }
    if (*citer1 == *citer2) {
      return false;
    }
  }
  return true;
}

template <class T>
void fill_intersection(
    const std::set<T>& set, const std::vector<T>& sorted_vect,
    std::set<T>& result) {
  fill_intersection<T, T>(
      set, sorted_vect, result, [](T value) { return value; },
      [](T value) { return value; });
}

template <class T, class Int>
void fill_intersection_ignoring_second_elements(
    const std::set<T>& set, const std::vector<std::pair<T, Int>>& sorted_vect,
    std::set<T>& result) {
  fill_intersection<T, std::pair<T, Int>>(
      set, sorted_vect, result,
      [](const std::pair<T, Int>& pair) { return pair.first; },
      [](T value) {
        return std::make_pair(value, std::numeric_limits<Int>::min());
      });
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
