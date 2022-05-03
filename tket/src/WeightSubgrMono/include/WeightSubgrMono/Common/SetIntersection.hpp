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
 * In the worst case, these algorithms use O(1) space and time
 *      O([min size].log [min size].log [max size]),
 * and so asymptotically are no worse than the naive
 * methods; but in practice they are usually faster.
 * They're especially good if the intersection is much smaller
 * than each individual set.
 *
 * It is faster with sorted vectors than std::sets
 * (e.g. std::lower_bound on sorted vectors can take a narrower range,
 * i.e., it takes a start iterator as well as an end iterator,
 * whereas std::set::lower_bound doesn't).
 */

template <class T>
bool disjoint(const std::set<T>& set1, const std::set<T>& set2);

template <class T>
void fill_intersection(
    const std::set<T>& set, const std::vector<T>& sorted_vect,
    std::vector<T>& result);

/** Assume that the ".first" T objects in the vector are distinct,
 * and that the vector is sorted lexicographically w.r.t. these T values.
 * Treat these T values as though they formed a set<T> object,
 * and fill "result" with the intersection.
 */
template <class T, class Number>
void fill_intersection_ignoring_second_elements(
    const std::set<T>& set,
    const std::vector<std::pair<T, Number>>& sorted_vect, std::set<T>& result);

/** Assume that "ExtraData" is some kind of object that contains a T value,
 * plus some kind of extra data, irrelevant (for the purposes of
 * this intersection, at least).
 * Assume that "ExtraData" objects can be ordered, in such a way that
 * distinct T-values within "ExtraData" uniquely determine
 * the order (e.g., std::pair<T, ...> has this property
 * with lexicographic ordering).
 *
 * Assume that we can convert back-and-forth between
 * T values and ExtraData objects (at least for the purposes of
 * this intersection; it doesn't have to be an invertible mapping,
 * but we assume that we can create a "fake" ExtraData object
 * which is good enough to find lower bounds in the vector accurately,
 * as long as the T-elements are DISTINCT.
 * Thus, pass in T -> ExtraData and ExtraData -> T functions.
 */
template <
    class T, class ExtraData, class GetTFromExtraData, class GetExtraDataFromT>
void fill_intersection(
    const std::set<T>& set, const std::vector<ExtraData>& sorted_vect,
    std::set<T>& result_set, GetTFromExtraData get_t,
    GetExtraDataFromT get_extra_data);

// Implementations below:

/** The general algorithm.
 * Assume that "ExtraData" is some kind of object that contains a T value,
 * plus some kind of extra data, irrelevant (for the purposes of
 * this intersection, at least); and that it can be sorted w.r.t. T.
 *
 * Assume that we can convert back-and-forth between
 * T values and ExtraData objects (at least for the purposes of
 * this intersection). Thus, pass in T -> ExtraData and ExtraData -> T
 * functions.
 */
template <
    class T, class ResultInserter, class ExtraData, class GetTFromExtraData,
    class GetExtraDataFromT>
void fill_intersection_using_inserter(
    const std::set<T>& set, const std::vector<ExtraData>& sorted_vect,
    ResultInserter& inserter, GetTFromExtraData get_t,
    GetExtraDataFromT get_extra_data) {
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
      inserter.insert(vect_t_value);
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

namespace internal {

template <class T>
struct SetInserter {
  std::set<T>& result;

  explicit SetInserter(std::set<T>& res) : result(res) {}

  void insert(T value) { result.insert(value); }
};

template <class T>
struct VectorInserter {
  std::vector<T>& result;

  explicit VectorInserter(std::vector<T>& res) : result(res) {}

  void insert(T value) { result.push_back(value); }
};

}  // namespace internal

template <class T, class Number>
void fill_intersection_ignoring_second_elements(
    const std::set<T>& set,
    const std::vector<std::pair<T, Number>>& sorted_vect, std::set<T>& result) {
  result.clear();
  internal::SetInserter<T> inserter(result);

  fill_intersection_using_inserter(
      set, sorted_vect, inserter,
      [](const std::pair<T, Number>& pair) { return pair.first; },
      [](T value) {
        return std::make_pair(value, std::numeric_limits<Number>::min());
      });
}

template <
    class T, class ExtraData, class GetTFromExtraData, class GetExtraDataFromT>
void fill_intersection(
    const std::set<T>& set, const std::vector<ExtraData>& sorted_vect,
    std::set<T>& result, GetTFromExtraData get_t,
    GetExtraDataFromT get_extra_data) {
  result.clear();
  internal::SetInserter<T> inserter(result);

  fill_intersection_using_inserter(
      set, sorted_vect, inserter, get_t, get_extra_data);
}

template <class T>
void fill_intersection(
    const std::set<T>& set, const std::vector<T>& sorted_vect,
    std::vector<T>& result) {
  result.clear();
  internal::VectorInserter<T> inserter(result);

  fill_intersection_using_inserter(
      set, sorted_vect, inserter, [](T value) { return value; },
      [](T value) { return value; });
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
