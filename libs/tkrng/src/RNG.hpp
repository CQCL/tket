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
#include <random>
#include <stdexcept>
#include <vector>

namespace tket {

// Something like this is needed for proper random test data generation
// if you want to be platform-independent, as the C++ standard is stupid.
// (A major scandal, in my opinion).
// The random engines are mostly guaranteed by the standard,
// but the DISTRIBUTIONS, e.g. uniform_distribution, are NOT
// (i.e., the actual algorithm used to convert a string of bits to a number
// in the range {0,1,2,...,N} is not specified at all by the C++ standard).
// Thus, we are NOT guaranteed to get the same results, even with the same
// (1) engine; (2) initial seed; (3) distribution,
// if we use different platforms (or even different compilers
// on the SAME platform), or even different compiler VERSIONS!!!
//
// The C++ standard as far as I know does not specify ANY distribution
// implementations, not even optionally, so you HAVE to do this yourself,
// even for something as simple as a uniform distribution.
// The same applies to, e.g., std::random_shuffle.

/** A random number generator class.
 * Of course, this is only for random test data generation,
 * definitely NOT suitable for any kind of cryptography!
 * Note that there are no functions involving doubles anywhere!
 * Actually, double calculations can give very slightly different answers
 * across platforms, compilers, compiler optimisation settings;
 * the numerical difference is absolutely negligible,
 * but it's worth being ultra cautious! (And it's much easier for testing
 * to get IDENTICAL results across platforms).
 */
class RNG {
 public:
  /**
   * Return a random integer from 0 to N, inclusive.
   * Approximately uniform, if max_value is much less than
   * the max possible value that can be returned.
   * N << sqrt(max uint64) ~ 2^32 ~ 4e9 will work well.
   * See the comments in the cpp file implementation for more detail.
   *
   * @param max_value The value N which is the (inclusive) maximum value
   * which can be returned.
   * @return A size_t from the inclusive range {0,1,2,...,N}.
   */
  size_t get_size_t(size_t max_value);

  /**
   * Returns a number in the inclusive interval, including the endpoints.
   * @param min_value The smallest value (inclusive) that can be returned.
   * @param max_value The largest value (inclusive) that can be returned.
   * @return A size_t from the inclusive range {a, a+1, a+2, ... , b}.
   */
  size_t get_size_t(size_t min_value, size_t max_value);

  /**
   * The behaviour of the RAW BITS of the Mersenne twister random engine
   * is guaranteed by the C++ standard.
   * The standard specifies 5489u as the default initial seed.
   * @param seed A seed value, to alter the RNG state.
   *    By default, uses the value specified by the standard.
   */
  void set_seed(size_t seed = 5489);

  /** Return true p% of the time.
   * (Very quick and dirty, doesn't check for, e.g., 110% effort...)
   * As mentioned above, we deliberately DON'T have a function returning
   * a uniform double. Sticking to integer values is safest.
   * @param percentage The probability of returning true, expressed as
   *  a percentage.
   * @return A random bool, returns true with specified probability.
   */
  bool check_percentage(size_t percentage);

  /**
   * Simply shuffle the elements around at random.
   * Approximately uniform "in practice" over all possible permutations.
   * (Although of course, strictly speaking very far from uniform for larger
   * vectors. The number of possible permutations grows very rapidly
   * and quickly becomes larger than the total number of distinct states
   * any fixed engine can take, no matter which is used. Thus, for larger
   * vectors, only a small proportion of permutations are actually possible).
   * This is necessary because C++ random_shuffle is
   * implementation-dependent (see above comments).
   * @param elements The vector to be shuffled randomly.
   */
  template <class T>
  void do_shuffle(std::vector<T>& elements) {
    if (elements.size() < 2) {
      return;
    }
    m_shuffling_data.resize(elements.size());
    for (size_t i = 0; i < m_shuffling_data.size(); ++i) {
      m_shuffling_data[i].first = m_engine();
      // Tricky subtle point: without this extra entry to break ties,
      // std::sort could give DIFFERENT results across platforms and compilers,
      // if the object T allows unequal elements comparing equal.
      m_shuffling_data[i].second = i;
    }
    std::sort(
        m_shuffling_data.begin(), m_shuffling_data.end(),
        [](const std::pair<std::uintmax_t, size_t>& lhs,
           const std::pair<std::uintmax_t, size_t>& rhs) {
          return lhs.first < rhs.first ||
                 (lhs.first == rhs.first && lhs.second < rhs.second);
        });
    // Don't need to make a copy of "elements"! Just do repeated swaps...
    for (size_t i = 0; i < m_shuffling_data.size(); ++i) {
      const size_t& j = m_shuffling_data[i].second;
      if (i != j) {
        std::swap(elements[i], elements[j]);
      }
    }
  }

  /** Return a random element from the vector.
   *  @param elements The vector to be sampled from.
   *  @return A reference to a random element, approximately uniform.
   */
  template <class T>
  const T& get_element(const std::vector<T>& elements) {
    if (elements.empty()) {
      throw std::runtime_error("RNG: get_element called on empty vector");
    }
    return elements[get_size_t(elements.size() - 1)];
  }

  /**
   * Pick out a random element from the vector, copy and return it,
   * but also remove that element from the vector (swapping with
   * the back for efficiency, i.e. the ordering changes).
   * @param elements The vector to be sampled from.
   *  Decreases size by one each time.
   *  Time O(1) because the ordering is allowed to change.
   * @return A copy of the removed element.
   */
  template <class T>
  T get_and_remove_element(std::vector<T>& elements) {
    if (elements.empty()) {
      throw std::runtime_error(
          "RNG: get_and_remove_element called on empty vector");
    }
    size_t index = get_size_t(elements.size() - 1);
    const T copy = elements[index];
    elements[index] = elements.back();
    elements.pop_back();
    return copy;
  }

  /** Returns the numbers {0,1,2,...,N-1} in some random order.
   *  @param size The size of the returned vector.
   *  @return An interval of nonnegative numbers, starting at zero,
   *    but rearranged randomly.
   */
  std::vector<size_t> get_permutation(size_t size);

 private:
  std::mt19937_64 m_engine;

  // Avoids repeated memory reallocation.
  std::vector<std::pair<std::uintmax_t, size_t>> m_shuffling_data;
};

}  // namespace tket
