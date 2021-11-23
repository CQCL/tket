// Copyright 2019-2021 Cambridge Quantum Computing
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
namespace graphs {
namespace tests {

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

/**
 * TODO: move this to a better place, once decided where.
 * A random number generator class.
 * Of course, this is only for random test data generation,
 * definitely NOT suitable for any kind of cryptography!
 * Note that there are no functions involving doubles anywhere!
 * Actually, double calculations can give very slightly different answers
 * across platforms, compilers, compiler optimisation settings;
 * the numerical difference is absolutely negligible,
 * but it's worth being ultra cautious!
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
  std::size_t get_size_t(std::size_t max_value);

  /**
   * Returns a number in the inclusive interval, including the endpoints.
   *
   * @return A size_t from the inclusive range {a, a+1, a+2, ... , b}.
   */
  std::size_t get_size_t(std::size_t min_value, std::size_t max_value);

  /**
   * I believe that the behaviour on the Mersenne twister random engine
   * is guaranteed by the C++ standard, although I'm not 100% sure.
   * The standard specifies 5489u as the default initial seed, so it would
   * be rather pointless to do that if the bits generated
   * were still implementation-dependent.
   */
  void set_seed(std::size_t seed);

  /** Return true p% of the time.
   * (Very quick and dirty, doesn't check for, e.g., 110% effort...)
   * As mentioned above, we deliberately DON'T have a function returning
   * a uniform double. Sticking to integer values is safest.
   */
  bool check_percentage(std::size_t percentage);

  /**
   * Simply shuffle the elements around at random.
   * Approximately uniform over all possible permutations.
   * This is necessary because C++ random_shuffle is
   * implementation-dependent (see above comments).
   */
  template <class T>
  void do_shuffle(std::vector<T>& elements) {
    if (elements.size() < 2) {
      return;
    }
    m_shuffling_data.resize(elements.size());
    for (std::size_t i = 0; i < m_shuffling_data.size(); ++i) {
      m_shuffling_data[i].first = m_engine();
      m_shuffling_data[i].second = i;
    }
    std::sort(
        m_shuffling_data.begin(), m_shuffling_data.end(),
        [](const std::pair<std::uintmax_t, std::size_t>& lhs,
           const std::pair<std::uintmax_t, std::size_t>& rhs) {
          return lhs.first < rhs.first ||
                 (lhs.first == rhs.first && lhs.second < rhs.second);
        });
    // Don't need to make a copy of "elements"! Just do repeated swaps...
    for (std::size_t i = 0; i < m_shuffling_data.size(); ++i) {
      const std::size_t& j = m_shuffling_data[i].second;
      if (i != j) {
        std::swap(elements[i], elements[j]);
      }
    }
  }

  /** Return a random element from the vector. */
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
   */
  template <class T>
  T get_and_remove_element(std::vector<T>& elements) {
    if (elements.empty()) {
      throw std::runtime_error(
          "RNG: get_and_remove_element called on empty vector");
    }
    std::size_t index = get_size_t(elements.size() - 1);
    const T copy = elements[index];
    elements[index] = elements.back();
    elements.pop_back();
    return copy;
  }

  /** Returns the numbers {0,1,2,...,N-1} in some random order. */
  std::vector<std::size_t> get_permutation(std::size_t size);

 private:
  std::mt19937_64 m_engine;

  // Avoids repeated memory reallocation.
  std::vector<std::pair<std::uintmax_t, std::size_t>> m_shuffling_data;
};

}  // namespace tests
}  // namespace graphs
}  // namespace tket
