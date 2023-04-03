// Copyright 2019-2023 Cambridge Quantum Computing
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
#include <cstdint>
#include <string>

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** For representing fractions n.2^p for n>=0, with n,p integral.
 * Used for comparing ratios of integers,
 * without overflow in all reasonable cases.
 * We do this because we use integral weights and want
 * to take ratios, but don't want to use doubles and
 * need to avoid overflow.
 *
 * This could be regarded as a very partial software implementation
 * of wider-range doubles, using only integer operations
 * and allowing only multiplication and positive numbers
 * (no addition, subtraction, division).
 * Close to best possible in terms of accuracy (discards
 * the least signficant bits first, and retains almost the
 * maximum number of bits).
 */
class DyadicFraction {
 public:
  typedef std::uint64_t UInt;

  /** Store the value x. */
  explicit DyadicFraction(UInt x = 0);

  /** Multiply our value by x.
   * @param x An integer.
   * @return this object, for chaining, after multiplying by x.
   */
  DyadicFraction& mult(UInt x);

  /** Multiply our value by the value in the other fraction.
   * Note that some bits of accuracy may be lost
   * (just as for doubles); however, it will shift bits sensibly
   * to try to minimise the accuracy loss.
   * @param other The other number.
   * @return this object, for chaining, after multiplying by the other number.
   */
  DyadicFraction& mult(const DyadicFraction& other);

  /** Multiply by n/K = n/1024, where K=1024.
   * @param n A positive integer.
   * @return this object, for chaining, after multiplying by n/1024.
   */
  DyadicFraction& mult_n_over_k(UInt n);

  /** Only for testing. Note that double values are NOT exactly portable;
   * different compilers, platforms, optimisation/fast-math
   * settings can produce very slightly different results!
   * @return The approximate value of this number, as an ordinary double.
   */
  double get_double() const;

  /** Only for testing.
   * @return The approximate value of log(this number), as an ordinary double.
   * This may be accurate even when the double returned by get_double() is
   * inaccurate, due to being too small or large.
   */
  double get_log() const;

  /** Only for testing.
   * @return a string representation of this object. However, it is portable and
   * exact.
   */
  std::string str() const;

  /** Like doubles, we suffer from "roundoff error" (loss of bits),
   * which we try to minimise; to remind us,
   * we deliberately do NOT include <= or == operators.
   * @param other The other number
   * @return True if this number is less than the other number
   */
  bool operator<(const DyadicFraction& other) const;

 private:
  /** The value n, such that x = n.2^p. */
  UInt m_value;

  /** The value p, such that x = n.2^p.
   * The point is that p WILL NOT overflow, until we start getting to
   * ridiculously small/big numbers like 2^{4 billion},
   * which will never happen in the intended applications.
   */
  int m_exponent;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
