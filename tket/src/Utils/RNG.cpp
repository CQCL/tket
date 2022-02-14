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

#include "Utils/RNG.hpp"

using std::vector;

namespace tket {

size_t RNG::get_size_t(size_t max_value) {
  if (max_value == 0) {
    return 0;
  }
  // Raw data; now must convert to a value to return!
  const std::uint64_t random_int = m_engine();

  if (max_value > m_engine.max() / 4) {
    // If choosing such a large potential number of values,
    // the bias will unavoidably be very bad,
    // if only generating a single random int.
    // Surely no deterministic function
    //    f : {0,1,...,N} -> {0,1,...,M}
    // can be close to giving a uniform distribution,
    // if N != M are both large and nearly equal.
    // (Should be a theorem in here somewhere!)
    if (max_value >= m_engine.max()) {
      // Care! Maybe max_value+1 == 0 by wraparound,
      // so we cannot do division by max_value+1 !
      return random_int;
    }
    return random_int % (max_value + 1);
  }

  // NOW we know that max_value+1 won't overflow.

  // Mathematical note on the below: let:
  //    m = maximum possible value of "random_int"
  //    w = interval_width
  //    v = max possible value to return.
  //
  // Thus, random_int could be one of {0,1,2,...,m},
  // and we must return one of {0,1,2,...,v}.
  //
  // With int arithmetic, we get w = int((m+1)/(v+1)).
  //
  // e.g., if m=5, v=2, then w = int(6/3) = 2,
  // the possible random_int values are {0,1,2,3,4,5},
  // and this is partitioned into 3 sets:
  // {0,1},  {2,3},  {4,5}.
  //
  // [Since, with int arithmetic,
  //     int(0/2) = int(1/2) = 0,  int(2/2) = int(3/2) = 1, ...]
  //
  // Because these sets have equal size 2, each of the values 0,1,2
  // has equal probability 1/3 of being returned.
  // BUT, what if (m+1)/(v+1) is not an integer?
  //
  // e.g., m=5, v=3.
  // Now, we must partition the set {0,1,2,3,4,5} into 4 subsets.
  // With the below algorithm, w=int((5+1)/(3+1)) = 1, so the partition is
  // {0}, {1}, {2}, {3,4,5}.
  // Notice that 0,1,2 have probabilities 1/6 of being returned,
  // but v=3 has probability 3/6 of being returned, quite a large bias.
  //
  // How bad can it be? In general:
  //
  //       (m+1)/(v+1) - 1 < w <= (m+1)/(v+1).
  // Thus
  //       m-v+1 <= w(v+1) <= m+1.
  //
  // Now, the random_int sets causing the values 0,1,...,v to be returned are
  //
  //  { 0, 1, ...,    w-1}  -->  returns 0
  //  { w, w+1, ...,  2w-1}  -->  returns 1
  //  {2w, 2w+1, ..., 3w-1}  -->  returns 2
  //      ....
  //  {vw, vw+1, ..., (v+1)w - 1, ... , m } -->   returns v
  //
  // Notice that all sets except the last have size w.
  // The last set has size  m-vw+1. So, the final value v has
  // more ways than the other values 0,1,... to be returned, by a factor of
  //
  //    U = (m-vw+1)/w = (m+1)/w - v.
  //
  // U is the "bias factor" which we want to be as close to 1 as possible.
  // Always,  U >= (m+1)/[(m+1)/(v+1)] - v = v+1-v = 1,
  // as we already know. Also,
  //
  //    U <= (m+1)(v+1)/(m-v+1) - v.
  //
  // Let's assume that v << m.
  // Then we can expand with a geometric series:
  //  (m+1)(v+1)/(m-v+1) = (v+1).[1-v/(m+1)]^{-1}
  //                     = (v+1).[1 + v/(m+1) + A]
  //                     = v+1 + v(v+1)/(m+1) + (v+1)A,
  //
  // where A ~ (v/m)^2,  with ~ here meaning
  // "roughly equal size, up to constant factors".
  // Thus, U <= 1 + v(v+1)/(m+1) + (v+1)A.
  //
  // So, finally, assume also that v(v+1) << m+1.
  // [This is the same as saying v^2 << m, since m = 2^64-1 is very large,
  // and thus m+1~m, sqrt(m+1)~sqrt(m)].
  //
  // ...then:   A ~ (v^2/m) / m << 1/m,  (v+1)A << v/m << 1,
  // and so U = 1 + C where C << 1.
  //
  // Thus, the bias towards the max value v is negligible, as required.

  // Divide range into approximately equal widths.
  // Notice, we can't do m_engine.max()+1 because it overflows to 0.
  // But the chance of getting m_engine.max() is negligibly small anyway.
  const std::uint64_t interval_width =
      m_engine.max() /
      // Doesn't overflow, because of the above checks.
      (static_cast<std::uint64_t>(max_value) + 1);

  // interval_width cannot be zero, because we ensured above that
  // max_value + 1 <= m_engine.max().
  const size_t result = random_int / interval_width;

  // Modulo arithmetic shouldn't be necessary, but be paranoid,
  // in case there are mistakes in the above analysis (very likely!)
  return result % (max_value + 1);
}

size_t RNG::get_size_t(size_t min_value, size_t max_value) {
  if (min_value > max_value) {
    std::swap(min_value, max_value);
  }
  return min_value + get_size_t(max_value - min_value);
}

vector<size_t> RNG::get_permutation(size_t size) {
  vector<size_t> numbers(size);
  for (size_t i = 0; i < size; ++i) {
    numbers[i] = i;
  }
  do_shuffle(numbers);
  return numbers;
}

void RNG::set_seed(size_t seed) { m_engine.seed(seed); }

bool RNG::check_percentage(size_t percentage) {
  // e.g. the numbers {0,1,2,3,4} are 5%
  // of the numbers {0,1,...,99}.
  return get_size_t(99) < percentage;
}

}  // namespace tket
