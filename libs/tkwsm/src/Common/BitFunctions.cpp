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

#include "tkwsm/Common/BitFunctions.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

unsigned BitFunctions::get_number_of_rightmost_zero_bits(std::uint64_t x) {
  // Many similar things to be found in
  // http://graphics.stanford.edu/~seander/bithacks.html
  if (x == 0) {
    return 64;
  }
  // Will be taken away at the end, if necessary.
  unsigned result = 1;

  // A similar argument is given in get_bit_length below.
  //
  // At each "if" statement, for a shift of S=2^k bits:
  // if it returns TRUE, then x has >= S trailing zeros BEFORE,
  // and <S zeros AFTER.
  // Otherwise, x has <S zeros.
  // So IN ALL CASES, x has <S trailing zeros just AFTER the "if" statement.
  // Also,  Trailing_zeros(x) + result  is UNCHANGED by each if statement.

  // We only ever destroy rightmost zeros.
  if ((x & 0xffffffff) == 0) {
    result += 32;
    x >>= 32;
  }
  if ((x & 0xffff) == 0) {
    result += 16;
    x >>= 16;
  }
  if ((x & 0xff) == 0) {
    result += 8;
    x >>= 8;
  }
  if ((x & 0xf) == 0) {
    result += 4;
    x >>= 4;
  }
  if ((x & 3) == 0) {
    result += 2;
    x >>= 2;
  }
  // See the above arguments: now
  //
  //    Trailing_zeros(x) + result  =  Trailing_zeros(orig x) + (orig result)
  //                                =  Trailing_zeros(orig x) + 1,
  //   Trailing_zeros(x) < 2.
  //
  // So Trailing_zeros(x) = 0 or 1.
  result -= (x & 1);
  return result;
}

unsigned BitFunctions::get_bit_length(std::uint64_t x) {
  if (x == 0) {
    return 0;
  }
  unsigned result = 1;

  // To prove correctness of the following:
  //
  // It's clear that x cannot be reduced to zero
  // by the "if" statements, since they
  // never destroy the leftmost bit of x.
  //
  // Thus, before/after each if statement,
  //
  //    "Bitlength(x) + result" is unchanged.
  //
  // Thus, after each "if" statement,
  //
  //    Bitlength(new x) + result = Bitlength(original x) + (original result)
  //                              = Bitlength(original x) + 1.
  //
  // Finally, by induction (the base case uses the fact that it's 64 bits!):
  //
  // For each "if" statement which possibly shifts x by k bits:
  //  Bitlength(old x) <= 2k,   Bitlength(new x) <= k,
  // whether "if" returned true or false.
  //
  if (x > 0xffffffff) {
    result += 32;
    x >>= 32;
  }
  if (x > 0xffff) {
    result += 16;
    x >>= 16;
  }
  if (x > 0xff) {
    result += 8;
    x >>= 8;
  }
  if (x > 0xf) {
    result += 4;
    x >>= 4;
  }
  if (x > 3) {
    result += 2;
    x >>= 2;
  }
  // By the above arguments,
  //    Bitlength(x) <= 2,
  // so x=1,2,3 (because x>0).
  //
  // Thus (at this stage), x=1 <==> Bitlength(x) = 1,
  //
  // and   result = Bitlength(original x) + 1 - Bitlength(x).
  //
  if (x > 1) {
    ++result;
  }
  return result;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
