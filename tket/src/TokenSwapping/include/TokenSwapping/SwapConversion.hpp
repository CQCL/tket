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

#include <cstdint>

#include "TokenSwapping/SwapFunctions.hpp"

namespace tket {
namespace tsa_internal {

/*
NOTE on ENCODING: with 6 vertices, there are 15 possible edges or swaps.
Thus, we can encode a single swap by a number in the range 0-15 (using 0 to
denote "no swap").

This fits into 4 bits exactly.

Thus, a single 64-bit unsigned int can store any swap sequence of length <= 16.
We also have the added benefit that ints written in hexadecimal are easier for a
human to read, since each hex digit 0-9 or A-F corresponds to a single swap.

An obvious optimisation is that adjacent swaps should be different;
and also, blocks of four zeros cannot occur within the encoding.
However, this would still only reduce the total number to about 30%,
so we'd still need 62 or 63 bits to represent all sequences of length <= 16.
So it's not worth trying fancy encodings to store more possible sequences in
fewer bits, without a good theoretical breakthrough to come up with a really
good way to encode and search through only optimal or "near optimal" sequences.

If we desire in future to increase the number of vertices, we'd have to use at
least 5 bits per swap, so could only fit sequences of length <= 12 in a 64-bit
int. Of course, (8*7)/2 = 28 < 31, so we could store swaps on <= 8 vertices
instead of 6.

*/

// Generally no checks on the input values, it's assumed that the caller
// knows how the table encoding works.
// The possible swaps (01), (02), (03), ..., (45) on vertices {0,1,2,3,4,5}
// are listed in a global vector, so with values 0,1,...,14.
// Adding 1 to the index gives possible values 1,2,...,15 for the swaps,
// and 0 means no swap. Thus a sequence of swaps is encoded by storing the bits
// in a uint, with first swap at the least significant bits, and so on with
// leftward shifts by 4 bits each time.

struct SwapConversion {
  /** Encodes a sequence of <=16 swaps, each swap being one of
   * the 15 possible swaps on vertices {0,1,2,3,4,5}, and hence encoded by 4
   * bits. Zero represents the empty sequence. */
  typedef std::uint64_t SwapHash;

  /** Encodes a set of swaps, each one taken from the 15 possibilities. With
   * each swap given a numerical value from 1 to 15, we simply shift 1u by that
   * amount (minus one), and OR them together. Thus, when looking up in a table,
   * we only allow swap sequences whose edge bitset is a SUBSET of a given edges
   * bitset (corresponding to the edges in the graph, i.e. allowed swaps).
   */
  typedef std::uint_fast16_t EdgesBitset;

  /** Given a valid number x, return the actual swap on vertices {0,1,2,3,4,5}
   * which it represents.
   * @param x A code number representing a single swap.
   * @return A single swap on vertices {0,1,2,3,4,5}.
   */
  static const Swap& get_swap_from_hash(SwapHash x);

  /** The opposite of get_swap_from_hash.
   * @param swap A swap on {0,1,2,3,4,5}. (Must be in standard order, i.e. (i,j)
   * with 0 <= i < j <= 5).
   * @return A number 1-15 which encodes that swap in the table.
   */
  static SwapHash get_hash_from_swap(const Swap& swap);

  /** Converting swaps to bitsets, which swaps are used in the code?
   * @param swaps_code An integer representing a sequence of swaps.
   * @return The set of swaps used in the sequence, encoded as a binary number.
   */
  static EdgesBitset get_edges_bitset(SwapHash swaps_code);

  /** The number of swaps in a sequence.
   * @param swaps_code An integer representing a sequence of swaps.
   * @return The length of the swap sequence.
   */
  static unsigned get_number_of_swaps(SwapHash swaps_code);
};

}  // namespace tsa_internal
}  // namespace tket
