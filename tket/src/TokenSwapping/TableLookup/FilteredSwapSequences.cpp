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

#include "TokenSwapping/FilteredSwapSequences.hpp"

#include <algorithm>
#include <limits>

#include "TokenSwapping/GeneralFunctions.hpp"
#include "TokenSwapping/SwapSequenceTable.hpp"
#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

/*
NOTE: the problem is: given a bitset, i.e. an unsigned int representing a set,
design a map-type data structure whose keys are unsigned integers representing
bitsets, and values are a collection of entries using that bitset (i.e., only
using swaps whose index in a global vector of allowed swaps has a "one" in the
appropriate position in the binary expansion of the bitset).

We must be able to look up all entries whose key is a SUBSET of the given set.
(And then, search further through those values).

We tried various things, e.g. sorting by key, using the fact that

(X is a subset of Y) ==> (x <= y)

where X,Y are subsets and x,y are the integers representing them; thus you can
do a kind of binary search.

If you want the SMALLEST value for a given key, you can sort them also and do a
kind of double binary search. (Another crucial point: when searching between two
key ranges in a sorted VECTOR of keys, you can determine how many keys exist in
the range in O(log N) time, rather than O(N) time for a map).

These fancy algorithms are all asymptotically much more efficient than the
obvious O(N) lookup, which just goes through EVERY key and checks if it's a
subset or not, then goes through every element.

HOWEVER, experiments showed that the fancy algorithms are actually quite a bit
slower than the obvious algorithm for the table size we care about.

*/

FilteredSwapSequences::SingleSequenceData::SingleSequenceData()
    : edges_bitset(0),
      swaps_code(0),
      number_of_swaps(std::numeric_limits<unsigned>::max()) {}

/*
If the entries are distributed "randomly" and fairly uniformly amongst the
bitset keys, i.e. given a bitset, look up all keys which are a subset of that,
then asymptotically using many bits in the keys is good.

For our table sizes, experiments suggested that it's worth having 1 bit in each
bitset key (2 min for 1 bit vs. 2 min 20 sec for no bits in one test), rather
then no keys at all, BUT not worth more than 1 bit in each key.

e.g., for 15 bits in each bitset, each of the 15 keys being one of the bits
(we have no empty keys - pointless trying to look up swap sequences if the graph
has no edges!), assume that an average lookup query contains 5 bits. Then 10/15
= 2/3 of the keys are disjoint from it, and so most of the keys immediately can
be ruled out.

However, it's a balancing act: if you have too many keys, then the lists for
each key become so short then you're effectively almost doing a linear search
through all entries.
*/

void FilteredSwapSequences::initialise(
    std::vector<SwapConversion::SwapHash> codes) {
  // Can only initialise once.
  TKET_ASSERT(m_internal_data.empty());
  std::sort(codes.begin(), codes.end());
  TKET_ASSERT(!codes.empty());
  TKET_ASSERT(codes[0] != 0);
  TrimmedSingleSequenceData datum;

  for (size_t ii = 0; ii < codes.size(); ++ii) {
    if (ii != 0 && codes[ii] == codes[ii - 1]) {
      // Filter out duplicate entries.
      continue;
    }
    datum.swaps_code = codes[ii];
    datum.edges_bitset = SwapConversion::get_edges_bitset(datum.swaps_code);
    push_back(datum);
  }
}

void FilteredSwapSequences::push_back(TrimmedSingleSequenceData datum) {
  auto bitset_copy = datum.edges_bitset;
  TKET_ASSERT(bitset_copy != 0);
  SwapConversion::EdgesBitset bit_to_use = 0;

  // We want to add to the smallest list, to keep the data balanced.
  // Tests showed that this works well; the entries are distributed
  // very close to uniformly amongst the 15 possible keys.
  //
  // This is maybe surprising, because you'd expect
  // more bias: you'd expect, due to the relabelling scheme, the table to have
  // swaps like (0,1), (0,2) much more frequently than higher-numbered
  // vertices like (4,5). This may or may not be the case, but whatever
  // the truth, there are still enough bits available overall to break
  // the entries up well enough).
  size_t list_size_to_use = std::numeric_limits<size_t>::max();

  while (bitset_copy != 0) {
    const auto new_bit = get_rightmost_bit(bitset_copy);
    // If the key does not exist, the newly created empty list will
    // immediately be filled; so no key is wasted. (They're not wasted anyway,
    // the table entries are very close to uniformly distributed
    // amongst all 15 keys).
    const auto list_size = m_internal_data[new_bit].size();

    if (list_size < list_size_to_use) {
      list_size_to_use = list_size;
      bit_to_use = new_bit;
      if (list_size == 0) {
        break;
      }
    }
  }
  TKET_ASSERT(bit_to_use != 0);
  m_internal_data[bit_to_use].push_back(datum);
}

FilteredSwapSequences::SingleSequenceData
FilteredSwapSequences::get_lookup_result(
    SwapConversion::EdgesBitset edges_bitset, unsigned max_num_swaps) const {
  // NOTE: this algorithm is quite crude, BUT it's so simple that
  // apparently clever algorithms, although asymptotically more efficient,
  // appear to be slower.
  // The clever algorithms seem only worth doing if the table becomes
  // much larger, >> 100 codes for each bit at least.

  max_num_swaps = std::min(max_num_swaps, 16u);

  // Value 0xFFF...F will never occur,
  // because this would be 16 consecutive equal swaps...!
  const auto impossible_max_code =
      std::numeric_limits<SwapConversion::SwapHash>::max();

  // Stop as soon as the swaps code gets too big.
  SwapConversion::SwapHash max_code;
  if (max_num_swaps == 16) {
    max_code = impossible_max_code;
  } else {
    max_code = 1;
    max_code <<= (4 * max_num_swaps);
    --max_code;
  }
  TrimmedSingleSequenceData best_datum;
  best_datum.swaps_code = impossible_max_code;

  for (const auto& entry : m_internal_data) {
    if (entry.first > edges_bitset) {
      // The swaps used by a sequence must be a SUBSET of the allowable edges.
      // Therefore, the swaps bitset must be <= the edges bitset.
      // Of course, it's a MAP, so the swaps bitsets are already in increasing
      // order.
      break;
    }
    if ((entry.first & edges_bitset) != entry.first) {
      // Every swap sequence in this entry contains ALL of the given edges
      // in the bitset key (as well as others), and thus it MUST be a subset
      // of the given edges_bitset, otherwise the entire entry
      // can be skipped.
      continue;
    }
    const auto& list = entry.second;
    for (const auto& single_entry : list) {
      if (single_entry.swaps_code > max_code ||
          single_entry.swaps_code >= best_datum.swaps_code) {
        // Because they're sorted by code value,
        // all subsequent entries will be too big also.
        break;
      }
      if ((single_entry.edges_bitset & edges_bitset) !=
          single_entry.edges_bitset) {
        // The EXACT set of edges used must be a subset of edges_bitset,
        // otherwise it's unsuitable - it uses a swap not allowed.
        continue;
      }
      best_datum = single_entry;
    }
  }

  SingleSequenceData result;
  if (best_datum.swaps_code < impossible_max_code) {
    // We actually got a result.
    result.edges_bitset = best_datum.edges_bitset;
    result.swaps_code = best_datum.swaps_code;
    result.number_of_swaps =
        SwapConversion::get_number_of_swaps(result.swaps_code);
  }
  return result;
}

size_t FilteredSwapSequences::get_total_number_of_entries() const {
  size_t total = 0;
  for (const auto& entry : m_internal_data) {
    total += entry.second.size();
  }
  return total;
}

// Convert the raw SwapSequenceTable object into
// FilteredSwapSequences-compatible data. The key is the permutation hash; the
// value is the lookup object which can find solutions to given problems.
static std::map<unsigned, FilteredSwapSequences>
construct_and_return_full_table() {
  std::map<unsigned, FilteredSwapSequences> result;
  const auto raw_table = SwapSequenceTable::get_table();
  for (const auto& entry : raw_table) {
    // The simplest nontrivial permutation arises from a single swap (a,b),
    // which under the canonical relabelling is converted to (01),
    // which has hash 2.
    TKET_ASSERT(entry.first >= 2);
    // The largest possible hash comes from (01)(23)(45).
    TKET_ASSERT(entry.first <= 222);
    result[entry.first].initialise(entry.second);
  }
  return result;
}

static const std::map<unsigned, FilteredSwapSequences>& get_full_table() {
  static const auto full_table(construct_and_return_full_table());
  return full_table;
}

FilteredSwapSequences::SingleSequenceData::SingleSequenceData(
    unsigned permutation_hash, SwapConversion::EdgesBitset edges_bitset,
    unsigned max_number_of_swaps)
    : SingleSequenceData() {
  if (permutation_hash == 0) {
    // The identity mapping, always possible.
    number_of_swaps = 0;
    return;
  }
  if (edges_bitset == 0) {
    // No swaps at all! This CAN happen...it just means that
    // we haven't seen enough vertices to connect up the given ones;
    // all solutions involve swaps using other vertices not yet seen
    // (i.e., not in this subgraph).
    // But it's not the identity, therefore it's impossible.
    return;
  }

  const auto& table = get_full_table();
  const auto citer = table.find(permutation_hash);
  if (citer == table.cend()) {
    // No result in the table.
    return;
  }
  *this = citer->second.get_lookup_result(edges_bitset, max_number_of_swaps);
}

}  // namespace tsa_internal
}  // namespace tket
