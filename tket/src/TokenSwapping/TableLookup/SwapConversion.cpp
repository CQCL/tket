#include "SwapConversion.hpp"

#include "Utils/Assert.hpp"

;
using std::vector;

namespace tket {
namespace tsa_internal {

static vector<Swap> get_swaps_fixed_vector() {
  vector<Swap> swaps;
  for (unsigned ii = 0; ii < 6; ++ii) {
    for (unsigned jj = ii + 1; jj < 6; ++jj) {
      swaps.push_back(get_swap(ii, jj));
    }
  }
  TKET_ASSERT(swaps.size() == 15);
  return swaps;
}

static const vector<Swap>& get_swaps_global() {
  static const auto swaps_vect(get_swaps_fixed_vector());
  return swaps_vect;
}

const Swap& SwapConversion::get_swap_from_hash(SwapHash x) {
  return get_swaps_global().at(x - 1);
}

static std::map<Swap, SwapConversion::SwapHash> get_swap_to_hash() {
  const auto swaps = get_swaps_fixed_vector();
  std::map<Swap, SwapConversion::SwapHash> map;
  for (unsigned ii = 0; ii < swaps.size(); ++ii) {
    map[swaps[ii]] = ii + 1;
  }
  return map;
}

static const std::map<Swap, SwapConversion::SwapHash>&
get_swap_to_hash_global() {
  static const auto map(get_swap_to_hash());
  return map;
}

SwapConversion::SwapHash SwapConversion::get_hash_from_swap(const Swap& swap) {
  return get_swap_to_hash_global().at(swap);
}

unsigned SwapConversion::get_number_of_swaps(
    SwapConversion::SwapHash swaps_code) {
  unsigned num_swaps = 0;
  while (swaps_code != 0) {
    ++num_swaps;
    const auto swap_hash = swaps_code & 0xF;
    swaps_code >>= 4;
    TKET_ASSERT(swap_hash > 0);
    TKET_ASSERT(swap_hash <= 15);
  }
  return num_swaps;
}

SwapConversion::EdgesBitset SwapConversion::get_edges_bitset(
    SwapHash swaps_code) {
  EdgesBitset edges_bitset = 0;
  while (swaps_code != 0) {
    const auto swap_hash = swaps_code & 0xF;
    TKET_ASSERT(swap_hash > 0);
    edges_bitset |= (1u << (swap_hash - 1));
    swaps_code >>= 4;
  }
  return edges_bitset;
}

}  // namespace tsa_internal
}  // namespace tket
