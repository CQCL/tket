#include "GeneralFunctions.hpp"

#include <numeric>
#include <stdexcept>

;

namespace tket {
namespace tsa_internal {

std::set<size_t> get_random_set(
    RNG& rng, size_t sample_size, size_t population_size) {
  if (sample_size > population_size) {
    throw std::runtime_error("get_random_set: sample too large");
  }
  std::set<size_t> result;
  if (sample_size == 0 || population_size == 0) {
    return result;
  }
  if (sample_size < population_size / 2) {
    while (result.size() < sample_size) {
      result.insert(rng.get_size_t(population_size - 1));
    }
    return result;
  }
  std::vector<size_t> elems(population_size);
  std::iota(elems.begin(), elems.end(), 0);
  rng.do_shuffle(elems);
  for (const auto& elem : elems) {
    result.insert(elem);
    if (result.size() == sample_size) {
      return result;
    }
  }
  throw std::runtime_error("get_random_set: dropped out of loop");
}

}  // namespace tsa_internal
}  // namespace tket
