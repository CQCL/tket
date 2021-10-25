#include "SwapFunctions.hpp"

#include <sstream>
#include <stdexcept>

;

namespace tket {
namespace tsa_internal {

Swap get_swap(size_t v1, size_t v2) {
  if (v1 == v2) {
    std::stringstream ss;
    ss << "get_swap : for equal vertices v1 = v2 = v_" << v1;
    throw std::runtime_error(ss.str());
  }
  if (v1 < v2) {
    return std::make_pair(v1, v2);
  }
  return std::make_pair(v2, v1);
}

bool disjoint(const Swap& s1, const Swap& s2) {
  return s1.first != s2.first && s1.first != s2.second &&
         s1.second != s2.first && s1.second != s2.second;
}

}  // namespace tsa_internal
}  // namespace tket
