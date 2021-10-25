#include "DebugFunctions.hpp"

#include <sstream>

namespace tket {
namespace tsa_internal {

std::string str(const VertexMapping& vertex_mapping) {
  std::stringstream ss;
  ss << "VM:";
  for (const auto& entry : vertex_mapping) {
    ss << " " << entry.first << "->" << entry.second << " ";
  }
  return ss.str();
}

std::string str(const SwapList& swaps) { return str(swaps.to_vector()); }

std::string str(const std::vector<Swap>& swaps) {
  std::stringstream ss;
  for (auto swap : swaps) {
    ss << " (" << swap.first << "," << swap.second << ") ";
  }
  return ss.str();
}

}  // namespace tsa_internal
}  // namespace tket
