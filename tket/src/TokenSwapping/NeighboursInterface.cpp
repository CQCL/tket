#include "NeighboursInterface.hpp"

#include "Utils/Exceptions.hpp"

namespace tket {
namespace tsa_internal {

const std::vector<size_t>& NeighboursInterface::operator()(size_t) {
  throw NotImplemented("NeighboursInterface::get_neighbours: not implemented");
}

NeighboursInterface::~NeighboursInterface() {}

}  // namespace tsa_internal
}  // namespace tket
