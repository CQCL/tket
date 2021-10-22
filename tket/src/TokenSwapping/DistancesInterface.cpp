#include "DistancesInterface.hpp"

;
using std::vector;

namespace tket {
namespace tsa_internal {

void DistancesInterface::register_shortest_path(
    const vector<size_t>& /*path*/) {}

void DistancesInterface::register_neighbours(
    size_t vertex, const vector<size_t>& neighbours) {
  for (size_t nv : neighbours) {
    register_edge(vertex, nv);
  }
}

void DistancesInterface::register_edge(size_t /*vertex1*/, size_t /*vertex2*/) {
}

DistancesInterface::~DistancesInterface() {}

}  // namespace tsa_internal
}  // namespace tket
