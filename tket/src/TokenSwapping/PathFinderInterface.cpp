#include "PathFinderInterface.hpp"

#include "Utils/Exceptions.hpp"

namespace tket {
namespace tsa_internal {

PathFinderInterface::PathFinderInterface() : m_name("Empty") {}

PathFinderInterface::~PathFinderInterface() {}

const std::vector<size_t>& PathFinderInterface::operator()(
    size_t /*vertex1*/, size_t /*vertex2*/) {
  throw NotImplemented("PathFinderInterface: get path");
}

const std::string& PathFinderInterface::name() const { return m_name; }

void PathFinderInterface::reset() {}

void PathFinderInterface::register_edge(
    size_t /*vertex1*/, size_t /*vertex2*/) {}

bool PathFinderInterface::edge_registration_has_effect() const { return false; }

}  // namespace tsa_internal
}  // namespace tket
