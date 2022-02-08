// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "PathFinderInterface.hpp"

#include "Utils/Exceptions.hpp"

namespace tket {
namespace tsa_internal {

PathFinderInterface::PathFinderInterface() : m_name("Empty") {}

PathFinderInterface::~PathFinderInterface() {}

// GCOVR_EXCL_START
const std::vector<size_t>& PathFinderInterface::operator()(
    size_t /*vertex1*/, size_t /*vertex2*/) {
  throw NotImplemented("PathFinderInterface: get path");
}

const std::string& PathFinderInterface::name() const { return m_name; }

void PathFinderInterface::reset() {}

void PathFinderInterface::register_edge(
    size_t /*vertex1*/, size_t /*vertex2*/) {}

bool PathFinderInterface::edge_registration_has_effect() const { return false; }
// GCOVR_EXCL_STOP

}  // namespace tsa_internal
}  // namespace tket
