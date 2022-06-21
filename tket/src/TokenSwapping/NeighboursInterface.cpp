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

#include "NeighboursInterface.hpp"

#include <stdexcept>

namespace tket {

const std::vector<size_t>& NeighboursInterface::operator()(size_t) {
  throw std::logic_error("NeighboursInterface::operator() not implemented");
}

NeighboursInterface::~NeighboursInterface() {}

}  // namespace tket
