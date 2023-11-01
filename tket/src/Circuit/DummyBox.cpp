// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tket/Circuit/DummyBox.hpp"

#include "Circuit/ResourceData.hpp"

namespace tket {

DummyBox::DummyBox(const ResourceData &resource_data_)
    : Box(OpType::DummyBox), resource_data(resource_data_) {}

DummyBox::DummyBox(const DummyBox &other)
    : Box(other), resource_data(other.resource_data) {}

ResourceData DummyBox::get_resource_data() const { return resource_data; }
void DummyBox::generate_circuit() const { throw DummyBoxNotDecomposable(); }

}  // namespace tket
