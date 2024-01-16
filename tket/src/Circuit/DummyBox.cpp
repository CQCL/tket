// Copyright 2019-2024 Cambridge Quantum Computing
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

#include <boost/lexical_cast.hpp>

#include "Circuit/ResourceData.hpp"
#include "tket/Ops/OpJsonFactory.hpp"

namespace tket {

DummyBox::DummyBox(
    unsigned n_qubits_, unsigned n_bits_, const ResourceData &resource_data_)
    : Box(OpType::DummyBox),
      n_qubits(n_qubits_),
      n_bits(n_bits_),
      resource_data(resource_data_) {}

DummyBox::DummyBox(const DummyBox &other)
    : Box(other),
      n_qubits(other.n_qubits),
      n_bits(other.n_bits),
      resource_data(other.resource_data) {}

bool DummyBox::is_equal(const Op &op_other) const {
  const DummyBox &other = dynamic_cast<const DummyBox &>(op_other);
  if (id_ == other.get_id()) return true;
  return resource_data == other.resource_data;
}

unsigned DummyBox::get_n_qubits() const { return n_qubits; }
unsigned DummyBox::get_n_bits() const { return n_bits; }
ResourceData DummyBox::get_resource_data() const { return resource_data; }

op_signature_t DummyBox::get_signature() const {
  op_signature_t sig(n_qubits, EdgeType::Quantum);
  op_signature_t bits(n_bits, EdgeType::Classical);
  sig.insert(sig.end(), bits.begin(), bits.end());
  return sig;
}

nlohmann::json DummyBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const DummyBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["n_qubits"] = box.get_n_qubits();
  j["n_bits"] = box.get_n_bits();
  j["resource_data"] = box.get_resource_data();
  return j;
}

Op_ptr DummyBox::from_json(const nlohmann::json &j) {
  DummyBox box = DummyBox(
      j.at("n_qubits").get<unsigned>(), j.at("n_bits").get<unsigned>(),
      j.at("resource_data").get<ResourceData>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(DummyBox, DummyBox)

void DummyBox::generate_circuit() const { throw DummyBoxNotDecomposable(); }

}  // namespace tket
