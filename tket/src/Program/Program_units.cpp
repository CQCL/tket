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

#include "Program.hpp"

namespace tket {

qubit_vector_t Program::all_qubits() const {
  qubit_vector_t all_qbs;
  for (auto [it, end] = units_.get<TagType>().equal_range(UnitType::Qubit);
       it != end; ++it) {
    all_qbs.push_back(Qubit(*it));
  }
  return all_qbs;
}

bit_vector_t Program::all_bits() const {
  bit_vector_t all_bs;
  for (auto [it, end] = units_.get<TagType>().equal_range(UnitType::Bit);
       it != end; it++) {
    all_bs.push_back(Bit(*it));
  }
  return all_bs;
}

unit_vector_t Program::all_units() const {
  unit_vector_t all_us;
  for (const UnitID& u : units_.get<TagID>()) {
    all_us.push_back(u);
  }
  return all_us;
}

std::map<Bit, unsigned> Program::bit_readout() const {
  std::map<Bit, unsigned> res;

  // Order bits to generate indices
  bit_vector_t all_bs = all_bits();
  std::sort(all_bs.begin(), all_bs.end());
  unsigned i = 0;
  for (const Bit& b : all_bs) {
    res.insert({b, i});
    i++;
  }

  return res;
}

std::map<Qubit, unsigned> Program::qubit_readout() const {
  std::map<Bit, unsigned> bit_ro = bit_readout();
  FGVertVec finals = get_predecessors(exit_);
  if (finals.size() == 1) {
    // Circuit may not contain every unit from the full program, so need to map
    // indices
    const Circuit& circ = get_circuit_ref(finals.front());
    std::map<Qubit, unsigned> circ_ro = circ.qubit_readout();
    bit_vector_t circ_bits = circ.all_bits();
    std::map<Qubit, unsigned> result;
    for (const std::pair<const Qubit, unsigned>& pair : circ_ro) {
      result.insert({pair.first, bit_ro.at(circ_bits.at(pair.second))});
    }
    return result;
  } else
    return {};
}

opt_reg_info_t Program::get_reg_info(std::string reg_name) const {
  unit_lookup_t::index<TagReg>::type::iterator found =
      units_.get<TagReg>().find(reg_name);
  if (found == units_.get<TagReg>().end())
    return std::nullopt;
  else
    return found->reg_info();
}

register_t Program::get_reg(std::string reg_name) const {
  register_t reg;
  for (auto [it, end] = units_.get<TagReg>().equal_range(reg_name); it != end;
       it++) {
    if (it->reg_dim() != 1)
      throw CircuitInvalidity("Cannot linearise register " + reg_name);
    reg.insert({it->index()[0], *it});
  }
  return reg;
}

void Program::add_qubit(const Qubit& id, bool reject_dups) {
  unit_lookup_t::index<TagID>::type::iterator found =
      units_.get<TagID>().find(id);
  if (found != units_.get<TagID>().end()) {
    if (reject_dups) {
      throw CircuitInvalidity(
          "A unit with ID \"" + id.repr() + "\" already exists");
    } else if (found->type() == UnitType::Qubit) {
      return;
    } else {
      throw CircuitInvalidity(
          "A bit with ID \"" + id.repr() + "\" already exists");
    }
  }
  opt_reg_info_t reg_info = get_reg_info(id.reg_name());
  register_info_t correct_info = {UnitType::Qubit, id.reg_dim()};
  if (reg_info && !(reg_info.value() == correct_info))
    throw CircuitInvalidity(
        "Cannot add qubit with ID \"" + id.repr() +
        "\" as register is not compatible");
  units_.insert(id);
}

void Program::add_bit(const Bit& id, bool reject_dups) {
  unit_lookup_t::index<TagID>::type::iterator found =
      units_.get<TagID>().find(id);
  if (found != units_.get<TagID>().end()) {
    if (reject_dups) {
      throw CircuitInvalidity(
          "A unit with ID \"" + id.repr() + "\" already exists");
    } else if (found->type() == UnitType::Bit) {
      return;
    } else {
      throw CircuitInvalidity(
          "A qubit with ID \"" + id.repr() + "\" already exists");
    }
  }
  opt_reg_info_t reg_info = get_reg_info(id.reg_name());
  register_info_t correct_info = {UnitType::Bit, id.reg_dim()};
  if (reg_info && !(reg_info.value() == correct_info))
    throw CircuitInvalidity(
        "Cannot add bit with ID \"" + id.repr() +
        "\" as register is not compatible");
  units_.insert(id);
}

register_t Program::add_q_register(std::string reg_name, unsigned size) {
  if (get_reg_info(reg_name))
    throw CircuitInvalidity(
        "A register with name \"" + reg_name + "\" already exists");
  register_t ids;
  for (unsigned i = 0; i < size; i++) {
    Qubit id(reg_name, i);
    units_.insert(id);
    ids.insert({i, id});
  }
  return ids;
}

register_t Program::add_c_register(std::string reg_name, unsigned size) {
  if (get_reg_info(reg_name))
    throw CircuitInvalidity(
        "A register with name \"" + reg_name + "\" already exists");
  register_t ids;
  for (unsigned i = 0; i < size; i++) {
    Bit id(reg_name, i);
    units_.insert(id);
    ids.insert({i, id});
  }
  return ids;
}

}  // namespace tket
