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

#pragma once

/**
 * @file
 * @brief Named registers of arrays of (quantum or classical) nodes
 */

#include <boost/functional/hash.hpp>
#include <map>
#include <memory>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <string>

#include "BiMapHeaders.hpp"
#include "Json.hpp"
#include "TketLog.hpp"

namespace tket {

/** Type of information held */
enum class UnitType { Qubit, Bit };

/** The type and dimension of a register */
typedef std::pair<UnitType, unsigned> register_info_t;

typedef std::optional<register_info_t> opt_reg_info_t;

const std::string &q_default_reg();
const std::string &c_default_reg();
const std::string &node_default_reg();
const std::string &c_debug_zero_prefix();
const std::string &c_debug_one_prefix();
const std::string &c_debug_default_name();

/** Conversion invalid */
class InvalidUnitConversion : public std::logic_error {
 public:
  InvalidUnitConversion(const std::string &name, const std::string &new_type)
      : std::logic_error("Cannot convert " + name + " to " + new_type) {}
};

/**
 * Location holding a bit or qubit of information
 *
 * Each location has a name (signifying the 'register' to which it belongs) and
 * an index within that register (which may be multi-dimensional).
 */
class UnitID {
 public:
  UnitID() : data_(std::make_shared<UnitData>()) {}

  /** String representation including name and index */
  std::string repr() const;

  /** Register name */
  std::string reg_name() const { return data_->name_; }

  /** Index dimension */
  unsigned reg_dim() const { return data_->index_.size(); }

  /** Index */
  std::vector<unsigned> index() const { return data_->index_; }

  /** Unit type */
  UnitType type() const { return data_->type_; }

  /** Register dimension and type */
  register_info_t reg_info() const { return {type(), reg_dim()}; }

  bool operator<(const UnitID &other) const {
    int n = data_->name_.compare(other.data_->name_);
    if (n > 0) return false;
    if (n < 0) return true;
    return data_->index_ < other.data_->index_;
  }
  bool operator>(const UnitID &other) const { return other < *this; }
  bool operator==(const UnitID &other) const {
    return (this->data_->name_ == other.data_->name_) &&
           (this->data_->index_ == other.data_->index_);
  }
  bool operator!=(const UnitID &other) const { return !(*this == other); }

  friend std::size_t hash_value(UnitID const &unitid) {
    size_t seed = 0;
    boost::hash_combine(seed, unitid.data_->name_);
    boost::hash_combine(seed, unitid.data_->index_);
    boost::hash_combine(seed, unitid.data_->type_);
    return seed;
  }

 protected:
  UnitID(
      const std::string &name, const std::vector<unsigned> &index,
      UnitType type)
      : data_(std::make_shared<UnitData>(name, index, type)) {}

 private:
  struct UnitData {
    std::string name_;
    std::vector<unsigned> index_;
    UnitType type_;

    UnitData() : name_(), index_(), type_(UnitType::Qubit) {}
    UnitData(
        const std::string &name, const std::vector<unsigned> &index,
        UnitType type)
        : name_(name), index_(index), type_(type) {
      static const std::string id_regex_str = "[a-z][A-Za-z0-9_]*";
      static const std::regex id_regex(id_regex_str);
      if (!name.empty() && !std::regex_match(name, id_regex)) {
        std::stringstream msg;
        msg << "UnitID name '" << name << "' does not match '" << id_regex_str
            << "', as required for QASM conversion.";
        tket_log()->warn(msg.str());
      }
    }
  };
  std::shared_ptr<UnitData> data_;
};

template <class Unit_T>
void unitid_to_json(nlohmann::json &j, const Unit_T &unit) {
  static_assert(std::is_base_of<UnitID, Unit_T>::value);
  j.push_back(unit.reg_name());
  j.push_back(unit.index());
}

template <class T>
void json_to_unitid(const nlohmann::json &j, T &unit) {
  unit = T(j.at(0).get<std::string>(), j.at(1).get<std::vector<unsigned>>());
}

/** Location holding a qubit */
class Qubit : public UnitID {
 public:
  Qubit() : UnitID("", {}, UnitType::Qubit) {}

  /** Qubit in default register */
  explicit Qubit(unsigned index)
      : UnitID(q_default_reg(), {index}, UnitType::Qubit) {}

  /** Named register with no index */
  explicit Qubit(const std::string &name) : UnitID(name, {}, UnitType::Qubit) {}

  /** Named register with a one-dimensional index */
  Qubit(const std::string &name, unsigned index)
      : UnitID(name, {index}, UnitType::Qubit) {}

  /** Named register with a two-dimensional index */
  Qubit(const std::string &name, unsigned row, unsigned col)
      : UnitID(name, {row, col}, UnitType::Qubit) {}

  /** Named register with a three-dimensional index */
  Qubit(const std::string &name, unsigned row, unsigned col, unsigned layer)
      : UnitID(name, {row, col, layer}, UnitType::Qubit) {}

  /** Named register with a multi-dimensional index */
  Qubit(const std::string &name, std::vector<unsigned> index)
      : UnitID(name, index, UnitType::Qubit) {}

  /** Copy constructor */
  explicit Qubit(const UnitID &other) : UnitID(other) {
    if (other.type() != UnitType::Qubit) {
      throw InvalidUnitConversion(other.repr(), "Qubit");
    }
  }
};

JSON_DECL(Qubit)

/** Location holding a bit */
class Bit : public UnitID {
 public:
  Bit() : UnitID("", {}, UnitType::Bit) {}

  /** Bit in default register */
  explicit Bit(unsigned index)
      : UnitID(c_default_reg(), {index}, UnitType::Bit) {}

  /** Named register with no index */
  explicit Bit(const std::string &name) : UnitID(name, {}, UnitType::Bit) {}

  /** Named register with a one-dimensional index */
  Bit(const std::string &name, unsigned index)
      : UnitID(name, {index}, UnitType::Bit) {}

  /** Named register with a two-dimensional index */
  Bit(const std::string &name, unsigned row, unsigned col)
      : UnitID(name, {row, col}, UnitType::Bit) {}

  /** Named register with a three-dimensional index */
  Bit(const std::string &name, unsigned row, unsigned col, unsigned layer)
      : UnitID(name, {row, col, layer}, UnitType::Bit) {}

  /** Named register with a multi-dimensional index */
  Bit(const std::string &name, std::vector<unsigned> index)
      : UnitID(name, index, UnitType::Bit) {}

  explicit Bit(const UnitID &other) : UnitID(other) {
    if (other.type() != UnitType::Bit) {
      throw InvalidUnitConversion(other.repr(), "Bit");
    }
  }
};

JSON_DECL(Bit)

/** Architectural qubit location */
class Node : public Qubit {
 public:
  Node() : Qubit() {}

  /** Qubit in default register */
  explicit Node(unsigned index) : Qubit(node_default_reg(), index) {}

  /** Named register with a one-dimensional index */
  Node(const std::string &name, unsigned index) : Qubit(name, index) {}

  /** Named register with a two-dimensional index */
  Node(const std::string &name, unsigned row, unsigned col)
      : Qubit(name, row, col) {}

  /** Named register with a three-dimensional index */
  Node(const std::string &name, unsigned row, unsigned col, unsigned layer)
      : Qubit(name, row, col, layer) {}

  /** Named register with a multi-dimensional index */
  Node(const std::string &name, std::vector<unsigned> index)
      : Qubit(name, index) {}

  explicit Node(const UnitID &other) : Qubit(other) {}
};

JSON_DECL(Node)

/** A correspondence between two sets of unit IDs */
typedef boost::bimap<UnitID, UnitID> unit_bimap_t;

/** A pair of ("initial" and "final") correspondences between unit IDs */
typedef struct {
  unit_bimap_t initial;
  unit_bimap_t final;
} unit_bimaps_t;

typedef std::vector<UnitID> unit_vector_t;
typedef std::map<UnitID, UnitID> unit_map_t;
typedef std::set<UnitID> unit_set_t;

typedef std::vector<Qubit> qubit_vector_t;
typedef std::map<Qubit, Qubit> qubit_map_t;
JSON_DECL(qubit_map_t)

typedef std::vector<Bit> bit_vector_t;
typedef std::map<Bit, Bit> bit_map_t;

typedef std::set<Node> node_set_t;
typedef std::vector<Node> node_vector_t;

/** A register of locations sharing the same name */
typedef std::map<unsigned, UnitID> register_t;

template <typename UnitA, typename UnitB>
static bool update_map(unit_bimap_t &m, const std::map<UnitA, UnitB> &um) {
  unit_map_t new_m;
  bool changed = false;
  for (const std::pair<const UnitA, UnitB> &pair : um) {
    const auto &it = m.right.find(pair.first);
    if (it == m.right.end()) {
      continue;
    }
    new_m.insert({it->second, pair.second});
    changed |= (m.right.erase(pair.first) > 0);
  }
  for (const std::pair<const UnitID, UnitID> &pair : new_m) {
    changed |= m.left.insert(pair).second;
  }
  return changed;
}

/**
 * Update a pair of "initial" and "final" correspondences.
 *
 * If \p maps is null then the function does nothing and returns false.
 *
 * @param[in,out] maps maps to be updated
 * @param[in] um_initial new correspondences added to initial map
 * @param[in] um_final new correspondences added to final map
 *
 * @tparam UnitA first unit type
 * @tparam UnitB second unit type
 *
 * @return whether any changes were made to the maps
 */
template <typename UnitA, typename UnitB>
bool update_maps(
    std::shared_ptr<unit_bimaps_t> maps,
    const std::map<UnitA, UnitB> &um_initial,
    const std::map<UnitA, UnitB> &um_final) {
  if (!maps) return false;
  // Can only work for Unit classes
  static_assert(std::is_base_of<UnitID, UnitA>::value);
  static_assert(std::is_base_of<UnitID, UnitB>::value);
  // Unit types must be related, so cannot rename e.g. Bits to Qubits
  static_assert(
      std::is_base_of<UnitA, UnitB>::value ||
      std::is_base_of<UnitB, UnitA>::value);

  bool changed = false;
  changed |= update_map(maps->initial, um_initial);
  changed |= update_map(maps->final, um_final);
  return changed;
}

}  // namespace tket
