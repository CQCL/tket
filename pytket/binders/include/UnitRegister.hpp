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

#include "Utils/UnitID.hpp"

namespace tket {
/*
Bit registers which can be interpreted as unsigned integers follow the
conventions defined here:
registers are up to _TKET_REG_WIDTH wide in bits and are interpreted as
equivalent to the C++ type _tket_uint_t
*/
#define _TKET_REG_WIDTH 32
typedef uint32_t _tket_uint_t;

template <typename T>
class UnitRegister {
 public:
  UnitRegister(const std::string &name, const std::size_t size)
      : name_(name), size_(size){};

  std::string name() const { return name_; }
  std::size_t size() const { return size_; }

  void set_name(const std::string &new_name) { name_ = new_name; }
  void set_size(const std::size_t &new_size) { size_ = new_size; }

  bool contains(const T &unit) const {
    const std::vector<unsigned> index = unit.index();
    return (unit.reg_name() == name_) && (index.size() == 1) &&
           (index[0] < size_);
  }

  const T operator[](std::size_t index) const {
    if (index >= size_) {
      throw std::out_of_range("Index out of range of UnitRegister.");
    };
    return T(name_, index);
  }
  bool operator==(const UnitRegister<T> &other) const {
    return (name_ == other.name_) && (size_ == other.size_);
  }
  bool operator<(const UnitRegister<T> &other) const {
    if (name_ == other.name_) {
      return size_ < other.size_;
    } else {
      return name_ < other.name_;
    }
  }

 private:
  std::string name_;
  std::size_t size_;
};

typedef UnitRegister<Bit> BitRegister;
typedef UnitRegister<Qubit> QubitRegister;
}  // namespace tket
