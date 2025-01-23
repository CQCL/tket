// Copyright Quantinuum
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

#include <pybind11/pybind11.h>
namespace py = pybind11;

template <class T>
int deletedHash(const T&) {
  py::type T_py = py::type::of<T>();
  auto T_name = T_py.attr("__name__").cast<std::string>();
  auto error_message = std::string("unhashable type: '") + T_name + "'";
  throw py::type_error(error_message);
}

constexpr auto deletedHashDocstring =
    "Hashing is not implemented for this class, attempting to hash an object "
    "will raise a type error";
