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

/*
 * This function template implements python equality for any class that
 * implements ==
 * -> Using py::object as the parameter solves mypy error in the stubs generated
 * when equality is defined only for a subset of the python objects (which is
 * frowned upon due to pythons lack of strict typing)
 * -> The typical way to do this in pybind11 leads to this error
 */

template <typename T>
bool py_equals(const T& self, const py::object& pyObject) {
  if (py::isinstance<T>(pyObject)) {
    T* cppObject = pyObject.cast<T*>();
    return *cppObject == self;
  }
  return false;
}

template <typename T>
bool py_not_equals(const T& self, const py::object& pyObject) {
  if (py::isinstance<T>(pyObject)) {
    T* cppObject = pyObject.cast<T*>();
    return *cppObject != self;
  }
  return true;
}
