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

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>

#include <complex>

#include "Utils/Json.hpp"
#include "binder_json.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

PYBIND11_MODULE(utils_serialization, m) {
  m.def(
      "complex_to_list",
      [](std::complex<double> c) {
        json j = c;
        return j;
      },
      "Convert complex number to serializable list [real, imag].");
  m.def(
      "list_to_complex", [](json j) { return j.get<std::complex<double>>(); },
      "Convert serializable list as output by `complex_to_list` to complex "
      "number.");
}

}  // namespace tket
