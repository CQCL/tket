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

#include <catch2/catch_test_macros.hpp>

#include "Characterisation/DeviceCharacterisation.hpp"
#include "OpType/OpDesc.hpp"

namespace tket {
namespace test_DeviceCharacterisation {

SCENARIO("Does Device correctly hold and return error information") {
  GIVEN("Errors for only two qubits, CX & X") {
    // Architecture object for Device object
    Node n0{0};
    Node n1{1};
    using Connections = std::vector<Architecture::Connection>;
    Architecture test_architecture(Connections{{n0, n1}});
    // Create test node errors & link errors for creating Device object
    gate_error_t single_gate_error(0.3);
    op_errors_t test_node_error({{OpType::X, single_gate_error}});

    op_node_errors_t ne{{n0, test_node_error}, {n1, test_node_error}};

    gate_error_t double_gate_error(0.2);
    op_errors_t test_link_error({{OpType::CX, double_gate_error}});

    op_link_errors_t le{{{n0, n1}, test_link_error}};

    DeviceCharacterisation characterisation(ne, le);

    gate_error_t ne_0_info = characterisation.get_error(n0, OpType::X);
    gate_error_t ne_1_info = characterisation.get_error(n1, OpType::X);

    // gate information preserved correctly
    REQUIRE(ne_0_info == 0.3);
    REQUIRE(ne_1_info == 0.3);

    // link information
    Architecture::Connection test_link = {n0, n1};
    gate_error_t cx_info_ge = characterisation.get_error(test_link, OpType::CX);
    double cx_info = cx_info_ge;
    REQUIRE(cx_info == 0.2);
  }
}

}  // namespace test_DeviceCharacterisation
}  // namespace tket
