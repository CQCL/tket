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

#include "Transform.hpp"

namespace tket {

namespace Transforms {

class ControlDecompError : public std::logic_error {
 public:
  explicit ControlDecompError(const std::string& message)
      : std::logic_error(message) {}
};

Circuit incrementer_borrow_1_qubit(unsigned n);

Circuit incrementer_borrow_n_qubits(unsigned n);

Circuit cnx_normal_decomp(unsigned n);

Circuit cnx_gray_decomp(unsigned n);

// does not use ancillae
// Expects: CCX + any other gates
// returns CX, H, T, Tdg + any previous gates
Transform decomp_CCX();

Circuit decomposed_CnRy(const Op_ptr op, unsigned arity);

// converts arbitrarily controlled Ry gates. It does not use ancillae, so it
// is not very depth-efficient Expects: CRys and any other gates returns Ry,
// CX, H, T, Tdg + whatever other gates were there before
Transform decomp_controlled_Rys();

// does not use ancillae
// Expects: any CnRys + CnXs + any other gates
// returns Ry, CX, H, T, Tdg + any previous gates
Transform decomp_arbitrary_controlled_gates();

}  // namespace Transforms

}  // namespace tket
