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

#include "CircPool.hpp"

#include "Utils/Expression.hpp"

namespace tket {

namespace CircPool {

const Circuit &BRIDGE_using_CX_0() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(3);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    return c;
  }());
  return *C;
}

const Circuit &BRIDGE_using_CX_1() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(3);
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    return c;
  }());
  return *C;
}

const Circuit &CX_using_TK2() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(tket::OpType::V, {0});
    c.add_op<unsigned>(tket::OpType::S, {0});
    c.add_op<unsigned>(tket::OpType::V, {1});
    c.add_op<unsigned>(tket::OpType::Z, {1});
    c.add_op<unsigned>(tket::OpType::TK2, {-0.5, 0, 0}, {0, 1});
    c.add_op<unsigned>(tket::OpType::H, {0});
    c.add_op<unsigned>(tket::OpType::Y, {1});
    return c;
  }());
  return *C;
}

const Circuit &CX_using_flipped_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(tket::OpType::H, {0});
    c.add_op<unsigned>(tket::OpType::H, {1});
    c.add_op<unsigned>(tket::OpType::CX, {1, 0});
    c.add_op<unsigned>(tket::OpType::H, {0});
    c.add_op<unsigned>(tket::OpType::H, {1});
    return c;
  }());
  return *C;
}

const Circuit &CX_using_ECR() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::U3, {-1, -1, -1.5}, {0});
    c.add_op<unsigned>(OpType::Rx, 0.5, {1});
    c.add_op<unsigned>(OpType::ECR, {0, 1});
    return c;
  }());
  return *C;
}

const Circuit &CX_using_ZZMax() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::Rz, 1.5, {0});
    c.add_op<unsigned>(OpType::Rx, 0.5, {1});
    c.add_op<unsigned>(OpType::Rz, 1.5, {1});
    c.add_op<unsigned>(OpType::Rx, 1.5, {1});
    c.add_op<unsigned>(OpType::ZZMax, {0, 1});
    c.add_op<unsigned>(OpType::Rx, 1.5, {1});
    c.add_op<unsigned>(OpType::Rz, 1.5, {1});
    c.add_phase(0.75);
    return c;
  }());
  return *C;
}

const Circuit &CX_using_XXPhase_0() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::Ry, 0.5, {0});
    c.add_op<unsigned>(OpType::XXPhase, 0.5, {0, 1});
    c.add_op<unsigned>(OpType::Ry, -0.5, {0});
    c.add_op<unsigned>(OpType::Rz, -0.5, {0});
    c.add_op<unsigned>(OpType::Rx, -0.5, {1});
    c.add_phase(-0.25);
    return c;
  }());
  return *C;
}

const Circuit &CX_using_XXPhase_1() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::Rx, 1.5, {1});
    c.add_op<unsigned>(OpType::Rx, 0.5, {0});
    c.add_op<unsigned>(OpType::Rz, 0.5, {0});
    c.add_op<unsigned>(OpType::XXPhase, 0.5, {0, 1});
    c.add_op<unsigned>(OpType::Rz, 0.5, {0});
    c.add_op<unsigned>(OpType::Rx, 0.5, {0});
    c.add_op<unsigned>(OpType::Rz, 0.5, {0});
    c.add_phase(-0.25);
    return c;
  }());
  return *C;
}

const Circuit &CX_VS_CX_reduced() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::Z, {0});
    c.add_op<unsigned>(OpType::X, {1});
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::V, {1});
    c.add_op<unsigned>(OpType::V, {0});
    c.add_op<unsigned>(OpType::S, {1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::V, {1});
    c.add_op<unsigned>(OpType::SWAP, {0, 1});
    c.add_phase(0.5);
    return c;
  }());
  return *C;
}

const Circuit &CX_V_CX_reduced() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::X, {0});
    c.add_op<unsigned>(OpType::V, {0});
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::V, {0});
    c.add_op<unsigned>(OpType::V, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::V, {0});
    c.add_op<unsigned>(OpType::S, {0});
    c.add_phase(0.25);
    return c;
  }());
  return *C;
}

const Circuit &CX_S_CX_reduced() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::Z, {1});
    c.add_op<unsigned>(OpType::S, {1});
    c.add_op<unsigned>(OpType::V, {1});
    c.add_op<unsigned>(OpType::S, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::S, {1});
    c.add_op<unsigned>(OpType::V, {1});
    return c;
  }());
  return *C;
}

const Circuit &CX_V_S_XC_reduced() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::Z, {1});
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::S, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::V, {0});
    c.add_op<unsigned>(OpType::S, {0});
    c.add_op<unsigned>(OpType::S, {1});
    return c;
  }());
  return *C;
}

const Circuit &CX_S_V_XC_reduced() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::X, {0});
    c.add_op<unsigned>(OpType::V, {0});
    c.add_op<unsigned>(OpType::V, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::V, {0});
    c.add_op<unsigned>(OpType::S, {1});
    c.add_op<unsigned>(OpType::V, {1});
    c.add_phase(0.75);
    return c;
  }());
  return *C;
}

const Circuit &CX_XC_reduced() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::SWAP, {0, 1});
    return c;
  }());
  return *C;
}

const Circuit &SWAP_using_CX_0() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    return c;
  }());
  return *C;
}

const Circuit &SWAP_using_CX_1() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    return c;
  }());
  return *C;
}

const Circuit &two_Rz1() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    Op_ptr z = get_op_ptr(OpType::Rz, 1.);
    c.add_op<unsigned>(z, {0});
    c.add_op<unsigned>(z, {1});
    return c;
  }());
  return *C;
}

const Circuit &X1_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::X, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    return c;
  }());
  return *C;
}

const Circuit &Z0_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::Z, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    return c;
  }());
  return *C;
}

const Circuit &CCX_modulo_phase_shift() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(3);
    c.add_op<unsigned>(OpType::Ry, 0.25, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Ry, 0.25, {2});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::Ry, -0.25, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Ry, -0.25, {2});
    return c;
  }());
  return *C;
}

const Circuit &CCX_normal_decomp() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(3);
    c.add_op<unsigned>(OpType::H, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Tdg, {2});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Tdg, {2});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::H, {2});
    c.add_op<unsigned>(OpType::T, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::Tdg, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    return c;
  }());
  return *C;
}

const Circuit &C3X_normal_decomp() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(4);
    c.add_op<unsigned>(OpType::H, {3});
    c.add_op<unsigned>(OpType::U1, 0.125, {0});
    c.add_op<unsigned>(OpType::U1, 0.125, {1});
    c.add_op<unsigned>(OpType::U1, 0.125, {2});
    c.add_op<unsigned>(OpType::U1, 0.125, {3});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::U1, -0.125, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::U1, -0.125, {2});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::U1, 0.125, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::U1, -0.125, {2});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::CX, {2, 3});
    c.add_op<unsigned>(OpType::U1, -0.125, {3});
    c.add_op<unsigned>(OpType::CX, {1, 3});
    c.add_op<unsigned>(OpType::U1, 0.125, {3});
    c.add_op<unsigned>(OpType::CX, {2, 3});
    c.add_op<unsigned>(OpType::U1, -0.125, {3});
    c.add_op<unsigned>(OpType::CX, {0, 3});
    c.add_op<unsigned>(OpType::U1, 0.125, {3});
    c.add_op<unsigned>(OpType::CX, {2, 3});
    c.add_op<unsigned>(OpType::U1, -0.125, {3});
    c.add_op<unsigned>(OpType::CX, {1, 3});
    c.add_op<unsigned>(OpType::U1, 0.125, {3});
    c.add_op<unsigned>(OpType::CX, {2, 3});
    c.add_op<unsigned>(OpType::U1, -0.125, {3});
    c.add_op<unsigned>(OpType::CX, {0, 3});
    c.add_op<unsigned>(OpType::H, {3});
    return c;
  }());
  return *C;
}

// Implementing C3X up to a relative phase
// https://arxiv.org/pdf/1508.03273.pdf figure 4
static const Circuit &RC3X_normal_decomp() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(4);
    c.add_op<unsigned>(OpType::U2, {0, 1}, {3});
    c.add_op<unsigned>(OpType::U1, 0.25, {3});
    c.add_op<unsigned>(OpType::CX, {2, 3});
    c.add_op<unsigned>(OpType::U1, -0.25, {3});
    c.add_op<unsigned>(OpType::U2, {0, 1}, {3});
    c.add_op<unsigned>(OpType::CX, {0, 3});
    c.add_op<unsigned>(OpType::U1, 0.25, {3});
    c.add_op<unsigned>(OpType::CX, {1, 3});
    c.add_op<unsigned>(OpType::U1, -0.25, {3});
    c.add_op<unsigned>(OpType::CX, {0, 3});
    c.add_op<unsigned>(OpType::U1, 0.25, {3});
    c.add_op<unsigned>(OpType::CX, {1, 3});
    c.add_op<unsigned>(OpType::U1, -0.25, {3});
    c.add_op<unsigned>(OpType::U2, {0, 1}, {3});
    c.add_op<unsigned>(OpType::U1, 0.25, {3});
    c.add_op<unsigned>(OpType::CX, {2, 3});
    c.add_op<unsigned>(OpType::U1, -0.25, {3});
    c.add_op<unsigned>(OpType::U2, {0, 1}, {3});
    return c;
  }());
  return *C;
}

// Implementing 3-controlled SX gate
// https://arxiv.org/pdf/quant-ph/9503016.pdf page 17
static const Circuit &C3SX_normal_decomp() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(4);
    c.add_op<unsigned>(OpType::H, {3});
    c.append_qubits(CU1_using_CX(-0.125), {0, 3});
    c.add_op<unsigned>(OpType::H, {3});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {3});
    c.append_qubits(CU1_using_CX(0.125), {1, 3});
    c.add_op<unsigned>(OpType::H, {3});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {3});
    c.append_qubits(CU1_using_CX(-0.125), {1, 3});
    c.add_op<unsigned>(OpType::H, {3});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::H, {3});
    c.append_qubits(CU1_using_CX(0.125), {2, 3});
    c.add_op<unsigned>(OpType::H, {3});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::H, {3});
    c.append_qubits(CU1_using_CX(-0.125), {2, 3});
    c.add_op<unsigned>(OpType::H, {3});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::H, {3});
    c.append_qubits(CU1_using_CX(0.125), {2, 3});
    c.add_op<unsigned>(OpType::H, {3});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::H, {3});
    c.append_qubits(CU1_using_CX(-0.125), {2, 3});
    c.add_op<unsigned>(OpType::H, {3});
    return c;
  }());
  return *C;
}

// https://arxiv.org/pdf/quant-ph/9503016.pdf lemma 7.5
const Circuit &C4X_normal_decomp() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(5);
    c.add_op<unsigned>(OpType::H, {4});
    c.append_qubits(CU1_using_CX(-0.5), {3, 4});
    c.add_op<unsigned>(OpType::H, {4});

    c.append_qubits(RC3X_normal_decomp(), {0, 1, 2, 3});
    c.add_op<unsigned>(OpType::H, {4});
    c.append_qubits(CU1_using_CX(0.5), {3, 4});
    c.add_op<unsigned>(OpType::H, {4});
    c.append_qubits(RC3X_normal_decomp().dagger(), {0, 1, 2, 3});
    c.append_qubits(C3SX_normal_decomp(), {0, 1, 2, 4});
    return c;
  }());
  return *C;
}

const Circuit &ladder_down() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(3);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    return c;
  }());
  return *C;
}

const Circuit &ladder_down_2() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(3);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::X, {0});
    c.add_op<unsigned>(OpType::X, {2});
    c.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    return c;
  }());
  return *C;
}

const Circuit &ladder_up() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(3);
    c.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::CX, {2, 1});
    return c;
  }());
  return *C;
}

const Circuit &X() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(1);
    c.add_op<unsigned>(OpType::X, {0});
    return c;
  }());
  return *C;
}

const Circuit &CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    return c;
  }());
  return *C;
}

const Circuit &CCX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(3);
    c.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    return c;
  }());
  return *C;
}

const Circuit &BRIDGE() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(3);
    c.add_op<unsigned>(OpType::BRIDGE, {0, 1, 2});
    return c;
  }());
  return *C;
}

const Circuit &H_CZ_H() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::CZ, {0, 1});
    c.add_op<unsigned>(OpType::H, {1});
    return c;
  }());
  return *C;
}

const Circuit &CZ_using_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {1});
    return c;
  }());
  return *C;
}

const Circuit &CY_using_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::Sdg, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::S, {1});
    return c;
  }());
  return *C;
}

const Circuit &CH_using_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::Sdg, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::T, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::T, {1});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::S, {1});
    c.add_op<unsigned>(OpType::X, {1});
    c.add_op<unsigned>(OpType::S, {0});
    c.add_phase(-0.25);
    return c;
  }());
  return *C;
}

const Circuit &CV_using_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c = CRx_using_CX(0.5);
    return c;
  }());
  return *C;
}

const Circuit &CVdg_using_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c = CRx_using_CX(-0.5);
    return c;
  }());
  return *C;
}

const Circuit &CSX_using_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {1});
    Circuit c2 = CU1_using_CX(0.5);
    c.append(c2);
    c.add_op<unsigned>(OpType::H, {1});
    return c;
  }());
  return *C;
}

const Circuit &CSXdg_using_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {1});
    Circuit c2 = CU1_using_CX(-0.5);
    c.append(c2);
    c.add_op<unsigned>(OpType::H, {1});
    return c;
  }());
  return *C;
}

const Circuit &CSWAP_using_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(3);
    c.add_op<unsigned>(OpType::CX, {2, 1});
    c.add_op<unsigned>(OpType::H, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Tdg, {2});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::Tdg, {2});
    c.add_op<unsigned>(OpType::CX, {0, 2});
    c.add_op<unsigned>(OpType::T, {1});
    c.add_op<unsigned>(OpType::T, {2});
    c.add_op<unsigned>(OpType::H, {2});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::Tdg, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {2, 1});
    return c;
  }());
  return *C;
}

const Circuit &ECR_using_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::Rx, -0.5, {1});
    c.add_op<unsigned>(OpType::U3, {1, 1.5, 1}, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    return c;
  }());
  return *C;
}

const Circuit &ZZMax_using_CX() {
  static std::unique_ptr<const Circuit> C = std::make_unique<Circuit>([]() {
    Circuit c(2);
    c.add_op<unsigned>(OpType::Rz, 0.5, {0});
    c.add_op<unsigned>(OpType::U3, {0.5, 0, 0}, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::U3, {0.5, -0.5, 1}, {1});
    return c;
  }());
  return *C;
}

Circuit CRz_using_TK2(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK1, {0.5, 0.5, 1}, {0});
  c.add_op<unsigned>(OpType::TK1, {0.5, 0.5, 0}, {1});
  c.add_op<unsigned>(OpType::TK2, {-0.5 * alpha, 0, 0}, {0, 1});
  c.add_op<unsigned>(OpType::TK1, {0, 0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::TK1, {-1 + 0.5 * alpha, 0.5, 0.5}, {1});
  c.add_phase(1);
  return c;
}

Circuit CRz_using_CX(const Expr &alpha) {
  Circuit c(2);
  if (equiv_expr(alpha, 1.)) {
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {1});
    if (equiv_expr(alpha, 1., 4)) {
      c.add_op<unsigned>(OpType::Sdg, {0});
    } else {
      c.add_op<unsigned>(OpType::S, {0});
    }
  } else {
    c.add_op<unsigned>(OpType::Rz, alpha / 2, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::Rz, -alpha / 2, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
  }
  return c;
}

Circuit CRx_using_TK2(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK1, {-0.5, 0.5, 0}, {0});
  c.add_op<unsigned>(OpType::TK1, {1, 0.5, 0}, {1});
  c.add_op<unsigned>(OpType::TK2, {-0.5 * alpha, 0, 0}, {0, 1});
  c.add_op<unsigned>(OpType::TK1, {-1, 0.5, -0.5}, {0});
  c.add_op<unsigned>(OpType::TK1, {1, 0.5 - 0.5 * alpha, 0}, {1});
  return c;
}

Circuit CRx_using_CX(const Expr &alpha) {
  Circuit c(2);
  if (equiv_expr(alpha, 1.)) {
    c.add_op<unsigned>(OpType::CX, {0, 1});
    if (equiv_expr(alpha, 1., 4)) {
      c.add_op<unsigned>(OpType::Sdg, {0});
    } else {
      c.add_op<unsigned>(OpType::S, {0});
    }
  } else {
    c.add_op<unsigned>(OpType::Rx, alpha / 2, {1});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::Rx, -alpha / 2, {1});
    c.add_op<unsigned>(OpType::H, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::H, {1});
  }
  return c;
}

Circuit CRy_using_TK2(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK1, {0.5, 0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::TK1, {0, 0.5, -0.5}, {1});
  c.add_op<unsigned>(OpType::TK2, {-0.5 * alpha, 0, 0}, {0, 1});
  c.add_op<unsigned>(OpType::TK1, {0.5, 0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::TK1, {-0.5, 0.5 - 0.5 * alpha, 1}, {1});
  c.add_phase(-1);
  return c;
}

Circuit CRy_using_CX(const Expr &alpha) {
  Circuit c(2);
  if (equiv_expr(alpha, 1.)) {
    c.add_op<unsigned>(OpType::Sdg, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::S, {1});
    if (equiv_expr(alpha, 1., 4)) {
      c.add_op<unsigned>(OpType::Sdg, {0});
    } else {
      c.add_op<unsigned>(OpType::S, {0});
    }
  } else {
    c.add_op<unsigned>(OpType::Ry, alpha / 2, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::Ry, -alpha / 2, {1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
  }
  return c;
}

Circuit CU1_using_TK2(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK1, {0.5, 0.5, 1}, {0});
  c.add_op<unsigned>(OpType::TK1, {0.5, 0.5, 0}, {1});
  c.add_op<unsigned>(OpType::TK2, {-0.5 * alpha, 0, 0}, {0, 1});
  c.add_op<unsigned>(OpType::TK1, {0.5 * alpha, 0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::TK1, {-1 + 0.5 * alpha, 0.5, 0.5}, {1});
  c.add_phase(-1 + 0.25 * alpha);
  return c;
}

Circuit CU1_using_CX(const Expr &lambda) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::U1, lambda / 2, {0});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U1, -lambda / 2, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U1, lambda / 2, {1});
  return c;
}

Circuit CU3_using_CX(const Expr &theta, const Expr &phi, const Expr &lambda) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::U1, (lambda + phi) / 2, {0});
  c.add_op<unsigned>(OpType::U1, (lambda - phi) / 2, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {-theta / 2, 0., -(lambda + phi) / 2}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {theta / 2, phi, 0.}, {1});
  return c;
}

Circuit ISWAP_using_TK2(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK2, {-0.5 * alpha, -0.5 * alpha, 0}, {0, 1});
  return c;
}

Circuit ISWAP_using_CX(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::U3, {0.5, -0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::U3, {0.5, -0.5, 0.5}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {-0.5 * alpha, -0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::Rz, -0.5 * alpha, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {-0.5, -0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::U3, {-0.5, -0.5, 0.5}, {1});
  return c;
}

Circuit XXPhase_using_TK2(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK2, {alpha, 0, 0}, {0, 1});
  return c;
}

Circuit XXPhase_using_CX(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {alpha, -0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  return c;
}

Circuit YYPhase_using_TK2(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK2, {0, alpha, 0}, {0, 1});
  return c;
}

Circuit YYPhase_using_CX(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::U3, {0.5, -0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::U3, {0.5, -0.5, 0.5}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::Rz, alpha, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {-0.5, -0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::U3, {-0.5, -0.5, 0.5}, {1});
  return c;
}

Circuit ZZPhase_using_TK2(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK2, {0, 0, alpha}, {0, 1});
  return c;
}

Circuit ZZPhase_using_CX(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::Rz, alpha, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  return c;
}

Circuit XXPhase_using_ZZPhase(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::H, {0});
  c.add_op<unsigned>(OpType::H, {1});
  c.add_op<unsigned>(OpType::ZZPhase, alpha, {0, 1});
  c.add_op<unsigned>(OpType::H, {0});
  c.add_op<unsigned>(OpType::H, {1});
  return c;
}

Circuit YYPhase_using_ZZPhase(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::Vdg, {0});
  c.add_op<unsigned>(OpType::Vdg, {1});
  c.add_op<unsigned>(OpType::ZZPhase, alpha, {0, 1});
  c.add_op<unsigned>(OpType::V, {0});
  c.add_op<unsigned>(OpType::V, {1});
  return c;
}

Circuit approx_TK2_using_1xCX() {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK1, {0.5, 1.5, 1.5}, {0});
  c.add_op<unsigned>(OpType::TK1, {0., 0.5, 0}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::TK1, {0.5, 0.5, 0}, {0});
  c.add_phase(0.25);
  return c;
}

Circuit approx_TK2_using_2xCX(const Expr &alpha, const Expr &beta) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK1, {0.5, 0.5, 0}, {0});
  c.add_op<unsigned>(OpType::TK1, {0., 1.5, 0}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::TK1, {0., 1 + alpha, 1.5}, {0});
  c.add_op<unsigned>(OpType::TK1, {0., 1.5, 2 - beta}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::TK1, {0., 0.5, 0}, {0});
  c.add_phase(0.5);
  return c;
}

Circuit TK2_using_3xCX(const Expr &alpha, const Expr &beta, const Expr &gamma) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK1, {0.5, 1.5, 1}, {0});
  c.add_op<unsigned>(OpType::TK1, {0, 0.5, 0}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::TK1, {3.5 + alpha, 3.5, 0.}, {0});
  c.add_op<unsigned>(OpType::TK1, {0.5, 1, 0.5 + beta}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::TK1, {0.5, 0.5, 0.}, {0});
  c.add_op<unsigned>(OpType::TK1, {0, 0, gamma}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_phase(0.75);
  return c;
}

Circuit TK2_using_CX(const Expr &alpha, const Expr &beta, const Expr &gamma) {
  // only handle TK2 if normalised to Weyl chamber
  if (equiv_0(alpha, 4) && equiv_0(beta, 4) && equiv_0(gamma, 4)) {
    return Circuit(2);
  } else if (
      equiv_expr(alpha, 0.5, 4) && equiv_0(beta, 4) && equiv_0(gamma, 4)) {
    return approx_TK2_using_1xCX();
  } else if (equiv_0(gamma, 4)) {
    return approx_TK2_using_2xCX(alpha, beta);
  } else {
    return TK2_using_3xCX(alpha, beta, gamma);
  }
}

Circuit approx_TK2_using_1xZZPhase(const Expr &alpha) {
  return XXPhase_using_ZZPhase(alpha);
}

Circuit approx_TK2_using_2xZZPhase(const Expr &alpha, const Expr &beta) {
  Circuit c(2);
  c.append(XXPhase_using_ZZPhase(alpha));
  c.append(YYPhase_using_ZZPhase(beta));
  return c;
}

Circuit TK2_using_ZZPhase(
    const Expr &alpha, const Expr &beta, const Expr &gamma) {
  Circuit c(2);
  c.append(XXPhase_using_ZZPhase(alpha));
  c.append(YYPhase_using_ZZPhase(beta));
  c.add_op<unsigned>(OpType::ZZPhase, gamma, {0, 1});
  return c;
}

Circuit XXPhase3_using_TK2(const Expr &alpha) {
  Circuit c(3);
  c.add_op<unsigned>(OpType::TK2, {alpha, 0, 0}, {0, 1});
  c.add_op<unsigned>(OpType::TK2, {alpha, 0, 0}, {1, 2});
  c.add_op<unsigned>(OpType::TK2, {alpha, 0, 0}, {0, 2});
  return c;
}

Circuit XXPhase3_using_CX(const Expr &alpha) {
  Circuit c(3);
  Circuit rep1 = XXPhase_using_CX(alpha);
  c.append_qubits(rep1, {0, 1});
  c.append_qubits(rep1, {1, 2});
  c.append_qubits(rep1, {0, 2});
  return c;
}

Circuit ESWAP_using_TK2(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK1, {0.5, 0.5, -0.5}, {0});
  c.add_op<unsigned>(OpType::TK1, {-0.5, 0.5, -0.5}, {1});
  c.add_op<unsigned>(
      OpType::TK2, {-0.5 * alpha, -0.5 * alpha, 0.5 * alpha}, {0, 1});
  c.add_op<unsigned>(OpType::TK1, {-0.5, 0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::TK1, {-0.5, 0.5, -0.5}, {1});
  c.add_phase(1 - 0.25 * alpha);
  return c;
}

Circuit ESWAP_using_CX(const Expr &alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::S, {0});
  c.add_op<unsigned>(OpType::X, {1});
  c.add_op<unsigned>(OpType::CX, {1, 0});
  c.add_op<unsigned>(OpType::U1, 0.5 - 0.5 * alpha, {0});
  c.add_op<unsigned>(OpType::Ry, -0.5 + 0.5 * alpha, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::Ry, 0.5 + 0.5 * alpha, {1});
  c.add_op<unsigned>(OpType::CX, {1, 0});
  c.add_op<unsigned>(OpType::X, {1});
  c.add_op<unsigned>(OpType::S, {1});
  c.add_phase(-0.5);
  return c;
}

Circuit FSim_using_TK2(const Expr &alpha, const Expr &beta) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::TK1, {-0.5, 0.5, -1}, {0});
  c.add_op<unsigned>(OpType::TK1, {0.5, 0.5, 1}, {1});
  c.add_op<unsigned>(OpType::TK2, {-0.5 * beta, -alpha, alpha}, {0, 1});
  c.add_op<unsigned>(OpType::TK1, {-0.5 * beta, 0.5, -0.5}, {0});
  c.add_op<unsigned>(OpType::TK1, {-0.5 * beta, 0.5, 0.5}, {1});
  c.add_phase(-0.25 * beta);
  return c;
}

Circuit FSim_using_CX(const Expr &alpha, const Expr &beta) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::U3, {0.5, 0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::U3, {0.5, 0., 1.5}, {1});
  c.add_op<unsigned>(OpType::CX, {1, 0});
  c.add_op<unsigned>(OpType::U1, 0.5 - alpha, {0});
  c.add_op<unsigned>(OpType::U3, {-0.5 + alpha, 0, 0}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {0.5 + 0.5 * beta, 0, 0}, {1});
  c.add_op<unsigned>(OpType::CX, {1, 0});
  c.add_op<unsigned>(OpType::U3, {0.5, 0.5 - 0.5 * beta, 1}, {0});
  c.add_op<unsigned>(OpType::U3, {0.5, -0.5 - 0.5 * beta, 0.5}, {1});
  c.add_phase(0.5 * alpha + 0.25 * beta);
  return c;
}

Circuit PhasedISWAP_using_TK2(const Expr &p, const Expr &t) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::Rz, p, {0});
  c.add_op<unsigned>(OpType::Rz, -p, {1});
  c.add_op<unsigned>(OpType::TK2, {-0.5 * t, -0.5 * t, 0}, {0, 1});
  c.add_op<unsigned>(OpType::Rz, -p, {0});
  c.add_op<unsigned>(OpType::Rz, p, {1});
  return c;
}

Circuit PhasedISWAP_using_CX(const Expr &p, const Expr &t) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::U3, {0.5, -0.5, 0.5 + p}, {0});
  c.add_op<unsigned>(OpType::U3, {0.5, -0.5, 0.5 - p}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {-0.5 * t, -0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::Rz, -0.5 * t, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {-0.5, -0.5 - p, 0.5}, {0});
  c.add_op<unsigned>(OpType::U3, {-0.5, -0.5 + p, 0.5}, {1});
  return c;
}

Circuit NPhasedX_using_PhasedX(
    unsigned int number_of_qubits, const Expr &alpha, const Expr &beta) {
  Circuit c(number_of_qubits);
  for (unsigned int i = 0; i < number_of_qubits; ++i) {
    c.add_op<unsigned>(OpType::PhasedX, {alpha, beta}, {i});
  }
  return c;
}

static unsigned int_half(const Expr &angle) {
  // Assume angle is an even integer
  double eval = eval_expr(angle).value();
  return lround(eval / 2);
}

Circuit tk1_to_rzsx(const Expr &alpha, const Expr &beta, const Expr &gamma) {
  Circuit c(1);
  Expr correction_phase = 0;
  if (equiv_0(beta)) {
    // b = 2k, if k is odd, then Rx(b) = -I
    c.add_op<unsigned>(OpType::Rz, alpha + gamma, {0});
    correction_phase = int_half(beta);
  } else if (equiv_0(beta + 1)) {
    // Use Rx(2k-1) = i(-1)^{k}SxSx
    correction_phase = -0.5 + int_half(beta - 1);
    if (equiv_0(alpha - gamma)) {
      // a - c = 2m
      // overall operation is (-1)^{m}Rx(2k -1)
      c.add_op<unsigned>(OpType::SX, {0});
      c.add_op<unsigned>(OpType::SX, {0});
      correction_phase += int_half(alpha - gamma);
    } else {
      c.add_op<unsigned>(OpType::Rz, gamma, {0});
      c.add_op<unsigned>(OpType::SX, {0});
      c.add_op<unsigned>(OpType::SX, {0});
      c.add_op<unsigned>(OpType::Rz, alpha, {0});
    }
  } else if (equiv_0(beta - 0.5) && equiv_0(alpha) && equiv_0(gamma)) {
    // a = 2k, b = 2m+0.5, c = 2n
    // Rz(2k)Rx(2m + 0.5)Rz(2n) = (-1)^{k+m+n}e^{-i \pi /4} SX
    c.add_op<unsigned>(OpType::SX, {0});
    correction_phase =
        int_half(beta - 0.5) + int_half(alpha) + int_half(gamma) - 0.25;
  } else if (equiv_0(beta - 0.5)) {
    // SX.Rz(2m-0.5).SX = (-1)^{m}e^{i \pi /4} Rz(-0.5).SX.Rz(-0.5)
    c.add_op<unsigned>(OpType::Rz, gamma, {0});
    c.add_op<unsigned>(OpType::SX, {0});
    c.add_op<unsigned>(OpType::Rz, alpha, {0});
    correction_phase = int_half(beta - 0.5) - 0.25;
  } else if (equiv_0(beta + 0.5)) {
    // SX.Rz(2m+0.5).SX = (-1)^{m}e^{i \pi /4} Rz(0.5).SX.Rz(0.5)
    c.add_op<unsigned>(OpType::Rz, gamma + 1, {0});
    c.add_op<unsigned>(OpType::SX, {0});
    c.add_op<unsigned>(OpType::Rz, alpha + 1, {0});
    correction_phase = int_half(beta - 1.5) - 0.25;
  } else if (equiv_0(alpha - 0.5) && equiv_0(gamma - 0.5)) {
    // Rz(2k + 0.5)Rx(b)Rz(2m + 0.5) = -i(-1)^{k+m}SX.Rz(1-b).SX
    c.add_op<unsigned>(OpType::SX, {0});
    c.add_op<unsigned>(OpType::Rz, 1 - beta, {0});
    c.add_op<unsigned>(OpType::SX, {0});
    correction_phase = int_half(alpha - 0.5) + int_half(gamma - 0.5) - 0.5;
  } else {
    c.add_op<unsigned>(OpType::Rz, gamma + 0.5, {0});
    c.add_op<unsigned>(OpType::SX, {0});
    c.add_op<unsigned>(OpType::Rz, beta - 1, {0});
    c.add_op<unsigned>(OpType::SX, {0});
    c.add_op<unsigned>(OpType::Rz, alpha + 0.5, {0});
    correction_phase = -0.5;
  }
  c.add_phase(correction_phase);
  return c;
}

Circuit tk1_to_rzh(const Expr &alpha, const Expr &beta, const Expr &gamma) {
  Circuit c(1);
  std::optional<unsigned> cliff = equiv_Clifford(beta, 4);
  if (cliff) {
    switch (*cliff % 4) {
      case 0: {
        c.add_op<unsigned>(OpType::Rz, gamma + alpha, {0});
        break;
      }
      case 1: {
        c.add_op<unsigned>(OpType::Rz, gamma - 0.5, {0});
        c.add_op<unsigned>(OpType::H, {0});
        c.add_op<unsigned>(OpType::Rz, alpha - 0.5, {0});
        c.add_phase(-0.5);
        break;
      }
      case 2: {
        c.add_op<unsigned>(OpType::Rz, gamma - alpha, {0});
        c.add_op<unsigned>(OpType::H, {0});
        c.add_op<unsigned>(OpType::Rz, 1., {0});
        c.add_op<unsigned>(OpType::H, {0});
        break;
      }
      case 3: {
        c.add_op<unsigned>(OpType::Rz, gamma + 0.5, {0});
        c.add_op<unsigned>(OpType::H, {0});
        c.add_op<unsigned>(OpType::Rz, alpha + 0.5, {0});
        c.add_phase(-0.5);
        break;
      }
    }
    if (cliff >= 4u) c.add_phase(1.);
  } else {
    c.add_op<unsigned>(OpType::Rz, gamma, {0});
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Rz, beta, {0});
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Rz, alpha, {0});
  }
  return c;
}

Circuit tk1_to_tk1(const Expr &alpha, const Expr &beta, const Expr &gamma) {
  Circuit c(1);
  c.add_op<unsigned>(OpType::TK1, {alpha, beta, gamma}, {0});
  return c;
}

Circuit tk1_to_rzrx(const Expr &alpha, const Expr &beta, const Expr &gamma) {
  Circuit c(1);
  c.add_op<unsigned>(OpType::Rz, gamma, {0});
  c.add_op<unsigned>(OpType::Rx, beta, {0});
  c.add_op<unsigned>(OpType::Rz, alpha, {0});
  return c;
}

Circuit tk1_to_PhasedXRz(
    const Expr &alpha, const Expr &beta, const Expr &gamma) {
  Circuit c(1);
  if (equiv_expr(beta, 1)) {
    // Angles β ∈ {π, 3π}
    c.add_op<unsigned>(OpType::PhasedX, {beta, (alpha - gamma) / 2.}, {0});
  } else if (equiv_expr(beta, 0)) {
    // Angle β ∈ {0, 2π}
    c.add_op<unsigned>(OpType::Rz, alpha + beta + gamma, {0});
  } else {
    c.add_op<unsigned>(OpType::Rz, alpha + gamma, {0});
    c.add_op<unsigned>(OpType::PhasedX, {beta, alpha}, {0});
  }
  return c;
}

}  // namespace CircPool

}  // namespace tket
