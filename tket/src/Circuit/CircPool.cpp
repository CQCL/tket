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

Circuit CRz_using_CX(Expr alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::Rz, alpha / 2, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::Rz, -alpha / 2, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  return c;
}

Circuit CRx_using_CX(Expr alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::Rx, alpha / 2, {1});
  c.add_op<unsigned>(OpType::H, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::H, {1});
  c.add_op<unsigned>(OpType::Rx, -alpha / 2, {1});
  c.add_op<unsigned>(OpType::H, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::H, {1});
  return c;
}

Circuit CRy_using_CX(Expr alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::Ry, alpha / 2, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::Ry, -alpha / 2, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  return c;
}

Circuit CU1_using_CX(Expr lambda) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::U1, lambda / 2, {0});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U1, -lambda / 2, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U1, lambda / 2, {1});
  return c;
}

Circuit CU3_using_CX(Expr theta, Expr phi, Expr lambda) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::U1, (lambda + phi) / 2, {0});
  c.add_op<unsigned>(OpType::U1, (lambda - phi) / 2, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {-theta / 2, 0., -(lambda + phi) / 2}, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {theta / 2, phi, 0.}, {1});
  return c;
}

Circuit ISWAP_using_CX(Expr alpha) {
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

Circuit XXPhase_using_CX(Expr alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::U3, {alpha, -0.5, 0.5}, {0});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  return c;
}

Circuit YYPhase_using_CX(Expr alpha) {
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

Circuit ZZPhase_using_CX(Expr alpha) {
  Circuit c(2);
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::Rz, alpha, {1});
  c.add_op<unsigned>(OpType::CX, {0, 1});
  return c;
}

Circuit XXPhase3_using_CX(Expr alpha) {
  Circuit c(3);
  Circuit rep1 = XXPhase_using_CX(alpha);
  c.append_qubits(rep1, {0, 1});
  c.append_qubits(rep1, {1, 2});
  c.append_qubits(rep1, {0, 2});
  return c;
}

Circuit ESWAP_using_CX(Expr alpha) {
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

Circuit FSim_using_CX(Expr alpha, Expr beta) {
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

Circuit PhasedISWAP_using_CX(Expr p, Expr t) {
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

Circuit NPhasedX_using_CX(
    unsigned int number_of_qubits, Expr alpha, Expr beta) {
  Circuit c(number_of_qubits);
  for (unsigned int i = 0; i < number_of_qubits; ++i) {
    c.add_op<unsigned>(OpType::PhasedX, {alpha, beta}, {i});
  }
  return c;
}

}  // namespace CircPool

}  // namespace tket
