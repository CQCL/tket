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

#include "GateUnitaryMatrixImplementations.hpp"
#include "GateUnitaryMatrixUtils.hpp"
#include "Utils/Constants.hpp"

namespace tket {
namespace internal {

namespace {

struct FixedData {
  Eigen::Matrix2cd X;
  Eigen::Matrix2cd Y;
  Eigen::Matrix2cd Z;
  Eigen::Matrix2cd S;
  Eigen::Matrix2cd Sdg;
  Eigen::Matrix2cd T;
  Eigen::Matrix2cd Tdg;
  Eigen::Matrix2cd V;
  Eigen::Matrix2cd Vdg;
  Eigen::Matrix2cd H;
  Eigen::Matrix2cd SX;
  Eigen::Matrix2cd SXdg;
  Eigen::Matrix4cd CSX;
  Eigen::Matrix4cd CSXdg;
  Eigen::Matrix4cd CX;
  Eigen::Matrix4cd CY;
  Eigen::Matrix4cd CZ;
  Eigen::Matrix4cd CH;
  Eigen::Matrix4cd CV;
  Eigen::Matrix4cd CVdg;
  Matrix8cd CCX;
  Eigen::Matrix4cd SWAP;
  Matrix8cd CSWAP;
  Matrix8cd BRIDGE;
  Eigen::Matrix2cd noop;
  Eigen::Matrix4cd ECR;
  Eigen::Matrix4cd ZZMax;
  Eigen::Matrix4cd Sycamore;
  Eigen::Matrix4cd ISWAPMax;
  std::array<unsigned, 8> bridge_columns;
  std::array<unsigned, 8> cswap_columns;

  FixedData() {
    X << 0, 1, 1, 0;
    Y << 0, -i_, i_, 0;
    Z << 1, 0, 0, -1;

    S << 1, 0, 0, i_;
    Sdg = S.adjoint();

    T << 1, 0, 0, std::polar(1.0, 0.25 * PI);
    Tdg = T.adjoint();

    V << 1, -i_, -i_, 1;
    V *= std::sqrt(0.5);
    Vdg = V.adjoint();

    H << 1, 1, 1, -1;
    H *= std::sqrt(0.5);

    const Complex one_plus_i = i_ + 1.0;
    const Complex one_minus_i = -i_ + 1.0;
    SX << one_plus_i, one_minus_i, one_minus_i, one_plus_i;
    SX *= 0.5;
    SXdg = SX.adjoint();

    SWAP << 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1;

    bridge_columns = {0, 1, 2, 3, 5, 4, 7, 6};
    cswap_columns = {0, 1, 2, 3, 4, 6, 5, 7};

    BRIDGE = Matrix8cd::Zero();
    CSWAP = BRIDGE;
    for (unsigned ii = 0; ii < 8; ++ii) {
      BRIDGE(ii, bridge_columns[ii]) = 1.0;
      CSWAP(ii, cswap_columns[ii]) = 1.0;
    }
    noop = Eigen::Matrix2cd::Identity();

    ECR << 0, 0, 1, i_, 0, 0, i_, 1, 1, -i_, 0, 0, -i_, 1, 0, 0;
    ECR *= std::sqrt(0.5);

    CX = GateUnitaryMatrixUtils::get_controlled_gate_unitary(X);

    CCX =
        GateUnitaryMatrixUtils::get_multi_controlled_gate_dense_unitary(CX, 3);

    CY = GateUnitaryMatrixUtils::get_controlled_gate_unitary(Y);
    CZ = GateUnitaryMatrixUtils::get_controlled_gate_unitary(Z);
    CH = GateUnitaryMatrixUtils::get_controlled_gate_unitary(H);
    CV = GateUnitaryMatrixUtils::get_controlled_gate_unitary(V);
    CVdg = GateUnitaryMatrixUtils::get_controlled_gate_unitary(Vdg);
    CSX = GateUnitaryMatrixUtils::get_controlled_gate_unitary(SX);
    CSXdg = GateUnitaryMatrixUtils::get_controlled_gate_unitary(SXdg);

    // Accuracy notes: std::sqrt is guaranteed by IEEE 754
    // but std::sin, std::cos are not (although C++ is not guaranteed
    // to follow IEEE 754 anyway), so where there is a choice
    // (e.g., cos(pi/4) vs. sqrt(2), cos(pi/6) vs. sqrt(3), etc.),
    // std::sqrt might be very slightly more accurate,
    // but it's not worth worrying about.
    // However, note that cos(pi/2), sin(pi), etc. need NOT be exact
    // so we should prefer exact values if it's easy to do.
    ZZMax = GateUnitaryMatrixImplementations::ZZPhase(0.5);
    Sycamore = GateUnitaryMatrixImplementations::FSim(0.5, 1.0 / 6.0);

    // Equals ISWAP(1) in theory, but ensure the exact result!
    ISWAPMax << 1, 0, 0, 0, 0, 0, i_, 0, 0, i_, 0, 0, 0, 0, 0, 1;
  }
};  // struct FixedData

const FixedData& get_fixed_data() {
  // Remember "Static Initialization Order Fiasco"
  static const FixedData data;
  return data;
}

}  // unnamed namespace

#ifdef GATE_FUNCTION_1Q
#error "Macro already defined!"
#endif
#define GATE_FUNCTION_1Q(funcname)                                       \
  const Eigen::Matrix2cd& GateUnitaryMatrixImplementations::funcname() { \
    return get_fixed_data().funcname;                                    \
  }
GATE_FUNCTION_1Q(X)
GATE_FUNCTION_1Q(Y)
GATE_FUNCTION_1Q(Z)
GATE_FUNCTION_1Q(S)
GATE_FUNCTION_1Q(Sdg)
GATE_FUNCTION_1Q(T)
GATE_FUNCTION_1Q(Tdg)
GATE_FUNCTION_1Q(V)
GATE_FUNCTION_1Q(Vdg)
GATE_FUNCTION_1Q(H)
GATE_FUNCTION_1Q(SX)
GATE_FUNCTION_1Q(SXdg)
GATE_FUNCTION_1Q(noop)
#undef GATE_FUNCTION_1Q

#ifdef GATE_FUNCTION_2Q
#error "Macro already defined!"
#endif
#define GATE_FUNCTION_2Q(funcname)                                       \
  const Eigen::Matrix4cd& GateUnitaryMatrixImplementations::funcname() { \
    return get_fixed_data().funcname;                                    \
  }
GATE_FUNCTION_2Q(SWAP)
GATE_FUNCTION_2Q(ECR)
GATE_FUNCTION_2Q(CX)
GATE_FUNCTION_2Q(CY)
GATE_FUNCTION_2Q(CZ)
GATE_FUNCTION_2Q(CH)
GATE_FUNCTION_2Q(CV)
GATE_FUNCTION_2Q(CVdg)
GATE_FUNCTION_2Q(CSX)
GATE_FUNCTION_2Q(CSXdg)
GATE_FUNCTION_2Q(ZZMax)
GATE_FUNCTION_2Q(Sycamore)
GATE_FUNCTION_2Q(ISWAPMax)
#undef GATE_FUNCTION_2Q

#ifdef GATE_FUNCTION_XQ
#error "Macro already defined!"
#endif
#define GATE_FUNCTION_3Q(funcname)                                \
  const Matrix8cd& GateUnitaryMatrixImplementations::funcname() { \
    return get_fixed_data().funcname;                             \
  }
GATE_FUNCTION_3Q(BRIDGE)
GATE_FUNCTION_3Q(CCX)
GATE_FUNCTION_3Q(CSWAP)
#undef GATE_FUNCTION_3Q

const std::array<unsigned, 8>&
GateUnitaryMatrixImplementations::get_bridge_columns() {
  return get_fixed_data().bridge_columns;
}

const std::array<unsigned, 8>&
GateUnitaryMatrixImplementations::get_cswap_columns() {
  return get_fixed_data().cswap_columns;
}

}  // namespace internal
}  // namespace tket
