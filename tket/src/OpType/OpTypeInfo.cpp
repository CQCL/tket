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

#include "OpTypeInfo.hpp"

#include <memory>

#include "OpType.hpp"

namespace tket {
const std::map<OpType, OpTypeInfo>& optypeinfo() {
  static const op_signature_t noargs;
  static const op_signature_t singleq(1, EdgeType::Quantum);
  static const op_signature_t doubleq(2, EdgeType::Quantum);
  static const op_signature_t tripleq(3, EdgeType::Quantum);

  // No memory leak since the objects are only constructed once and
  // the memory is reclaimed on program termination.

  // OpTypeInfo: name, latex_name, param_mod, signature
  static const std::map<OpType, OpTypeInfo> typeinfo{
      {OpType::Z, {"Z", "$Z$", {}, singleq}},
      {OpType::X, {"X", "$X$", {}, singleq}},
      {OpType::Y, {"Y", "$Y$", {}, singleq}},
      {OpType::S, {"S", "$S$", {}, singleq}},
      {OpType::Sdg, {"Sdg", "$S^\\dagger$", {}, singleq}},
      {OpType::T, {"T", "$T$", {}, singleq}},
      {OpType::Tdg, {"Tdg", "$T^\\dagger$", {}, singleq}},
      {OpType::V, {"V", "$R_X{\\frac12}$", {}, singleq}},
      {OpType::Vdg, {"Vdg", "$R_X(\\frac12)^\\dagger$", {}, singleq}},
      {OpType::SX, {"SX", "$\\sqrt{X}$", {}, singleq}},
      {OpType::SXdg, {"SXdg", "$\\sqrt{X}^\\dagger$", {}, singleq}},
      {OpType::H, {"H", "$H$", {}, singleq}},
      {OpType::Rx, {"Rx", "$R_X$", {4}, singleq}},
      {OpType::Ry, {"Ry", "$R_Y$", {4}, singleq}},
      {OpType::Rz, {"Rz", "$R_Z$", {4}, singleq}},
      {OpType::U3, {"U3", "U3", {4, 2, 2}, singleq}},
      {OpType::U2, {"U2", "U2", {2, 2}, singleq}},
      {OpType::U1, {"U1", "U1", {2}, singleq}},
      {OpType::CX, {"CX", "CX", {}, doubleq}},
      {OpType::CY, {"CY", "CY", {}, doubleq}},
      {OpType::CZ, {"CZ", "CZ", {}, doubleq}},
      {OpType::CH, {"CH", "CH", {}, doubleq}},
      {OpType::CV, {"CV", "CV", {}, doubleq}},
      {OpType::CVdg, {"CVdg", "$CV^\\dagger$", {}, doubleq}},
      {OpType::CSX, {"CSX", "CSX", {}, doubleq}},
      {OpType::CSXdg, {"CSXdg", "$CSX^\\dagger$", {}, doubleq}},
      {OpType::CRz, {"CRz", "CRz", {4}, doubleq}},
      {OpType::CRx, {"CRx", "CRx", {4}, doubleq}},
      {OpType::CRy, {"CRy", "CRy", {4}, doubleq}},
      {OpType::CU1, {"CU1", "CU1", {2}, doubleq}},
      {OpType::CU3, {"CU3", "CU3", {4, 2, 2}, doubleq}},
      {OpType::PhaseGadget,
       {"PhaseGadget", "$Z^{\\otimes n}$", {4}, std::nullopt}},
      {OpType::CCX, {"CCX", "CCX", {}, tripleq}},
      {OpType::SWAP, {"SWAP", "SWAP", {}, doubleq}},
      {OpType::CSWAP, {"CSWAP", "CSWAP", {}, tripleq}},
      {OpType::BRIDGE, {"BRIDGE", "BRIDGE", {}, tripleq}},
      {OpType::Input, {"Input", "Q IN", {}, singleq}},
      {OpType::Output, {"Output", "Q OUT", {}, singleq}},
      {OpType::Create, {"Create", "Q CREATE", {}, singleq}},
      {OpType::Discard, {"Discard", "Q DISCARD", {}, singleq}},
      {OpType::ClInput,
       {"ClInput", "C IN", {}, op_signature_t({EdgeType::Classical})}},
      {OpType::ClOutput,
       {"ClOutput", "C OUT", {}, op_signature_t({EdgeType::Classical})}},
      {OpType::Label, {"Label", "Label", {}, noargs}},
      {OpType::Branch,
       {"Branch", "Branch", {}, op_signature_t({EdgeType::Boolean})}},
      {OpType::Goto, {"Goto", "Goto", {}, noargs}},
      {OpType::Stop, {"Stop", "Stop", {}, noargs}},
      {OpType::noop, {"noop", "-", {}, singleq}},
      {OpType::CircBox, {"CircBox", "CircBox", {}, std::nullopt}},
      {OpType::Unitary1qBox, {"Unitary1qBox", "Unitary1qBox", {}, singleq}},
      {OpType::Unitary2qBox, {"Unitary2qBox", "Unitary2qBox", {}, doubleq}},
      {OpType::Unitary3qBox, {"Unitary3qBox", "Unitary3qBox", {}, tripleq}},
      {OpType::ExpBox, {"ExpBox", "ExpBox", {}, doubleq}},
      {OpType::PauliExpBox, {"PauliExpBox", "PauliExpBox", {}, std::nullopt}},
      {OpType::CustomGate, {"CustomGate", "CustomGate", {}, std::nullopt}},
      {OpType::Barrier, {"Barrier", "Barrier", {}, std::nullopt}},
      {OpType::Measure,
       {"Measure",
        "Measure",
        {},
        op_signature_t({EdgeType::Quantum, EdgeType::Classical})}},
      {OpType::Collapse, {"Collapse", "Collapse", {}, singleq}},
      {OpType::Reset, {"Reset", "Reset", {}, singleq}},
      {OpType::ECR, {"ECR", "ECR", {}, doubleq}},
      {OpType::ISWAP, {"ISWAP", "ISWAP", {4}, doubleq}},
      {OpType::PhasedX, {"PhasedX", "Ph$X$", {4, 2}, singleq}},
      {OpType::NPhasedX, {"NPhasedX", "n-Ph$X$", {4, 2}, std::nullopt}},
      {OpType::ZZMax, {"ZZMax", "$ZZ(\\frac{\\pi}{4})$", {}, doubleq}},
      {OpType::XXPhase, {"XXPhase", "$R_{XX}$", {4}, doubleq}},
      {OpType::YYPhase, {"YYPhase", "$R_{YY}$", {4}, doubleq}},
      {OpType::ZZPhase, {"ZZPhase", "$R_{ZZ}$", {4}, doubleq}},
      {OpType::XXPhase3,
       {"XXPhase3", "$R_{X_0X_1}R_{X_0X_2}R_{X_1X_2}$", {4}, tripleq}},
      {OpType::CnRy, {"CnRy", "CnRy", {4}, std::nullopt}},
      {OpType::CnX, {"CnX", "CnX", {}, std::nullopt}},
      {OpType::TK1, {"TK1", "TK1", {4, 4, 4}, singleq}},
      {OpType::TK2, {"TK2", "TK2", {4, 4, 4}, doubleq}},
      {OpType::ESWAP, {"ESWAP", "$\\mathrm{eSWAP}$", {4}, doubleq}},
      {OpType::FSim, {"FSim", "$\\mathrm{fSim}$", {2, 2}, doubleq}},
      {OpType::Sycamore, {"Sycamore", "\\mathrm{Syc}", {}, doubleq}},
      {OpType::ISWAPMax, {"ISWAPMax", "ISWAP", {}, doubleq}},
      {OpType::PhasedISWAP, {"PhasedISWAP", "PhasedISWAP", {1, 4}, doubleq}},
      {OpType::CliffBox, {"CliffBox", "Clifford", {}, std::nullopt}},
      {OpType::PhasePolyBox,
       {"PhasePolyBox", "PhasePolyBox", {}, std::nullopt}},
      {OpType::QControlBox, {"QControlBox", "Ctrl", {}, std::nullopt}},
      {OpType::Conditional, {"Conditional", "If", {}, std::nullopt}},
      {OpType::ProjectorAssertionBox,
       {"ProjectorAssertionBox", "ProjectorAssertionBox", {}, std::nullopt}},
      {OpType::StabiliserAssertionBox,
       {"StabiliserAssertionBox", "StabiliserAssertionBox", {}, std::nullopt}},
      {OpType::ClassicalTransform,
       {"ClassicalTransform", "ClassicalTransform", {}, std::nullopt}},
      {OpType::WASM, {"WASM", "WASM", {}, std::nullopt}},
      {OpType::SetBits, {"SetBits", "SetBits", {}, std::nullopt}},
      {OpType::CopyBits, {"CopyBits", "CopyBits", {}, std::nullopt}},
      {OpType::RangePredicate,
       {"RangePredicate", "RangePredicate", {}, std::nullopt}},
      {OpType::ExplicitPredicate,
       {"ExplicitPredicate", "ExplicitPredicate", {}, std::nullopt}},
      {OpType::ExplicitModifier,
       {"ExplicitModifier", "ExplicitModifier", {}, std::nullopt}},
      {OpType::ClassicalExpBox,
       {"ClassicalExpBox", "ClassicalExpBox", {}, std::nullopt}},
      {OpType::MultiBit, {"MultiBit", "MultiBit", {}, std::nullopt}},
      {OpType::UnitaryTableauBox,
       {"UnitaryTableauBox", "UnitaryTableauBox", {}, std::nullopt}}};
  static std::unique_ptr<const std::map<OpType, OpTypeInfo>> opinfo =
      std::make_unique<const std::map<OpType, OpTypeInfo>>(typeinfo);
  return *opinfo;
}

}  // namespace tket
