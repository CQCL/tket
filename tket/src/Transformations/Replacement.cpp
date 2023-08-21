// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tket/Transformations/Replacement.hpp"

#include "tket/Circuit/CircPool.hpp"
#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Gate/GatePtr.hpp"
#include "tket/Gate/GateUnitaryMatrix.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/OpType/OpTypeInfo.hpp"
#include "tket/Transformations/Decomposition.hpp"
#include "tket/Transformations/Transform.hpp"

namespace tket {

using namespace Transforms;

Circuit multi_controlled_to_2q(
    const Op_ptr op, const std::optional<OpType>& two_q_type) {
  unsigned n_qubits = op->n_qubits();
  OpType optype = op->get_type();
  Circuit c(n_qubits);
  switch (optype) {
    case OpType::CnRy:
      c = CircPool::CnRy_normal_decomp(op, n_qubits);
      break;
    case OpType::CnX:
    case OpType::CnZ:
    case OpType::CnY:
      if (n_qubits >= 6 && n_qubits <= 50) {
        // CnU_linear_depth_decomp performs better in this case
        OpType target_type =
            (optype == OpType::CnX)
                ? OpType::X
                : ((optype == OpType::CnZ) ? OpType::Z : OpType::Y);
        Eigen::Matrix2cd target_u =
            GateUnitaryMatrix::get_unitary(target_type, 1, {});
        c = CircPool::CnU_linear_depth_decomp(n_qubits - 1, target_u);
      } else {
        if (optype == OpType::CnZ) {
          c.add_op<unsigned>(OpType::H, {n_qubits - 1});
        } else if (optype == OpType::CnY) {
          c.add_op<unsigned>(OpType::Sdg, {n_qubits - 1});
        }
        c.append(CircPool::CnX_normal_decomp(n_qubits - 1));
        if (optype == OpType::CnZ) {
          c.add_op<unsigned>(OpType::H, {n_qubits - 1});
        } else if (optype == OpType::CnY) {
          c.add_op<unsigned>(OpType::S, {n_qubits - 1});
        }
      }
      break;
    default:
      throw BadOpType("The operation is not multi-controlled", optype);
  }

  if (two_q_type == std::nullopt) {
    return c;
  } else if (two_q_type.value() == OpType::CX) {
    decompose_multi_qubits_CX().apply(c);
  } else if (two_q_type.value() == OpType::TK2) {
    decompose_multi_qubits_TK2().apply(c);
  } else {
    throw BadOpType("The target 2-q gate can only be CX or TK2", optype);
  }
  return c;
}

Circuit TK2_circ_from_multiq(const Op_ptr op) {
  OpDesc desc = op->get_desc();
  std::vector<Expr> params = op->get_params();
  if (!desc.is_gate())
    throw BadOpType(
        "Can only build replacement circuits for basic gates", desc.type());
  switch (desc.type()) {
    case OpType::CnRy:
    case OpType::CnX:
    case OpType::CnZ:
    case OpType::CnY:
      // TODO We should be able to do better than this.
      return multi_controlled_to_2q(op, OpType::TK2);
    case OpType::XXPhase:
      return CircPool::XXPhase_using_TK2(params[0]);
    case OpType::YYPhase:
      return CircPool::YYPhase_using_TK2(params[0]);
    case OpType::ZZPhase:
      return CircPool::ZZPhase_using_TK2(params[0]);
    default:
      return with_TK2(as_gate_ptr(op));
  }
}

Circuit CX_circ_from_multiq(const Op_ptr op) {
  OpDesc desc = op->get_desc();
  if (!desc.is_gate())
    throw BadOpType(
        "Can only build replacement circuits for basic gates", desc.type());
  switch (desc.type()) {
    case OpType::CnRy:
    case OpType::CnX:
    case OpType::CnZ:
    case OpType::CnY:
      return multi_controlled_to_2q(op, OpType::CX);
    default:
      return with_CX(as_gate_ptr(op));
  }
}

Circuit CX_ZX_circ_from_op(const Op_ptr op) {
  OpDesc desc = op->get_desc();
  if (!desc.is_gate())
    throw BadOpType(
        "Can only build replacement circuits for basic gates", desc.type());
  switch (desc.type()) {
    case OpType::Phase: {
      Circuit replacement(0);
      replacement.add_phase(op->get_params()[0]);
      return replacement;
    }
    case OpType::Z: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rz, 1., {0});
      replacement.add_phase(0.5);
      return replacement;
    }
    case OpType::X: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rx, 1., {0});
      replacement.add_phase(0.5);
      return replacement;
    }
    case OpType::Y: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rz, 1., {0});
      replacement.add_op<unsigned>(OpType::Rx, 1., {0});
      replacement.add_phase(-0.5);
      return replacement;
    }
    case OpType::S: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rz, 0.5, {0});
      replacement.add_phase(0.25);
      return replacement;
    }
    case OpType::Sdg: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rz, -0.5, {0});
      replacement.add_phase(-0.25);
      return replacement;
    }
    case OpType::T: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rz, 0.25, {0});
      replacement.add_phase(0.125);
      return replacement;
    }
    case OpType::Tdg: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rz, -0.25, {0});
      replacement.add_phase(-0.125);
      return replacement;
    }
    case OpType::V: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rx, 0.5, {0});
      return replacement;
    }
    case OpType::Vdg: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rx, -0.5, {0});
      return replacement;
    }
    case OpType::SX: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rx, 0.5, {0});
      replacement.add_phase(0.25);
      return replacement;
    }
    case OpType::SXdg: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rx, -0.5, {0});
      replacement.add_phase(-0.25);
      return replacement;
    }
    case OpType::H: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rz, 0.5, {0});
      replacement.add_op<unsigned>(OpType::Rx, 0.5, {0});
      replacement.add_op<unsigned>(OpType::Rz, 0.5, {0});
      replacement.add_phase(0.5);
      return replacement;
    }
    case OpType::Ry: {
      Expr angle = op->get_params()[0];
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Rz, -0.5, {0});
      replacement.add_op<unsigned>(OpType::Rx, angle, {0});
      replacement.add_op<unsigned>(OpType::Rz, 0.5, {0});
      return replacement;
    }
    case OpType::Rx:
    case OpType::Rz:
    case OpType::Measure:
    case OpType::Collapse: {
      Circuit replacement(1);
      replacement.add_op<unsigned>(op, {0});
      return replacement;
    }
    case OpType::U3: {
      Circuit replacement(1);
      Expr angle_z1 = op->get_params()[2];
      Expr angle_y = op->get_params()[0];
      Expr angle_z2 = op->get_params()[1];
      replacement.add_op<unsigned>(OpType::Rz, angle_z1 - 0.5, {0});
      replacement.add_op<unsigned>(OpType::Rx, angle_y, {0});
      replacement.add_op<unsigned>(OpType::Rz, angle_z2 + 0.5, {0});
      replacement.add_phase((angle_z1 + angle_z2) / 2);
      return replacement;
    }
    case OpType::U2: {
      Circuit replacement(1);
      Expr angle_z1 = op->get_params()[1];
      Expr angle_z2 = op->get_params()[0];
      replacement.add_op<unsigned>(OpType::Rz, angle_z1 - 0.5, {0});
      replacement.add_op<unsigned>(OpType::Rx, 0.5, {0});
      replacement.add_op<unsigned>(OpType::Rz, angle_z2 + 0.5, {0});
      replacement.add_phase((angle_z1 + angle_z2) / 2);
      return replacement;
    }
    case OpType::U1: {
      Circuit replacement(1);
      Expr angle = op->get_params()[0];
      replacement.add_op<unsigned>(OpType::Rz, angle, {0});
      replacement.add_phase(angle / 2);
      return replacement;
    }
    case OpType::PhasedX: {
      Circuit replacement(1);
      Expr theta = op->get_params()[0];
      Expr phi = op->get_params()[1];
      replacement.add_op<unsigned>(OpType::Rz, -phi, {0});
      replacement.add_op<unsigned>(OpType::Rx, theta, {0});
      replacement.add_op<unsigned>(OpType::Rz, phi, {0});
      return replacement;
    }
    case OpType::CX: {
      Circuit replacement(2);
      replacement.add_op<unsigned>(op, {0, 1});
      return replacement;
    }
    case OpType::TK2:
    case OpType::CY:
    case OpType::CZ:
    case OpType::CH:
    case OpType::CV:
    case OpType::CVdg:
    case OpType::CSX:
    case OpType::CSXdg:
    case OpType::CRz:
    case OpType::CRx:
    case OpType::CRy:
    case OpType::CU1:
    case OpType::CU3:
    case OpType::PhaseGadget:
    case OpType::CCX:
    case OpType::SWAP:
    case OpType::CSWAP:
    case OpType::ECR:
    case OpType::ISWAP:
    case OpType::XXPhase:
    case OpType::XXPhase3:
    case OpType::ZZMax:
    case OpType::ZZPhase:
    case OpType::YYPhase:
    case OpType::CnRy:
    case OpType::CnX:
    case OpType::ESWAP:
    case OpType::FSim:
    case OpType::Sycamore:
    case OpType::ISWAPMax:
    case OpType::BRIDGE: {
      Circuit replacement = CX_circ_from_multiq(op);
      decompose_ZX().apply(replacement);
      return replacement;
    }
    case OpType::TK1: {
      Circuit replacement(1);
      std::vector<Expr> params = op->get_params();
      replacement.add_op<unsigned>(OpType::Rz, params[2], {0});
      replacement.add_op<unsigned>(OpType::Rx, params[1], {0});
      replacement.add_op<unsigned>(OpType::Rz, params[0], {0});
      return replacement;
    }
    default:
      throw BadOpType("Cannot find replacement circuit", desc.type());
  }
}

}  // namespace tket
