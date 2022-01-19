// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "Gate.hpp"

#include <vector>

#include "GateUnitaryMatrix.hpp"
#include "GateUnitaryMatrixError.hpp"
#include "OpPtrFunctions.hpp"
#include "OpType/OpType.hpp"
#include "OpType/OpTypeInfo.hpp"

namespace tket {

Op_ptr Gate::dagger() const {
  OpType optype = get_type();
  switch (optype) {
    case OpType::H:
    case OpType::X:
    case OpType::Y:
    case OpType::Z:
    case OpType::SWAP:
    case OpType::CH:
    case OpType::CX:
    case OpType::CY:
    case OpType::CZ:
    case OpType::CCX:
    case OpType::noop:
    case OpType::CSWAP:
    case OpType::ECR:
    case OpType::BRIDGE: {
      return get_op_ptr(optype);
    }
    case OpType::CnX: {
      return get_op_ptr(optype, std::vector<Expr>(), n_qubits_);
    }
    case OpType::S: {
      return get_op_ptr(OpType::Sdg);
    }
    case OpType::Sdg: {
      return get_op_ptr(OpType::S);
    }
    case OpType::T: {
      return get_op_ptr(OpType::Tdg);
    }
    case OpType::Tdg: {
      return get_op_ptr(OpType::T);
    }
    case OpType::V: {
      return get_op_ptr(OpType::Vdg);
    }
    case OpType::Vdg: {
      return get_op_ptr(OpType::V);
    }
    case OpType::CV: {
      return get_op_ptr(OpType::CVdg);
    }
    case OpType::CVdg: {
      return get_op_ptr(OpType::CV);
    }
    case OpType::SX: {
      return get_op_ptr(OpType::SXdg);
    }
    case OpType::SXdg: {
      return get_op_ptr(OpType::SX);
    }
    case OpType::CSX: {
      return get_op_ptr(OpType::CSXdg);
    }
    case OpType::CSXdg: {
      return get_op_ptr(OpType::CSX);
    }
    case OpType::CRz:
    case OpType::CRx:
    case OpType::CRy:
    case OpType::CU1:
    case OpType::U1:
    case OpType::Rz:
    case OpType::Ry:
    case OpType::Rx:
    case OpType::PhaseGadget:
    case OpType::CnRy:
    case OpType::XXPhase:
    case OpType::YYPhase:
    case OpType::ZZPhase:
    case OpType::XXPhase3:
    case OpType::ISWAP:
    case OpType::ESWAP: {
      return get_op_ptr(optype, -params_[0], n_qubits_);
    }
    case OpType::ZZMax: {
      // ZZMax.dagger = ZZPhase(-0.5)
      return get_op_ptr(OpType::ZZPhase, -0.5);
    }
    case OpType::FSim: {
      // FSim(a,b).dagger() == FSim(-a,-b)
      return get_op_ptr(OpType::FSim, {-params_[0], -params_[1]});
    }
    case OpType::Sycamore: {
      return get_op_ptr(OpType::FSim, {-0.5, -1. / 6.});
    }
    case OpType::ISWAPMax: {
      return get_op_ptr(OpType::ISWAP, 3.);
    }
    case OpType::U2: {
      // U2(a,b).dagger() == U3(-pi/2,-b,-a)
      return get_op_ptr(OpType::U3, {-0.5, -params_[1], -params_[0]});
    }
    case OpType::U3:
    case OpType::CU3:
      // U3(a,b,c).dagger() == U3(-a,-c.-b)
      { return get_op_ptr(optype, {-params_[0], -params_[2], -params_[1]}); }
    case OpType::TK1:
      // TK1(a,b,c).dagger() == TK1(-c,-b,-a)
      {
        return get_op_ptr(OpType::TK1, {-params_[2], -params_[1], -params_[0]});
      }
    case OpType::PhasedX:
    case OpType::NPhasedX:
      // PhasedX(a,b).dagger() == PhasedX(-a,b)
      { return get_op_ptr(optype, {-params_[0], params_[1]}, n_qubits_); }
    case OpType::PhasedISWAP:
      // PhasedISWAP(a,b).dagger() == PhasedISWAP(a,-b)
      { return get_op_ptr(OpType::PhasedISWAP, {params_[0], -params_[1]}); }
    default: {
      throw NotImplemented(
          "Cannot retrieve the dagger of OpType::" + get_desc().name());
    }
  }
}

Op_ptr Gate::transpose() const {
  OpType optype = get_type();
  switch (optype) {
    case OpType::H:
    case OpType::X:
    case OpType::Z:
    case OpType::SWAP:
    case OpType::CH:
    case OpType::CX:
    case OpType::CZ:
    case OpType::CV:
    case OpType::CVdg:
    case OpType::CSX:
    case OpType::CSXdg:
    case OpType::CCX:
    case OpType::noop:
    case OpType::CSWAP:
    case OpType::BRIDGE:
    case OpType::S:
    case OpType::T:
    case OpType::V:
    case OpType::Vdg:
    case OpType::SX:
    case OpType::SXdg:
    case OpType::CRz:
    case OpType::CRx:
    case OpType::CU1:
    case OpType::U1:
    case OpType::Rz:
    case OpType::Rx:
    case OpType::PhaseGadget:
    case OpType::XXPhase:
    case OpType::YYPhase:
    case OpType::ZZPhase:
    case OpType::XXPhase3:
    case OpType::ESWAP:
    case OpType::FSim: {
      return get_op_ptr(optype, params_);
    }
    case OpType::Y: {
      return get_op_ptr(OpType::U3, {3, 0.5, 0.5});
    }
    case OpType::Ry:
    case OpType::CRy:
    case OpType::CnRy: {
      return get_op_ptr(optype, -params_[0], n_qubits_);
    }
    case OpType::CnX: {
      return get_op_ptr(optype, std::vector<Expr>(), n_qubits_);
    }
    case OpType::U2: {
      // U2(a,b).transpose() == U2(b+1,a+1)
      return get_op_ptr(OpType::U2, {params_[1] + 1., params_[0] + 1.});
    }
    case OpType::U3:
    case OpType::CU3: {
      // U3(a,b,c).transpose() == U3(-a,c,b)
      return get_op_ptr(OpType::U3, {-params_[0], params_[2], params_[1]});
    }
    case OpType::TK1: {
      // TK1(a,b,c).transpose() == TK1(c,b,a)
      return get_op_ptr(OpType::TK1, {params_[2], params_[1], params_[0]});
    }
    case OpType::PhasedX:
    case OpType::NPhasedX: {
      // PhasedX(a,b).transpose() == PhasedX(a,-b)
      return get_op_ptr(optype, {params_[0], -params_[1]}, n_qubits_);
    }
    case OpType::PhasedISWAP: {
      // PhasedISWAP(a,b).transpose() == PhasedISWAP(-a,b)
      return get_op_ptr(OpType::PhasedISWAP, {-params_[0], params_[1]});
    }

    default: {
      throw NotImplemented(
          "Cannot retrieve the transpose of OpType::" + get_desc().name());
    }
  }
}

Op_ptr Gate::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  std::vector<Expr> new_params;
  for (const Expr& p : this->params_) {
    new_params.push_back(p.subs(sub_map));
  }
  return get_op_ptr(this->type_, new_params, this->n_qubits_);
}

std::optional<double> Gate::is_identity() const {
  static const std::optional<double> notid;
  const std::vector<Expr>& params = get_params();
  switch (get_type()) {
    case OpType::Rx:
    case OpType::Ry:
    case OpType::Rz:
    case OpType::PhasedX:
    case OpType::NPhasedX:
    case OpType::XXPhase:
    case OpType::YYPhase:
    case OpType::ZZPhase:
    case OpType::XXPhase3:
    case OpType::ESWAP: {
      Expr e = params[0];
      if (equiv_0(e, 4)) {
        return 0.;
      } else if (equiv_0(e + 2, 4)) {
        return 1.;
      } else
        return notid;
    }
    case OpType::U1:
    case OpType::CU1: {
      return equiv_0(params[0]) ? 0. : notid;
    }
    case OpType::U3: {
      Expr theta = params[0];
      if (equiv_0(params[1] + params[2])) {
        if (equiv_0(theta, 4)) {
          return 0.;
        } else if (equiv_0(theta + 2, 4)) {
          return 1.;
        } else
          return notid;
      } else
        return notid;
    }
    case OpType::CU3: {
      if (equiv_0(params[0], 4) && equiv_0(params[1] + params[2])) {
        return 0.;
      } else
        return notid;
    }
    case OpType::TK1: {
      Expr s = params[0] + params[2], t = params[1];
      if (equiv_0(s) && equiv_0(t)) {
        return (equiv_0(s, 4) ^ equiv_0(t, 4)) ? 1. : 0.;
      } else
        return notid;
    }
    case OpType::CRz:
    case OpType::CRx:
    case OpType::CRy:
    case OpType::PhaseGadget:
    case OpType::ISWAP:
    case OpType::CnRy: {
      return equiv_0(params[0], 4) ? 0. : notid;
    }
    case OpType::FSim: {
      return (equiv_0(params[0]) && equiv_0(params[1])) ? 0. : notid;
    }
    case OpType::PhasedISWAP: {
      return (equiv_0(params[0], 1) && equiv_0(params[1], 4)) ? 0. : notid;
    }
    default:
      return notid;
  }
}

std::string Gate::get_name(bool latex) const {
  OpDesc desc = get_desc();
  if (params_.size() > 0) {
    std::stringstream name;
    if (latex) {
      name << desc.latex() << "(";
    } else {
      name << desc.name() << "(";
    }
    for (unsigned i = 0; i < params_.size(); ++i) {
      std::optional<double> reduced =
          eval_expr_mod(params_[i], desc.param_mod(i));
      if (reduced) {
        name << reduced.value();
      } else {
        name << params_[i];
      }
      if (i < params_.size() - 1) name << ", ";
    }
    name << ")";
    return name.str();
  } else {
    return Op::get_name(latex);
  }
}

std::string Gate::get_command_str(const unit_vector_t& args) const {
  if (get_type() == OpType::Measure) {
    std::stringstream out;
    out << get_name() << " " << args[0].repr() << " --> " << args[1].repr()
        << ";";
    return out.str();
  } else
    return Op::get_command_str(args);
}

unsigned Gate::n_qubits() const {
  OptUInt n = desc_.n_qubits();
  if (n == any) {
    return n_qubits_;
  } else {
    return n.value();
  }
}

bool Gate::is_equal(const Op& op_other) const {
  const Gate& other = dynamic_cast<const Gate&>(op_other);

  OpDesc desc = get_desc();
  if (n_qubits() != other.n_qubits()) return false;
  std::vector<Expr> params1 = this->get_params();
  std::vector<Expr> params2 = other.get_params();
  unsigned param_count = params1.size();
  if (params2.size() != param_count) return false;
  for (unsigned i = 0; i < param_count; i++) {
    if (!equiv_expr(params1[i], params2[i], desc.param_mod(i))) return false;
  }
  return true;
}

std::vector<Expr> Gate::get_params_reduced() const {
  OpDesc desc = get_desc();
  unsigned n_params = desc.n_params();
  std::vector<Expr> params(n_params);
  for (unsigned i = 0; i < n_params; i++) {
    Expr param = params_[i];
    std::optional<double> e = eval_expr_mod(param, desc.param_mod(i));
    if (e) {
      params[i] = e.value();
    } else {
      params[i] = param;
    }
  }
  return params;
}

std::vector<Expr> Gate::get_tk1_angles() const {
  switch (get_type()) {
    case OpType::noop: {
      return {0., 0., 0., 0.};
    }
    case OpType::Z: {
      return {0., 0., 1., 0.5};
    }
    case OpType::X: {
      return {0., 1., 0., 0.5};
    }
    case OpType::Y: {
      return {0.5, 1., -0.5, 0.5};
    }
    case OpType::S: {
      return {0., 0., 0.5, 0.25};
    }
    case OpType::Sdg: {
      return {0., 0., -0.5, -0.25};
    }
    case OpType::T: {
      return {0., 0., 0.25, 0.125};
    }
    case OpType::Tdg: {
      return {0., 0., -0.25, -0.125};
    }
    case OpType::V: {
      return {0., 0.5, 0., 0.};
    }
    case OpType::Vdg: {
      return {0., -0.5, -0., 0.};
    }
    case OpType::SX: {
      return {0., 0.5, 0., 0.25};
    }
    case OpType::SXdg: {
      return {0., -0.5, 0., -0.25};
    }
    case OpType::H: {
      return {0.5, 0.5, 0.5, 0.5};
    }
    case OpType::Rx: {
      return {0., params_.at(0), 0., 0.};
    }
    case OpType::Ry: {
      return {0.5, params_.at(0), -0.5, 0.};
    }
    case OpType::Rz:
    case OpType::PhaseGadget: {
      return {0., 0., params_.at(0), 0.};
    }
    case OpType::U1: {
      return {0., 0., params_.at(0), params_.at(0) / 2};
    }
    case OpType::U2: {
      return {
          params_.at(0) + 0.5, 0.5, params_.at(1) - 0.5,
          (params_.at(0) + params_.at(1)) / 2};
    }
    case OpType::U3: {
      return {
          params_.at(1) + 0.5, params_.at(0), params_.at(2) - 0.5,
          (params_.at(1) + params_.at(2)) / 2};
    }
    case OpType::NPhasedX: {
      if (n_qubits_ != 1) {
        throw NotValid(
            "OpType::NPhasedX can only be decomposed into a TK1 "
            "if it acts on a single qubit");
      }
      return {params_.at(1), params_.at(0), -params_.at(1), 0.};
    }
    case OpType::PhasedX: {
      return {params_.at(1), params_.at(0), -params_.at(1), 0.};
    }
    case OpType::TK1: {
      return {params_.at(0), params_.at(1), params_.at(2), 0.};
    }
    default: {
      throw NotImplemented(
          "Cannot retrieve the TK1 angles of OpType::" + get_desc().name());
    }
  }
}

std::vector<Expr> Gate::get_params() const { return params_; }

SymSet Gate::free_symbols() const { return expr_free_symbols(get_params()); }

std::optional<Pauli> Gate::commuting_basis(port_t port) const {
  unsigned n_q = n_qubits();
  if (port >= n_q) throw NotValid();
  switch (get_type()) {
    case OpType::XXPhase:
    case OpType::XXPhase3: {
      return Pauli::X;
    }
    case OpType::YYPhase: {
      return Pauli::Y;
    }
    case OpType::Z:
    case OpType::S:
    case OpType::Sdg:
    case OpType::T:
    case OpType::Tdg:
    case OpType::Rz:
    case OpType::U1:
    case OpType::CZ:
    case OpType::CRz:
    case OpType::CU1:
    case OpType::PhaseGadget:
    case OpType::ZZMax:
    case OpType::ZZPhase: {
      return Pauli::Z;
    }
    case OpType::noop: {
      return Pauli::I;
    }
    case OpType::H:
    case OpType::U3:
    case OpType::U2:
    case OpType::PhasedX:
    case OpType::NPhasedX:
    case OpType::TK1: {
      return std::nullopt;
    }
    case OpType::CH:
    case OpType::CU3:
    case OpType::CSWAP: {
      if (port == 0) {
        return Pauli::Z;
      } else {
        return std::nullopt;
      }
    }
    case OpType::BRIDGE: {
      if (port == 0) {
        return Pauli::Z;
      } else if (port == 2) {
        return Pauli::X;
      } else {
        return Pauli::I;
      }
    }
    case OpType::X:
    case OpType::V:
    case OpType::Vdg:
    case OpType::CV:
    case OpType::CVdg:
    case OpType::SX:
    case OpType::SXdg:
    case OpType::CSX:
    case OpType::CSXdg:
    case OpType::Rx:
    case OpType::CRx:
    case OpType::CX:
    case OpType::CCX:
    case OpType::CnX: {
      if (port == n_q - 1) {
        return Pauli::X;
      } else {
        return Pauli::Z;
      }
    }
    case OpType::ECR: {
      if (port == 1) {
        return Pauli::X;
      } else {
        return std::nullopt;
      }
    }
    case OpType::Y:
    case OpType::Ry:
    case OpType::CY:
    case OpType::CRy:
    case OpType::CnRy: {
      if (port == n_q - 1) {
        return Pauli::Y;
      } else {
        return Pauli::Z;
      }
    }
    default: {
      return std::nullopt;
    }
  }
}

bool Gate::commutes_with_basis(
    const std::optional<Pauli>& colour, port_t port) const {
  if (colour == Pauli::I) return true;
  const std::optional<Pauli> my_colour = commuting_basis(port);
  if (!colour && !my_colour) return false;
  if (my_colour == Pauli::I || my_colour == colour) return true;
  return false;
}

op_signature_t Gate::get_signature() const {
  std::optional<op_signature_t> sig = desc_.signature();
  if (sig)
    return *sig;
  else
    return op_signature_t(n_qubits_, EdgeType::Quantum);
}

nlohmann::json Gate::serialize() const {
  nlohmann::json j;
  OpType optype = get_type();
  j["type"] = optype;
  // if type has a fixed signature, don't store number of qubits
  if (!optypeinfo().at(optype).signature) {
    j["n_qb"] = n_qubits();
  }
  std::vector<Expr> params = get_params();
  if (!params.empty()) {
    j["params"] = params;
  }
  return j;
}

Op_ptr Gate::deserialize(const nlohmann::json& j) {
  OpType optype = j.at("type").get<OpType>();
  std::vector<Expr> params;
  if (j.contains("params")) {
    params = j.at("params").get<std::vector<Expr>>();
  }
  // if type has fixed number of qubits use it, otherwise it should have been
  // stored
  const auto& sig = optypeinfo().at(optype).signature;
  unsigned n_qb;
  if (sig) {
    auto check_quantum = [](unsigned sum, const EdgeType& e) {
      return sum + (e == EdgeType::Quantum ? 1 : 0);
    };
    n_qb = std::accumulate(sig->begin(), sig->end(), 0, check_quantum);
  } else {
    n_qb = j.at("n_qb").get<unsigned>();
  }
  return get_op_ptr(optype, params, n_qb);
}

Eigen::MatrixXcd Gate::get_unitary() const {
  try {
    return GateUnitaryMatrix::get_unitary(*this);
  } catch (const GateUnitaryMatrixError& e) {
    switch (e.cause) {
      case GateUnitaryMatrixError::Cause::GATE_NOT_IMPLEMENTED:
        throw NotImplemented(e.what());
      // Note that symbolic parameters
      // are "not valid" rather than "not implemented"
      // because the returned matrix is clearly numerical, not symbolic.
      default:
        throw NotValid(e.what());
    }
  }
}

Gate::Gate(OpType type, const std::vector<Expr>& params, unsigned n_qubits)
    : Op(type), params_(params), n_qubits_(n_qubits) {
  if (!is_gate_type(type)) {
    throw NotValid();
  }
  if (params.size() != optypeinfo().at(type).n_params()) {
    throw InvalidParameterCount();
  }
}

Gate::Gate() : Op(OpType::noop), params_() {}

}  // namespace tket
