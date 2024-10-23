// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tket/Ops/ClExpr.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>
#include <string>
#include <tkassert/Assert.hpp>
#include <variant>
#include <vector>

#include "tket/OpType/OpType.hpp"

namespace tket {

std::ostream& operator<<(std::ostream& os, ClOp fn) {
  switch (fn) {
    case ClOp::INVALID:
      return os << "INVALID";
    case ClOp::BitAnd:
      return os << "and";
    case ClOp::BitOr:
      return os << "or";
    case ClOp::BitXor:
      return os << "xor";
    case ClOp::BitEq:
      return os << "eq";
    case ClOp::BitNeq:
      return os << "neq";
    case ClOp::BitNot:
      return os << "not";
    case ClOp::BitZero:
      return os << "zero";
    case ClOp::BitOne:
      return os << "one";
    case ClOp::RegAnd:
      return os << "and";
    case ClOp::RegOr:
      return os << "or";
    case ClOp::RegXor:
      return os << "xor";
    case ClOp::RegEq:
      return os << "eq";
    case ClOp::RegNeq:
      return os << "neq";
    case ClOp::RegNot:
      return os << "not";
    case ClOp::RegZero:
      return os << "zero";
    case ClOp::RegOne:
      return os << "one";
    case ClOp::RegLt:
      return os << "lt";
    case ClOp::RegGt:
      return os << "gt";
    case ClOp::RegLeq:
      return os << "leq";
    case ClOp::RegGeq:
      return os << "geq";
    case ClOp::RegAdd:
      return os << "add";
    case ClOp::RegSub:
      return os << "sub";
    case ClOp::RegMul:
      return os << "mul";
    case ClOp::RegDiv:
      return os << "div";
    case ClOp::RegPow:
      return os << "pow";
    case ClOp::RegLsh:
      return os << "lsh";
    case ClOp::RegRsh:
      return os << "rsh";
    case ClOp::RegNeg:
      return os << "neg";
  }
  throw std::logic_error("Invalid data");
}

std::ostream& operator<<(std::ostream& os, const ClBitVar& var) {
  return os << "b" << var.index;
}

std::ostream& operator<<(std::ostream& os, const ClRegVar& var) {
  return os << "r" << var.index;
}

std::ostream& operator<<(std::ostream& os, const ClExprVar& var) {
  if (const ClBitVar* bvar = std::get_if<ClBitVar>(&var)) {
    return os << *bvar;
  } else {
    ClRegVar rvar = std::get<ClRegVar>(var);
    return os << rvar;
  }
}

void to_json(nlohmann::json& j, const ClExprVar& var) {
  nlohmann::json inner_j;
  if (const ClBitVar* bvar = std::get_if<ClBitVar>(&var)) {
    j["type"] = "bit";
    to_json(inner_j, *bvar);
  } else {
    j["type"] = "reg";
    ClRegVar rvar = std::get<ClRegVar>(var);
    to_json(inner_j, rvar);
  }
  j["var"] = inner_j;
}

void from_json(const nlohmann::json& j, ClExprVar& var) {
  const std::string vartype = j.at("type").get<std::string>();
  if (vartype == "bit") {
    var = j.at("var").get<ClBitVar>();
  } else {
    TKET_ASSERT(vartype == "reg");
    var = j.at("var").get<ClRegVar>();
  }
}

std::ostream& operator<<(std::ostream& os, const ClExprTerm& term) {
  if (const int* n = std::get_if<int>(&term)) {
    return os << *n;
  } else {
    ClExprVar var = std::get<ClExprVar>(term);
    return os << var;
  }
}

void to_json(nlohmann::json& j, const ClExprTerm& term) {
  nlohmann::json inner_j;
  if (const int* n = std::get_if<int>(&term)) {
    j["type"] = "int";
    inner_j = *n;
  } else {
    j["type"] = "var";
    ClExprVar var = std::get<ClExprVar>(term);
    to_json(inner_j, var);
  }
  j["term"] = inner_j;
}

void from_json(const nlohmann::json& j, ClExprTerm& term) {
  const std::string termtype = j.at("type").get<std::string>();
  if (termtype == "int") {
    term = j.at("term").get<int>();
  } else {
    TKET_ASSERT(termtype == "var");
    term = j.at("term").get<ClExprVar>();
  }
}

std::ostream& operator<<(std::ostream& os, const ClExprArg& arg) {
  if (const ClExprTerm* term = std::get_if<ClExprTerm>(&arg)) {
    return os << *term;
  } else {
    ClExpr expr = std::get<ClExpr>(arg);
    return os << expr;
  }
}

void to_json(nlohmann::json& j, const ClExprArg& arg) {
  nlohmann::json inner_j;
  if (const ClExprTerm* term = std::get_if<ClExprTerm>(&arg)) {
    j["type"] = "term";
    to_json(inner_j, *term);
  } else {
    j["type"] = "expr";
    ClExpr expr = std::get<ClExpr>(arg);
    to_json(inner_j, expr);
  }
  j["input"] = inner_j;
}

void from_json(const nlohmann::json& j, ClExprArg& arg) {
  const std::string inputtype = j.at("type").get<std::string>();
  if (inputtype == "term") {
    arg = j.at("input").get<ClExprTerm>();
  } else {
    TKET_ASSERT(inputtype == "expr");
    ClExpr expr;
    from_json(j.at("input"), expr);
    arg = expr;
  }
}

ClExpr::ClExpr() : ClExpr(ClOp::INVALID, {}) {}

ClExpr::ClExpr(ClOp op, std::vector<ClExprArg> args)
    : op(op), args(args), all_bit_vars(), all_reg_vars() {
  for (const ClExprArg& input : args) {
    if (std::holds_alternative<ClExprTerm>(input)) {
      ClExprTerm basic_input = std::get<ClExprTerm>(input);
      if (std::holds_alternative<ClExprVar>(basic_input)) {
        ClExprVar var = std::get<ClExprVar>(basic_input);
        if (std::holds_alternative<ClBitVar>(var)) {
          ClBitVar bit_var = std::get<ClBitVar>(var);
          all_bit_vars.insert(bit_var.index);
        } else {
          ClRegVar reg_var = std::get<ClRegVar>(var);
          all_reg_vars.insert(reg_var.index);
        }
      }
    } else {
      ClExpr expr = std::get<ClExpr>(input);
      std::set<unsigned> expr_bit_vars = expr.all_bit_variables();
      std::set<unsigned> expr_reg_vars = expr.all_reg_variables();
      all_bit_vars.insert(expr_bit_vars.begin(), expr_bit_vars.end());
      all_reg_vars.insert(expr_reg_vars.begin(), expr_reg_vars.end());
    }
  }
}

bool ClExpr::operator==(const ClExpr& other) const {
  return op == other.op && args == other.args;
}

std::ostream& operator<<(std::ostream& os, const ClExpr& expr) {
  os << expr.get_op() << "(";
  const std::vector<ClExprArg>& args = expr.get_args();
  unsigned n_args = args.size();
  for (unsigned i = 0; i < n_args; i++) {
    os << args[i];
    if (i + 1 < n_args) {
      os << ", ";
    }
  }
  os << ")";
  return os;
}

ClOp ClExpr::get_op() const { return op; }

std::vector<ClExprArg> ClExpr::get_args() const { return args; }

std::set<unsigned> ClExpr::all_bit_variables() const { return all_bit_vars; }

std::set<unsigned> ClExpr::all_reg_variables() const { return all_reg_vars; }

void to_json(nlohmann::json& j, const ClExpr& expr) {
  nlohmann::json j_op = expr.get_op();
  nlohmann::json j_args = expr.get_args();
  j["op"] = j_op;
  j["args"] = j_args;
}

void from_json(const nlohmann::json& j, ClExpr& expr) {
  ClOp op = j.at("op").get<ClOp>();
  std::vector<ClExprArg> args = j.at("args").get<std::vector<ClExprArg>>();
  expr = ClExpr(op, args);
}

WiredClExpr::WiredClExpr() : WiredClExpr({}, {}, {}, {}) {}

WiredClExpr::WiredClExpr(
    const ClExpr& expr, const std::map<unsigned, unsigned>& bit_posn,
    const std::map<unsigned, std::vector<unsigned>>& reg_posn,
    const std::vector<unsigned> output_posn)
    : expr(expr),
      bit_posn(bit_posn),
      reg_posn(reg_posn),
      output_posn(output_posn) {
  std::set<unsigned> b;
  std::set<unsigned> r;
  std::set<unsigned> posns;
  for (const auto& pair : bit_posn) {
    b.insert(pair.first);
    unsigned bit_pos = pair.second;
    if (posns.contains(bit_pos)) {
      throw ClExprWiringError("Invalid maps constructing WiredClExpr");
    }
    posns.insert(bit_pos);
    all_bit_posns.insert(bit_pos);
  }
  for (const auto& pair : reg_posn) {
    r.insert(pair.first);
    for (unsigned bit_pos : pair.second) {
      if (posns.contains(bit_pos)) {
        throw ClExprWiringError("Invalid maps constructing WiredClExpr");
      }
      posns.insert(bit_pos);
    }
    all_reg_posns.insert(pair.second);
  }
  total_n_bits = posns.size();
  for (const unsigned& posn : output_posn) {
    if (!posns.contains(posn)) {
      total_n_bits++;
    }
  }
  if (output_posn.size() == 1) {
    // It mustn't be one of an input register of size > 1
    unsigned i = output_posn[0];
    for (const std::vector<unsigned>& reg : all_reg_posns) {
      if (reg.size() > 1 && std::any_of(
                                reg.begin(), reg.end(),
                                [&i](const unsigned& j) { return i == j; })) {
        throw ClExprWiringError(
            "Output bit contained in a larger input register");
      }
    }
  } else {
    // It must either be disjoint from everything or match one of the registers
    if (std::any_of(
            output_posn.begin(), output_posn.end(),
            [&posns](const unsigned& j) { return posns.contains(j); })) {
      if (!std::any_of(
              all_reg_posns.begin(), all_reg_posns.end(),
              [&output_posn](const std::vector<unsigned>& reg) {
                return output_posn == reg;
              })) {
        throw ClExprWiringError("Output register inconsistent with inputs");
      }
    }
  }
  if (b != expr.all_bit_variables()) {
    throw ClExprWiringError(
        "Mismatch of bit variables constructing WiredClExpr");
  }
  if (r != expr.all_reg_variables()) {
    throw ClExprWiringError(
        "Mismatch of register variables constructing WiredClExpr");
  }
}

bool WiredClExpr::operator==(const WiredClExpr& other) const {
  return expr == other.expr && bit_posn == other.bit_posn &&
         reg_posn == other.reg_posn && output_posn == other.output_posn;
}

std::ostream& operator<<(std::ostream& os, const WiredClExpr& expr) {
  os << expr.expr << " [";
  unsigned n_vars = expr.bit_posn.size() + expr.reg_posn.size();
  unsigned i = 0;
  for (const std::pair<unsigned, unsigned> pair : expr.bit_posn) {
    os << "b" << pair.first << ":" << pair.second;
    i++;
    if (i < n_vars) {
      os << ", ";
    }
  }
  for (const std::pair<unsigned, std::vector<unsigned>> pair : expr.reg_posn) {
    os << "r" << pair.first << ":(";
    unsigned reg_size = pair.second.size();
    for (unsigned j = 0; j < reg_size; j++) {
      os << pair.second[j];
      if (j + 1 < reg_size) {
        os << ",";
      }
    }
    os << ")";
    i++;
    if (i < n_vars) {
      os << ", ";
    }
  }
  os << " --> (";
  unsigned n_outs = expr.output_posn.size();
  for (unsigned i = 0; i < n_outs; i++) {
    os << expr.output_posn[i];
    if (i + 1 < n_outs) {
      os << ",";
    }
  }
  os << ")]";
  return os;
}

ClExpr WiredClExpr::get_expr() const { return expr; }

std::map<unsigned, unsigned> WiredClExpr::get_bit_posn() const {
  return bit_posn;
}

std::map<unsigned, std::vector<unsigned>> WiredClExpr::get_reg_posn() const {
  return reg_posn;
}

std::vector<unsigned> WiredClExpr::get_output_posn() const {
  return output_posn;
}

unsigned WiredClExpr::get_total_n_bits() const { return total_n_bits; }

void to_json(nlohmann::json& j, const WiredClExpr& expr) {
  nlohmann::json j_expr = expr.get_expr();
  nlohmann::json j_bit_posn = expr.get_bit_posn();
  nlohmann::json j_reg_posn = expr.get_reg_posn();
  nlohmann::json j_output_posn = expr.get_output_posn();
  j["expr"] = j_expr;
  j["bit_posn"] = j_bit_posn;
  j["reg_posn"] = j_reg_posn;
  j["output_posn"] = j_output_posn;
}

void from_json(const nlohmann::json& j, WiredClExpr& expr) {
  ClExpr e = j.at("expr").get<ClExpr>();
  std::map<unsigned, unsigned> bit_posn =
      j.at("bit_posn").get<std::map<unsigned, unsigned>>();
  std::map<unsigned, std::vector<unsigned>> reg_posn =
      j.at("reg_posn").get<std::map<unsigned, std::vector<unsigned>>>();
  std::vector<unsigned> output_posn =
      j.at("output_posn").get<std::vector<unsigned>>();
  expr = WiredClExpr(e, bit_posn, reg_posn, output_posn);
}

ClExprOp::ClExprOp(const WiredClExpr& expr) : Op(OpType::ClExpr), expr(expr) {}

Op_ptr ClExprOp::symbol_substitution(const SymEngine::map_basic_basic&) const {
  return std::make_shared<ClExprOp>(*this);
}

SymSet ClExprOp::free_symbols() const { return SymSet(); }

op_signature_t ClExprOp::get_signature() const {
  return op_signature_t(expr.get_total_n_bits(), EdgeType::Classical);
}

WiredClExpr ClExprOp::get_wired_expr() const { return expr; }

nlohmann::json ClExprOp::serialize() const {
  nlohmann::json j;
  j["type"] = get_type();
  j["expr"] = get_wired_expr();
  return j;
}

Op_ptr ClExprOp::deserialize(const nlohmann::json& j) {
  ClExprOp exprop{j.at("expr").get<WiredClExpr>()};
  return std::make_shared<const ClExprOp>(exprop);
}

}  // namespace tket
