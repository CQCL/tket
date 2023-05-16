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

#include "tket/Circuit/DiagonalBox.hpp"

#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/Multiplexor.hpp"
#include "tket/Gate/Rotation.hpp"
#include "tket/Ops/OpJsonFactory.hpp"
#include "tket/Utils/HelperFunctions.hpp"
#include "tket/Utils/Json.hpp"

namespace tket {

DiagonalBox::DiagonalBox(const Eigen::VectorXcd &diagonal, bool upper_triangle)
    : Box(OpType::DiagonalBox),
      diagonal_(diagonal),
      upper_triangle_(upper_triangle) {
  size_t length = diagonal.size();
  if (length < 2 || (length & (length - 1)) != 0) {
    throw std::invalid_argument(
        "The size of the diagonal operator passed to DiagonalBox is not a "
        "power of 2.");
  }
  // check the operator is unitary
  for (unsigned i = 0; i < diagonal.size(); i++) {
    if (std::abs(1 - std::abs(diagonal[i])) > EPS) {
      throw std::invalid_argument(
          "The input diagonal passed to DiagonalBox is not unitary.");
    }
  }
}

DiagonalBox::DiagonalBox(const DiagonalBox &other)
    : Box(other),
      diagonal_(other.diagonal_),
      upper_triangle_(other.upper_triangle_) {}

Op_ptr DiagonalBox::dagger() const {
  return std::make_shared<DiagonalBox>(diagonal_.conjugate(), upper_triangle_);
}

Op_ptr DiagonalBox::transpose() const {
  return std::make_shared<DiagonalBox>(diagonal_, upper_triangle_);
}

op_signature_t DiagonalBox::get_signature() const {
  op_signature_t qubits((unsigned)log2(diagonal_.size()), EdgeType::Quantum);
  return qubits;
}

Eigen::VectorXcd DiagonalBox::get_diagonal() const { return diagonal_; }
bool DiagonalBox::is_upper_triangle() const { return upper_triangle_; }

/**
 * @brief Construct a circuit that implements a diagonal operator
 * https://arxiv.org/abs/quant-ph/0406176 Theorem 7
 *
 * @param diagonal
 * @return Circuit
 */
static Circuit diagonal_circ(
    const Eigen::VectorXcd &diagonal, bool upper_triangle) {
  unsigned n_qubits = (unsigned)log2(diagonal.size());
  Circuit circ(n_qubits);
  Eigen::VectorXcd d(diagonal);
  unsigned n_ctrl_qubits = n_qubits - 1;
  while (d.size() != 1) {
    ctrl_op_map_t multip_rz;
    Eigen::VectorXcd new_d(d.size() / 2);
    for (unsigned i = 0; i < new_d.size(); i++) {
      Complex a;
      Complex b;
      if (upper_triangle) {
        a = d[2 * i];
        b = d[2 * i + 1];
      } else {
        a = d[i];
        b = d[i + (unsigned)d.size() / 2];
      }
      // convert diag[a,b] into a p*Rz(alpha)
      double a_phase = std::arg(a);
      double b_phase = std::arg(b);
      double alpha = (b_phase - a_phase) / PI;
      Complex p = std::exp((b_phase + a_phase) * 0.5 * i_);
      std::vector<bool> bitstr = dec_to_bin(i, n_ctrl_qubits);
      if (std::abs(alpha) > EPS) {
        multip_rz.insert({bitstr, get_op_ptr(OpType::Rz, alpha)});
      }
      new_d[i] = p;
    }
    if (!multip_rz.empty()) {
      std::vector<unsigned> args(n_ctrl_qubits);
      if (upper_triangle) {
        std::iota(std::begin(args), std::end(args), 0);
        args.push_back(n_ctrl_qubits);
      } else {
        std::iota(std::begin(args), std::end(args), n_qubits - n_ctrl_qubits);
        args.push_back(n_qubits - n_ctrl_qubits - 1);
      }
      circ.add_box(MultiplexedRotationBox(multip_rz), args);
    }
    d = new_d;
    n_ctrl_qubits--;
  }
  circ.add_phase(std::arg(d[0]) / PI);
  return circ;
};

void DiagonalBox::generate_circuit() const {
  circ_ = std::make_shared<Circuit>(diagonal_circ(diagonal_, upper_triangle_));
}

nlohmann::json DiagonalBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const DiagonalBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["diagonal"] = box.get_diagonal();
  j["upper_triangle"] = box.is_upper_triangle();
  return j;
}

Op_ptr DiagonalBox::from_json(const nlohmann::json &j) {
  DiagonalBox box = DiagonalBox(
      j.at("diagonal").get<Eigen::VectorXcd>(),
      j.at("upper_triangle").get<bool>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(DiagonalBox, DiagonalBox)

}  // namespace tket