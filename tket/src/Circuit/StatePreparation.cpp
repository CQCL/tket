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

#include "tket/Circuit/StatePreparation.hpp"

#include <boost/dynamic_bitset.hpp>

#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/Multiplexor.hpp"
#include "tket/Gate/Rotation.hpp"
#include "tket/Ops/OpJsonFactory.hpp"
#include "tket/Utils/HelperFunctions.hpp"
#include "tket/Utils/Json.hpp"

namespace tket {

StatePreparationBox::StatePreparationBox(
    const Eigen::VectorXcd &statevector, bool is_inverse)
    : Box(OpType::StatePreparationBox),
      statevector_(statevector),
      is_inverse_(is_inverse) {
  size_t length = statevector.size();
  if (length < 2 || (length & (length - 1)) != 0) {
    throw std::invalid_argument(
        "The length of the statevector is not a power of 2.");
  }
  // check the state is normalised
  if (std::abs(statevector.norm() - 1) > EPS) {
    throw std::invalid_argument("The input statevector is not normalised.");
  }
}

StatePreparationBox::StatePreparationBox(const StatePreparationBox &other)
    : Box(other),
      statevector_(other.statevector_),
      is_inverse_(other.is_inverse_) {}

Op_ptr StatePreparationBox::dagger() const {
  return std::make_shared<StatePreparationBox>(statevector_, !is_inverse_);
}

op_signature_t StatePreparationBox::get_signature() const {
  op_signature_t qubits((unsigned)log2(statevector_.size()), EdgeType::Quantum);
  return qubits;
}

Eigen::VectorXcd StatePreparationBox::get_statevector() const {
  return statevector_;
}

bool StatePreparationBox::is_inverse() const { return is_inverse_; }

/**
 * @brief Construct a circuit that prepares an arbitrary quantum state
 * https://arxiv.org/abs/quant-ph/0406176 Theorem 9
 * when is_inverse is set, this function returns the dagger of the state
 * preparation circuit
 *
 * @param statevector
 * @param is_inverse
 * @return Circuit
 */
static Circuit state_prep_circ(
    const Eigen::VectorXcd &statevector, bool is_inverse) {
  unsigned n_qubits = (unsigned)log2(statevector.size());
  Circuit circ(n_qubits);
  std::vector<std::optional<MultiplexedRotationBox>> mutip_ry_vec, mutip_rz_vec;
  Eigen::VectorXcd psi = statevector;
  for (unsigned step = 0; step < n_qubits; step++) {
    unsigned half_length = (unsigned)psi.size() / 2;
    // in each step, we disentangle one qubit from |psi> so that
    // multip_ry*multip_rz*|psi> = rs*exp(ts*i*pi)|0>
    // multip_ry and multip_rz are multiplexed Ry and Rz respectively
    // after this loop, we apply the multiplexors in reverse order to obtain the
    // state prep circuit (if is_inverse is false)
    ctrl_op_map_t multip_ry, multip_rz;
    Eigen::VectorXd rs = Eigen::VectorXd::Zero(half_length);
    Eigen::VectorXd ts = Eigen::VectorXd::Zero(half_length);
    for (unsigned c = 0; c < half_length; c++) {
      Complex a = psi[2 * c];
      Complex b = psi[2 * c + 1];
      if (std::abs(a) < EPS && std::abs(b) < EPS) {
        // if both a and b are zeros
        continue;
      }
      // normalise [a,b]
      double r = 1;
      double norm = a.real() * a.real() + a.imag() * a.imag() +
                    b.real() * b.real() + b.imag() * b.imag();
      if (std::abs(norm - 1) > EPS) {
        r = std::sqrt(norm);
        a = a / r;
        b = b / r;
      }
      double theta, phi, phase;
      std::tie(theta, phi, phase) = get_bloch_coordinate_from_state(a, b);
      ts[c] = phase + 0.5 * phi;
      rs[c] = r;
      unsigned n_controls = n_qubits - step - 1;
      std::vector<bool> bitstr;
      if (n_controls > 0) {
        bitstr = dec_to_bin(c, n_controls);
      }

      double y_angle = is_inverse ? -theta : theta;
      double z_angle = is_inverse ? -phi : phi;
      if (std::abs(y_angle) > EPS) {
        multip_ry.insert({bitstr, get_op_ptr(OpType::Ry, y_angle)});
      }
      if (std::abs(z_angle) > EPS) {
        multip_rz.insert({bitstr, get_op_ptr(OpType::Rz, z_angle)});
      }
    }
    if (!multip_ry.empty()) {
      mutip_ry_vec.push_back(MultiplexedRotationBox(multip_ry));
    } else {
      mutip_ry_vec.push_back(std::nullopt);
    }
    if (!multip_rz.empty()) {
      mutip_rz_vec.push_back(MultiplexedRotationBox(multip_rz));
    } else {
      mutip_rz_vec.push_back(std::nullopt);
    }
    if (step == n_qubits - 1) {
      TKET_ASSERT(ts.size() == 1);
      double circ_phase = is_inverse ? -ts[0] : ts[0];
      circ.add_phase(circ_phase);
      break;
    }
    psi = rs.array() * (ts.array() * i_ * PI).exp();
  }
  TKET_ASSERT(mutip_rz_vec.size() == n_qubits);
  // add multiplexors;
  if (!is_inverse) {
    for (unsigned i = n_qubits; i-- > 0;) {
      std::vector<unsigned> args(n_qubits - i);
      std::iota(std::begin(args), std::end(args), 0);
      if (mutip_ry_vec[i] != std::nullopt) {
        circ.add_box(mutip_ry_vec[i].value(), args);
      }
      if (mutip_rz_vec[i] != std::nullopt) {
        circ.add_box(mutip_rz_vec[i].value(), args);
      }
    }
  } else {
    for (unsigned i = 0; i < n_qubits; i++) {
      std::vector<unsigned> args(n_qubits - i);
      std::iota(std::begin(args), std::end(args), 0);
      if (mutip_rz_vec[i] != std::nullopt) {
        circ.add_box(mutip_rz_vec[i].value(), args);
      }
      if (mutip_ry_vec[i] != std::nullopt) {
        circ.add_box(mutip_ry_vec[i].value(), args);
      }
    }
  }
  return circ;
};

void StatePreparationBox::generate_circuit() const {
  circ_ = std::make_shared<Circuit>(state_prep_circ(statevector_, is_inverse_));
}

nlohmann::json StatePreparationBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const StatePreparationBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["statevector"] = box.get_statevector();
  j["is_inverse"] = box.is_inverse();
  return j;
}

Op_ptr StatePreparationBox::from_json(const nlohmann::json &j) {
  StatePreparationBox box = StatePreparationBox(
      j.at("statevector").get<Eigen::VectorXcd>(),
      j.at("is_inverse").get<bool>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(StatePreparationBox, StatePreparationBox)

}  // namespace tket