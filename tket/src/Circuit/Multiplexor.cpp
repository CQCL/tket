// Copyright Quantinuum
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

#include "tket/Circuit/Multiplexor.hpp"

#include <complex>

#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/DiagonalBox.hpp"
#include "tket/Gate/GatePtr.hpp"
#include "tket/Gate/GateUnitaryMatrix.hpp"
#include "tket/Gate/Rotation.hpp"
#include "tket/Ops/OpJsonFactory.hpp"
#include "tket/Utils/Constants.hpp"
#include "tket/Utils/HelperFunctions.hpp"

namespace tket {

// to limit the time to decompose MultiplexedRotationBox and MultiplexedU2Box.
// can be relaxed in the future
static const unsigned MAX_N_CONTROLS = 32;

/**
 * @brief Implement a multiplexor by sequentially applying QControlBoxes.
 * Assume all ops have width n_targets, and all bitstrings have size n_controls
 * @param op_map
 * @param n_controls
 * @param n_targets
 * @return Circuit
 */
static Circuit multiplexor_sequential_decomp(
    const ctrl_op_map_t &op_map, unsigned n_controls, unsigned n_targets) {
  Circuit c(n_controls + n_targets);
  std::vector<unsigned> qubits(n_controls + n_targets);
  std::iota(std::begin(qubits), std::end(qubits), 0);
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    std::vector<unsigned> zero_ctrls;
    for (unsigned i = 0; i < it->first.size(); i++) {
      if (!it->first[i]) {
        zero_ctrls.push_back(i);
        c.add_op<unsigned>(OpType::X, {i});
      }
    }
    QControlBox qcbox(it->second, n_controls);
    c.add_box(qcbox, qubits);
    for (unsigned i : zero_ctrls) {
      c.add_op<unsigned>(OpType::X, {i});
    }
  }
  return c;
}

/**
 * @brief Indicates whether a recursion step in recursive_demultiplex_rotation
 * is either a left child, a right child, or the root.
 *
 */
enum class RecursionNodeType { Left = 0, Right = 1, Root = 2 };

/**
 * @brief Implement multiplexed rotation gate (i.e.uniformly controlled
 * same-axis rotations (UCR)) with 2^ctrl_qubits SQ rotations, 2^ctrl_qubits
 * CXs, and 2 H gates for X-axis rotations.
 *
 * https://arxiv.org/abs/quant-ph/0410066
 * This is a special case derived from equation (3)
 * A UCR gate controlled by n qubits have the decomposition UCR = CX P CX Q
 * (multiplication order), where P and Q are themselves UCR gates controlled by
 * n-1 qubits.
 *
 * Also notice that CX P CX Q = Q CX P CX, therefore we can control the
 * direction of each decomposition to avoid adding adjacent CX gates.
 * e.g. UCR = CX P CX Q = (CX Q' CX P' CX*) CX (CX* P'' CX Q'')
 * The two CX* can be cancelled,
 * hence UCR = CX P CX Q = (CX Q' CX P') CX (P''CX Q'')
 *
 *
 * @param angles list of 2^ctrl_qubits angles, angles[i] is the angle activated
 * by bitstring binary(i)
 * @param axis can be either Ry or Rz
 * @param total_qubits the total number of qubits in the final output circuit
 * @param circ circuit to update
 * @param node_type controls the decomposition of each demultiplex
 * step. Will be implemented as P CX Q if the step is a left child, and Q CX P
 * if the step is a right child. CX P CX Q will be implemented if the step is
 * the root.
 *
 */
static void recursive_demultiplex_rotation(
    const std::vector<Expr> &angles, const OpType &axis, unsigned total_qubits,
    std::vector<GateSpec> &commands, const RecursionNodeType &node_type) {
  unsigned n_rotations = angles.size();
  unsigned n_qubits = (unsigned)log2(n_rotations) + 1;
  unsigned mid = (unsigned)(n_rotations / 2);
  std::vector<Expr> p_angles;
  std::vector<Expr> q_angles;
  for (unsigned i = 0; i < mid; i++) {
    p_angles.push_back((angles[i] - angles[mid + i]) / 2);
    q_angles.push_back((angles[i] + angles[mid + i]) / 2);
  }
  // UCR = CX P CX Q = Q CX P CX
  // the left recursion child implements P CX Q, and the
  // right recursion child implements Q CX P to cancel the CXs
  if (node_type == RecursionNodeType::Right) {
    std::swap(p_angles, q_angles);
  }
  if (q_angles.size() == 1) {
    // base step
    commands.push_back(GateSpec(axis, q_angles[0]));
  } else {
    recursive_demultiplex_rotation(
        q_angles, axis, total_qubits, commands, RecursionNodeType::Left);
  }
  commands.push_back(GateSpec(OpType::CX, total_qubits - n_qubits));
  if (p_angles.size() == 1) {
    // base step
    commands.push_back(GateSpec(axis, p_angles[0]));
  } else {
    recursive_demultiplex_rotation(
        p_angles, axis, total_qubits, commands, RecursionNodeType::Right);
  }
  if (node_type == RecursionNodeType::Root) {
    // for the root step, we implement UCR = CX P CX Q
    commands.push_back(GateSpec(OpType::CX, total_qubits - n_qubits));
  }
}

/**
 * @brief Decompose diag(a,b) using eq(3)
 * Returns the matrices u and v, and the multiplexed Z rotation matrix
 * R defined by the rotation angles a0 and a1 (in half-turns) activated by 0 and
 * 1 respectively. The matrix D is fixed to ZZPhase(-0.5)
 *
 * @param a 2x2 unitary
 * @param b 2x2 unitary
 * @return std::tuple<Eigen::Matrix2cd, Eigen::Matrix2cd, double>
 */
static std::tuple<Eigen::Matrix2cd, Eigen::Matrix2cd, double, double>
constant_demultiplex(const Eigen::Matrix2cd &a, const Eigen::Matrix2cd &b) {
  Eigen::Matrix2cd X = a * b.adjoint();
  // decompose X using eq(19)
  std::vector<double> tk1_params = tk1_angles_from_unitary(X);
  Eigen::Matrix2cd tk1_su2 = get_matrix_from_tk1_angles(
      {tk1_params[0], tk1_params[1], tk1_params[2], 0.0});
  Complex x0 = tk1_su2(0, 0);
  double phi = tk1_params[3] * PI * 2;
  // compute r matrix using eq(11) and eq(12)
  double a0 = -PI / 2 - phi / 2 - std::arg(x0);
  double a1 = PI / 2 - phi / 2 + std::arg(x0);
  Complex r0 = std::exp(0.5 * i_ * a0);
  Complex r1 = std::exp(0.5 * i_ * a1);
  Eigen::Matrix2cd r = Eigen::Matrix2cd::Zero();
  r(0, 0) = r0;
  r(1, 1) = r1;
  // a0 and a1 defines the R matrix in eq(3). And they are the z rotation angles
  // activated by 0 and 1 respectively.
  a0 = a0 / PI;
  a1 = a1 / PI;
  Eigen::ComplexEigenSolver<Eigen::Matrix2cd> eigen_solver(r * X * r);
  Eigen::Matrix2cd u = eigen_solver.eigenvectors();
  // the eigenvalues are guaranteed to be {i, -i}. We permute u so its
  // eigenvalues have the order [i, -i].
  if (std::abs(eigen_solver.eigenvalues()[0] + i_) < EPS) {
    u.col(0).swap(u.col(1));
  }
  Eigen::Matrix2cd d = Eigen::Matrix2cd::Zero();
  d(0, 0) = std::sqrt(i_);
  d(1, 1) = std::sqrt(-i_);
  Eigen::Matrix2cd v = d * u.adjoint() * r.adjoint() * b;
  return {u, v, a0, a1};
}

/**
 * @brief Given the angles for a UCRz gate, return its block diagonal matrix
 * representation as a vector of 2x2 matrices
 *
 * @param angles
 * @return std::vector<Eigen::Matrix2cd>
 */
static std::vector<Eigen::Matrix2cd> ucrz_angles_to_diagonal(
    const std::vector<double> &angles) {
  std::vector<Eigen::Matrix2cd> diag;
  unsigned mid = (unsigned)(angles.size() / 2);
  for (unsigned i = 0; i < mid; i++) {
    Eigen::Matrix2cd u = Eigen::Matrix2cd::Zero();
    u(0, 0) = std::exp(-0.5 * angles[i] * i_ * PI);
    u(1, 1) = std::exp(-0.5 * angles[i + mid] * i_ * PI);
    diag.push_back(u);
  }
  for (unsigned i = 0; i < mid; i++) {
    Eigen::Matrix2cd u = Eigen::Matrix2cd::Zero();
    u(0, 0) = std::exp(0.5 * angles[i] * i_ * PI);
    u(1, 1) = std::exp(0.5 * angles[i + mid] * i_ * PI);
    diag.push_back(u);
  }
  return diag;
}

static void op_map_validate(const ctrl_op_map_t &op_map) {
  unsigned n_controls = 0;
  unsigned n_targets = 0;
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    op_signature_t op_sig = it->second->get_signature();
    if ((unsigned long)std::count(
            op_sig.begin(), op_sig.end(), EdgeType::Quantum) != op_sig.size()) {
      throw BadOpType(
          "Multiplexed operations cannot have classical wires.",
          it->second->get_type());
    }
    if (it == op_map.begin()) {
      n_controls = (unsigned)it->first.size();
      n_targets = (unsigned)op_sig.size();
    } else {
      if (it->first.size() != n_controls) {
        throw std::invalid_argument(
            "The bitstrings passed to the multiplexor must have the same "
            "width.");
      }
      if (op_sig.size() != n_targets) {
        throw std::invalid_argument(
            "Multiplexed operations must have the same width.");
      }
    }
  }
}

/**
 * @brief Recursively decompose a multiplexed U2 gate. (i.e. uniformly
 * controlled U2)
 *
 * Generates 2^ctrl_qubits Unitary1qBox, 2^ctrl_qubits CXs and a ladder of
 * MultiplexedRotationBoxes https://arxiv.org/abs/quant-ph/0410066 eq(3)
 *
 * During each recursion step, the multiplexor with n qubits defined
 * using `unitaries` are decomposed into
 * UCU = R (I tensor U) ZZPhase(-0.5, [0, n-1]) (I tensor V)
 * R is a UCRz gate, U and V are multiplexors.
 * Replace ZZPhase with CX and local gates we have
 * UCU = (R+1.5)(I tensor U TK1(0.5,0.5,0.5)) CX(0,n-1)(I tensor
 * TK1(0.5,0.5,0)V) and a 1.75 phase. R+1.5 means adding 1.5 to every Rz
 * rotations.
 *
 * At each subsequent step, the R gate can be merged with the multiplexor
 * on the left (in terms of matrix composition).
 * In the end, we will have a ladder of R gates at the end of the circuit, which
 * the user can decide whether to implement.
 *
 * @param unitaries list of 2^ctrl_qubits 2x2 unitaries, unitaries[i] is the
 * unitary activated by bitstring binary(i)
 * @param total_qubits the total number of qubits in the final output circuit
 * @param commands a vector holding the minimum command spec for each CX or U1
 * gate added
 * @param phase total global phase added by the method
 * @param ucrzs keep track of the ladder of MultiplexedRotationBoxes (R
 * gates). ucrzs[i] stores the Rz rotations angles (in half-turns) for the
 * MultiplexedRotationBox with i+2 qubits.
 * @param left_compose 2x2 unitary to be absorbed to left half
 * @param right_compose 2x2 unitary to be absorbed to right half
 */
static void recursive_demultiplex_u2(
    std::vector<Eigen::Matrix2cd> &unitaries, unsigned total_qubits,
    std::vector<GateSpec> &commands, float &phase,
    std::vector<std::vector<double>> &ucrzs,
    const Eigen::Matrix2cd &left_compose,
    const Eigen::Matrix2cd &right_compose) {
  // The following two constant matrices are the bottom SQ unitaries resulted
  // from decomposing the D gate (i.e. ZZPhase(-0.5)) using CX
  const static Eigen::Matrix2cd U_MULT =
      get_matrix_from_tk1_angles({0.5, 0.5, 0.5, 0.0});
  const static Eigen::Matrix2cd V_MULT =
      get_matrix_from_tk1_angles({0.5, 0.5, 0, 0.0});
  unsigned n_unitaries = unitaries.size();
  unsigned n_qubits = (unsigned)log2(n_unitaries) + 1;
  unsigned mid = (unsigned)(n_unitaries / 2);
  // We generalise eq(3) for n controls, demultiplex the multiplexor
  // by demultiplexing all pairs {unitaries[i], unitaries[mid+i]} 0<=i<mid.
  // i.e. I tensor diag(u) = I tensor diag(u_list)
  // I tensor diag(v) = I tensor diag(v_list)
  // D = ZZPhase(-0.5)
  // R = UCRz(rz_list, [q_{n-1}, q_{1}, q_{2}, ...,  q_{n-2}, q_0])
  std::vector<Eigen::Matrix2cd> u_list;
  std::vector<Eigen::Matrix2cd> v_list;
  std::vector<double> rz_list(n_unitaries);

  // merge previous UCRz gate into the multiplexor
  std::vector<Eigen::Matrix2cd> ucrz_diag =
      ucrz_angles_to_diagonal(ucrzs[n_qubits - 2]);
  for (unsigned i = 0; i < unitaries.size(); i++) {
    unitaries[i] = unitaries[i] * ucrz_diag[i];
  }
  // demultiplex pairs (unitaries[i], unitaries[mid+i])
  for (unsigned i = 0; i < mid; i++) {
    auto [u, v, a0, a1] =
        constant_demultiplex(unitaries[i], unitaries[mid + i]);
    u_list.push_back(u);
    v_list.push_back(v);
    rz_list[i] = a0;
    rz_list[i + mid] = a1;
  }

  // update the ucrzs with the 1.5 angle resulted from decomposing ZZPhase(-0.5)
  std::for_each(rz_list.begin(), rz_list.end(), [](double &f) { f += 1.5; });
  ucrzs[n_qubits - 2] = rz_list;

  // adding gates to the circuit
  // add v
  if (v_list.size() == 1) {
    Eigen::Matrix2cd v_prime = V_MULT * v_list[0] * left_compose;
    commands.push_back(GateSpec(OpType::U1, v_prime));
  } else {
    recursive_demultiplex_u2(
        v_list, total_qubits, commands, phase, ucrzs, left_compose, V_MULT);
  }
  // add CX
  commands.push_back(GateSpec(OpType::CX, total_qubits - n_qubits));

  phase += 1.75;
  // add u
  if (u_list.size() == 1) {
    Eigen::Matrix2cd u_prime = right_compose * u_list[0] * U_MULT;
    commands.push_back(GateSpec(OpType::U1, u_prime));
  } else {
    recursive_demultiplex_u2(
        u_list, total_qubits, commands, phase, ucrzs, U_MULT, right_compose);
  }
  return;
}

static ctrl_op_map_t op_map_symbol_sub(
    const SymEngine::map_basic_basic &sub_map, const ctrl_op_map_t &op_map) {
  ctrl_op_map_t new_op_map;
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    new_op_map.insert({it->first, it->second->symbol_substitution(sub_map)});
  }
  return new_op_map;
}

static SymSet op_map_free_symbols(const ctrl_op_map_t &op_map) {
  SymSet all_symbols;
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    SymSet op_symbols = it->second->free_symbols();
    all_symbols.insert(op_symbols.begin(), op_symbols.end());
  }
  return all_symbols;
}

static ctrl_op_map_t op_map_dagger(const ctrl_op_map_t &op_map) {
  ctrl_op_map_t new_op_map;
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    new_op_map.insert({it->first, it->second->dagger()});
  }
  return new_op_map;
}

static ctrl_op_map_t op_map_transpose(const ctrl_op_map_t &op_map) {
  ctrl_op_map_t new_op_map;
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    new_op_map.insert({it->first, it->second->transpose()});
  }
  return new_op_map;
}

static bool opmap_it_equal(
    const std::pair<std::vector<bool>, Op_ptr> &lhs,
    const std::pair<std::vector<bool>, Op_ptr> &rhs) {
  return lhs.first == rhs.first && *lhs.second == *rhs.second;
}
static bool tensored_opmap_it_equal(
    const std::pair<std::vector<bool>, std::vector<Op_ptr>> &lhs,
    const std::pair<std::vector<bool>, std::vector<Op_ptr>> &rhs) {
  return lhs.first == rhs.first &&
         std::equal(
             lhs.second.begin(), lhs.second.end(), rhs.second.begin(),
             rhs.second.end(),
             [](const Op_ptr &a, const Op_ptr &b) { return *a == *b; });
}

// Check if two ctrl_op_map_t are semantically equal.
static bool opmap_compare(
    const ctrl_op_map_t &map1, const ctrl_op_map_t &map2) {
  return std::equal(
      map1.begin(), map1.end(), map2.begin(), map2.end(), opmap_it_equal);
}
// Check if two ctrl_tensored_op_map_t are semantically equal.
static bool opmap_compare(
    const ctrl_tensored_op_map_t &map1, const ctrl_tensored_op_map_t &map2) {
  return std::equal(
      map1.begin(), map1.end(), map2.begin(), map2.end(),
      tensored_opmap_it_equal);
}

MultiplexorBox::MultiplexorBox(const ctrl_op_map_t &op_map)
    : Box(OpType::MultiplexorBox), op_map_(op_map) {
  auto it = op_map.begin();
  if (it == op_map.end()) {
    throw std::invalid_argument(
        "The op_map argument passed to MultiplexorBox cannot be empty.");
  }
  n_controls_ = (unsigned)it->first.size();
  n_targets_ = it->second->n_qubits();
  op_map_validate(op_map);
}

MultiplexorBox::MultiplexorBox(const MultiplexorBox &other)
    : Box(other),
      n_controls_(other.n_controls_),
      n_targets_(other.n_targets_),
      op_map_(other.op_map_) {}

Op_ptr MultiplexorBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  ctrl_op_map_t new_op_map = op_map_symbol_sub(sub_map, op_map_);
  return std::make_shared<MultiplexorBox>(new_op_map);
}

SymSet MultiplexorBox::free_symbols() const {
  return op_map_free_symbols(op_map_);
}

Op_ptr MultiplexorBox::dagger() const {
  return std::make_shared<MultiplexorBox>(op_map_dagger(op_map_));
}

Op_ptr MultiplexorBox::transpose() const {
  return std::make_shared<MultiplexorBox>(op_map_transpose(op_map_));
}

op_signature_t MultiplexorBox::get_signature() const {
  op_signature_t qubits(n_controls_ + n_targets_, EdgeType::Quantum);
  return qubits;
}

bool MultiplexorBox::is_equal(const Op &op_other) const {
  const MultiplexorBox &other = dynamic_cast<const MultiplexorBox &>(op_other);
  if (id_ == other.get_id()) return true;
  return opmap_compare(op_map_, other.op_map_);
}

nlohmann::json MultiplexorBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const MultiplexorBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["op_map"] = box.get_op_map();
  return j;
}

Op_ptr MultiplexorBox::from_json(const nlohmann::json &j) {
  MultiplexorBox box = MultiplexorBox(j.at("op_map").get<ctrl_op_map_t>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

void MultiplexorBox::generate_circuit() const {
  circ_ = std::make_shared<Circuit>(
      multiplexor_sequential_decomp(op_map_, n_controls_, n_targets_));
}

MultiplexedRotationBox::MultiplexedRotationBox(const ctrl_op_map_t &op_map)
    : Box(OpType::MultiplexedRotationBox), op_map_(op_map) {
  auto it = op_map.begin();
  if (it == op_map.end()) {
    throw std::invalid_argument(
        "The op_map argument passed to MultiplexedRotationBox cannot be "
        "empty.");
  }
  for (; it != op_map.end(); it++) {
    if (it == op_map.begin()) {
      n_controls_ = (unsigned)it->first.size();
      if (n_controls_ > MAX_N_CONTROLS) {
        throw std::invalid_argument(
            "MultiplexedRotationBox only supports bitstrings up to " +
            std::to_string(MAX_N_CONTROLS) + " bits.");
      }
      axis_ = it->second->get_type();
      if (axis_ != OpType::Rx && axis_ != OpType::Ry && axis_ != OpType::Rz) {
        throw BadOpType(
            "Ops passed to MultiplexedRotationBox must be either Rx, Ry, or "
            "Rz.",
            axis_);
      }
    } else {
      if (it->second->get_type() != axis_) {
        throw std::invalid_argument(
            "Ops passed to MultiplexedRotationBox must have the same rotation "
            "type.");
      }
    }
  }
  op_map_validate(op_map);
}

MultiplexedRotationBox::MultiplexedRotationBox(
    const MultiplexedRotationBox &other)
    : Box(other),
      n_controls_(other.n_controls_),
      op_map_(other.op_map_),
      axis_(other.axis_) {}

Op_ptr MultiplexedRotationBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  ctrl_op_map_t new_op_map = op_map_symbol_sub(sub_map, op_map_);
  return std::make_shared<MultiplexedRotationBox>(new_op_map);
}

SymSet MultiplexedRotationBox::free_symbols() const {
  return op_map_free_symbols(op_map_);
}

Op_ptr MultiplexedRotationBox::dagger() const {
  return std::make_shared<MultiplexedRotationBox>(op_map_dagger(op_map_));
}

Op_ptr MultiplexedRotationBox::transpose() const {
  return std::make_shared<MultiplexedRotationBox>(op_map_transpose(op_map_));
}

op_signature_t MultiplexedRotationBox::get_signature() const {
  op_signature_t qubits(n_controls_ + 1, EdgeType::Quantum);
  return qubits;
}

bool MultiplexedRotationBox::is_equal(const Op &op_other) const {
  const MultiplexedRotationBox &other =
      dynamic_cast<const MultiplexedRotationBox &>(op_other);
  if (id_ == other.get_id()) return true;
  return opmap_compare(op_map_, other.op_map_);
}

nlohmann::json MultiplexedRotationBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const MultiplexedRotationBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["op_map"] = box.get_op_map();
  return j;
}

Op_ptr MultiplexedRotationBox::from_json(const nlohmann::json &j) {
  MultiplexedRotationBox box =
      MultiplexedRotationBox(j.at("op_map").get<ctrl_op_map_t>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

std::vector<GateSpec> MultiplexedRotationBox::decompose() const {
  unsigned long long n_rotations = 1ULL << n_controls_;
  std::vector<Expr> rotations(n_rotations);
  // convert op_map to a vector of 2^n_controls_ angles
  for (unsigned long long i = 0; i < n_rotations; i++) {
    auto it = op_map_.find(dec_to_bin(i, n_controls_));
    if (it == op_map_.end()) {
      rotations[i] = 0;
    } else {
      rotations[i] = it->second->get_params()[0];
    }
  }
  std::vector<GateSpec> commands;
  OpType axis = axis_;
  if (axis_ == OpType::Rx) {
    commands.push_back(GateSpec(OpType::H, n_controls_));
    axis = OpType::Rz;
  }
  recursive_demultiplex_rotation(
      rotations, axis, n_controls_ + 1, commands, RecursionNodeType::Root);
  if (axis_ == OpType::Rx) {
    commands.push_back(GateSpec(OpType::H, n_controls_));
  }
  return commands;
}

void MultiplexedRotationBox::generate_circuit() const {
  Circuit circ(n_controls_ + 1);
  if (n_controls_ == 0) {
    auto it = op_map_.begin();
    circ.add_op<unsigned>(it->second, {0});
    circ_ = std::make_shared<Circuit>(circ);
    return;
  }

  for (const GateSpec &gs : this->decompose()) {
    switch (gs.type) {
      case OpType::CX:
        circ.add_op<unsigned>(OpType::CX, {*gs.qubit, n_controls_});
        break;
      case OpType::Rx:
        circ.add_op<unsigned>(OpType::Rx, *gs.angle, {n_controls_});
        break;
      case OpType::Ry:
        circ.add_op<unsigned>(OpType::Ry, *gs.angle, {n_controls_});
        break;
      case OpType::Rz:
        circ.add_op<unsigned>(OpType::Rz, *gs.angle, {n_controls_});
        break;
      case OpType::H:
        circ.add_op<unsigned>(OpType::H, {n_controls_});
        break;
      default:
        // this should never be hit
        TKET_ASSERT(false);
    }
  }

  circ_ = std::make_shared<Circuit>(circ);
}

MultiplexedU2Box::MultiplexedU2Box(const ctrl_op_map_t &op_map, bool impl_diag)
    : Box(OpType::MultiplexedU2Box), op_map_(op_map), impl_diag_(impl_diag) {
  auto it = op_map.begin();
  if (it == op_map.end()) {
    throw std::invalid_argument(
        "The op_map argument passed to MultiplexedU2Box cannot be empty.");
  }
  n_controls_ = (unsigned)it->first.size();
  if (n_controls_ > MAX_N_CONTROLS) {
    throw std::invalid_argument(
        "MultiplexedU2Box only supports bitstrings up to " +
        std::to_string(MAX_N_CONTROLS) + " bits.");
  }
  for (; it != op_map.end(); it++) {
    OpType optype = it->second->get_type();
    if (!is_single_qubit_unitary_type(optype) &&
        optype != OpType::Unitary1qBox) {
      throw BadOpType(
          "Ops passed to MultiplexedU2Box must be single-qubit unitary gate "
          "types or Unitary1qBox.",
          optype);
    }
  }
  op_map_validate(op_map);
}

MultiplexedU2Box::MultiplexedU2Box(const MultiplexedU2Box &other)
    : Box(other),
      n_controls_(other.n_controls_),
      op_map_(other.op_map_),
      impl_diag_(other.impl_diag_) {}

Op_ptr MultiplexedU2Box::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  ctrl_op_map_t new_op_map = op_map_symbol_sub(sub_map, op_map_);
  return std::make_shared<MultiplexedU2Box>(new_op_map, impl_diag_);
}

SymSet MultiplexedU2Box::free_symbols() const {
  return op_map_free_symbols(op_map_);
}

Op_ptr MultiplexedU2Box::dagger() const {
  return std::make_shared<MultiplexedU2Box>(op_map_dagger(op_map_), impl_diag_);
}

Op_ptr MultiplexedU2Box::transpose() const {
  return std::make_shared<MultiplexedU2Box>(
      op_map_transpose(op_map_), impl_diag_);
}

op_signature_t MultiplexedU2Box::get_signature() const {
  op_signature_t qubits(n_controls_ + 1, EdgeType::Quantum);
  return qubits;
}

bool MultiplexedU2Box::is_equal(const Op &op_other) const {
  const MultiplexedU2Box &other =
      dynamic_cast<const MultiplexedU2Box &>(op_other);
  if (id_ == other.get_id()) return true;
  return impl_diag_ == other.impl_diag_ &&
         opmap_compare(op_map_, other.op_map_);
}

nlohmann::json MultiplexedU2Box::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const MultiplexedU2Box &>(*op);
  nlohmann::json j = core_box_json(box);
  j["op_map"] = box.get_op_map();
  j["impl_diag"] = box.get_impl_diag();
  return j;
}

Op_ptr MultiplexedU2Box::from_json(const nlohmann::json &j) {
  MultiplexedU2Box box = MultiplexedU2Box(
      j.at("op_map").get<ctrl_op_map_t>(), j.at("impl_diag").get<bool>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

MultiplexedU2Commands MultiplexedU2Box::decompose() const {
  unsigned long long n_unitaries = 1ULL << n_controls_;
  std::vector<Eigen::Matrix2cd> unitaries(n_unitaries);
  // convert op_map to a vector of 2^n_controls_ unitaries
  for (unsigned long long i = 0; i < n_unitaries; i++) {
    auto it = op_map_.find(dec_to_bin(i, n_controls_));
    if (it == op_map_.end()) {
      unitaries[i] = Eigen::Matrix2cd::Identity();
    } else {
      if (it->second->get_type() == OpType::Unitary1qBox) {
        std::shared_ptr<const Unitary1qBox> u1box =
            std::dynamic_pointer_cast<const Unitary1qBox>(it->second);
        unitaries[i] = u1box->get_matrix();
      } else {
        if (!it->second->free_symbols().empty()) {
          throw Unsupported("Can't decompose symbolic MultiplexedU2Box.");
        }
        unitaries[i] = GateUnitaryMatrix::get_unitary(*as_gate_ptr(it->second));
      }
    }
  }

  // initialise the ucrz list
  std::vector<std::vector<double>> ucrzs(n_controls_);
  for (unsigned i = 0; i < n_controls_; i++) {
    ucrzs[i] = std::vector<double>(1ULL << (i + 1), 0.0);
  }

  std::vector<GateSpec> commands = {};
  float phase = 0;
  recursive_demultiplex_u2(
      unitaries, n_controls_ + 1, commands, phase, ucrzs,
      Eigen::Matrix2cd::Identity(), Eigen::Matrix2cd::Identity());
  // convert the ucrzs to a diagonal matrix
  Eigen::VectorXcd diag =
      Eigen::VectorXcd::Constant(1ULL << (n_controls_ + 1), 1);
  for (unsigned i = 0; i < n_controls_; i++) {
    // ith ucrzs acts on i+2 qubits
    // which has n_controls_ + 1 - (i+2) identities its the tensor product
    // therefore (n_controls_ + 1 - (i+2))^2 copies in the diagonal
    for (unsigned long long offset = 0;
         offset < (1ULL << (n_controls_ + 1 - (i + 2))); offset++) {
      for (unsigned long long j = 0; j < (1ULL << (i + 1)); j++) {
        // the bitstrings in a ucrz are mapped to qubits not in the standard
        // order
        unsigned long long diag_idx =
            (j >= (1ULL << i)) ? (j - (1ULL << i)) * 2 + 1 : j * 2;
        diag[diag_idx + offset * (1ULL << (i + 2))] *=
            std::exp(-0.5 * i_ * PI * ucrzs[i][j]);
        diag[diag_idx + offset * (1ULL << (i + 2)) + (1ULL << (i + 1))] *=
            std::exp(0.5 * i_ * PI * ucrzs[i][j]);
      }
    }
  }
  return MultiplexedU2Commands(commands, diag, phase);
}

void MultiplexedU2Box::generate_circuit() const {
  Circuit circ(n_controls_ + 1);
  Eigen::VectorXcd diag_vec;

  if (n_controls_ == 0) {
    auto it = op_map_.begin();
    circ.add_op<unsigned>(it->second, {0});
    circ_ = std::make_shared<Circuit>(circ);
    return;
  }

  MultiplexedU2Commands decomp = this->decompose();
  for (unsigned i = 0; i < decomp.commands.size(); i++) {
    GateSpec gc = decomp.commands[i];
    // n.b. with zero indexing "n_controls" corresponds to the target qubit
    switch (gc.type) {
      case OpType::CX:
        circ.add_op<unsigned>(OpType::CX, {*gc.qubit, n_controls_});
        break;
      case OpType::U1:
        circ.add_box(Unitary1qBox(*gc.matrix), {n_controls_});
        break;
      default:
        // this should never be hit
        TKET_ASSERT(false);
    }
  }

  circ.add_phase(decomp.phase);

  if (impl_diag_ &&
      (decomp.diag - Eigen::VectorXcd::Constant(1ULL << circ.n_qubits(), 1))
              .cwiseAbs()
              .sum() > EPS) {
    std::vector<unsigned> args(circ.n_qubits());
    std::iota(std::begin(args), std::end(args), 0);
    circ.add_box(DiagonalBox(decomp.diag), args);
  }
  circ_ = std::make_shared<Circuit>(circ);
}

MultiplexedTensoredU2Box::MultiplexedTensoredU2Box(
    const ctrl_tensored_op_map_t &op_map)
    : Box(OpType::MultiplexedTensoredU2Box), op_map_(op_map) {
  auto it = op_map.begin();
  if (it == op_map.end()) {
    throw std::invalid_argument(
        "The op_map argument passed to MultiplexedTensoredU2Box cannot be "
        "empty.");
  }
  n_controls_ = (unsigned)it->first.size();
  n_targets_ = (unsigned)it->second.size();
  if (n_controls_ > MAX_N_CONTROLS) {
    throw std::invalid_argument(
        "MultiplexedTensoredU2Box only supports bitstrings up to " +
        std::to_string(MAX_N_CONTROLS) + " bits.");
  }
  for (; it != op_map.end(); it++) {
    if (it->first.size() != n_controls_) {
      throw std::invalid_argument(
          "The bitstrings passed to MultiplexedTensoredU2Box must have the "
          "same width.");
      ;
    }
    if (it->second.size() != n_targets_) {
      throw std::invalid_argument(
          "Each tensored operation passed to MultiplexedTensoredU2Box must "
          "have the same number of U2 components");
    }
    for (auto op : it->second) {
      OpType optype = op->get_type();
      if (!is_single_qubit_unitary_type(optype) &&
          optype != OpType::Unitary1qBox) {
        throw BadOpType(
            "Ops passed to MultiplexedTensoredU2Box must be single-qubit "
            "unitary gate types or Unitary1qBox.",
            optype);
      }
    }
  }
}

MultiplexedTensoredU2Box::MultiplexedTensoredU2Box(
    const MultiplexedTensoredU2Box &other)
    : Box(other),
      n_controls_(other.n_controls_),
      n_targets_(other.n_targets_),
      op_map_(other.op_map_) {}

Op_ptr MultiplexedTensoredU2Box::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  ctrl_tensored_op_map_t new_op_map;
  for (auto it = op_map_.begin(); it != op_map_.end(); it++) {
    std::vector<Op_ptr> ops;
    for (auto op : it->second) {
      ops.push_back(op->symbol_substitution(sub_map));
    }
    new_op_map.insert({it->first, ops});
  }
  return std::make_shared<MultiplexedTensoredU2Box>(new_op_map);
}

SymSet MultiplexedTensoredU2Box::free_symbols() const {
  SymSet all_symbols;
  for (auto it = op_map_.begin(); it != op_map_.end(); it++) {
    for (auto op : it->second) {
      SymSet op_symbols = op->free_symbols();
      all_symbols.insert(op_symbols.begin(), op_symbols.end());
    }
  }
  return all_symbols;
}

Op_ptr MultiplexedTensoredU2Box::dagger() const {
  ctrl_tensored_op_map_t new_op_map;
  for (auto it = op_map_.begin(); it != op_map_.end(); it++) {
    std::vector<Op_ptr> ops;
    for (auto op : it->second) {
      ops.push_back(op->dagger());
    }
    new_op_map.insert({it->first, ops});
  }
  return std::make_shared<MultiplexedTensoredU2Box>(new_op_map);
}

Op_ptr MultiplexedTensoredU2Box::transpose() const {
  ctrl_tensored_op_map_t new_op_map;
  for (auto it = op_map_.begin(); it != op_map_.end(); it++) {
    std::vector<Op_ptr> ops;
    for (auto op : it->second) {
      ops.push_back(op->transpose());
    }
    new_op_map.insert({it->first, ops});
  }
  return std::make_shared<MultiplexedTensoredU2Box>(new_op_map);
}

op_signature_t MultiplexedTensoredU2Box::get_signature() const {
  op_signature_t qubits(n_controls_ + n_targets_, EdgeType::Quantum);
  return qubits;
}

bool MultiplexedTensoredU2Box::is_equal(const Op &op_other) const {
  const MultiplexedTensoredU2Box &other =
      dynamic_cast<const MultiplexedTensoredU2Box &>(op_other);
  if (id_ == other.get_id()) return true;
  return opmap_compare(op_map_, other.op_map_);
}

nlohmann::json MultiplexedTensoredU2Box::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const MultiplexedTensoredU2Box &>(*op);
  nlohmann::json j = core_box_json(box);
  j["op_map"] = box.get_op_map();
  return j;
}

Op_ptr MultiplexedTensoredU2Box::from_json(const nlohmann::json &j) {
  MultiplexedTensoredU2Box box =
      MultiplexedTensoredU2Box(j.at("op_map").get<ctrl_tensored_op_map_t>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

void add_cx_u1(
    Circuit &circ, const std::vector<MultiplexedU2Commands> &m_u2_decomps,
    unsigned n_controls_, unsigned n_targets_) {
  TKET_ASSERT(m_u2_decomps.size() == n_targets_);
  // Each Multiplexor Decomposition will be up to some global phase
  // First we add these phase contributions to the circuit
  unsigned reference_size = m_u2_decomps[0].commands.size();
  for (unsigned i = 0; i < m_u2_decomps.size(); i++) {
    circ.add_phase(m_u2_decomps[i].phase);
    // We also confirm that each Multiplexor decomposition has the
    // same number of commands
    TKET_ASSERT(reference_size == m_u2_decomps[i].commands.size());
  }

  // we now iterate through all the commands, adding them to the circuit
  // in an interleaved manner
  for (unsigned i = 0; i < reference_size; i++) {
    for (unsigned target = 0; target < m_u2_decomps.size(); target++) {
      GateSpec gate = m_u2_decomps[target].commands[i];
      unsigned rotated_index;
      switch (gate.type) {
        case OpType::CX:
          // we also need to map gate.qubit to the correct qubit
          // we know that the bitstrings for the "target"th target have been
          // left rotated by "target", so:
          rotated_index = (*gate.qubit + (target % n_targets_)) % n_controls_;
          TKET_ASSERT(i % 2 == 1);
          circ.add_op<unsigned>(
              OpType::CX, {rotated_index, n_controls_ + target});
          break;
        case OpType::U1:
          TKET_ASSERT(i % 2 == 0);
          circ.add_box(Unitary1qBox(*gate.matrix), {n_controls_ + target});
          break;
        default:
          // this should never be hit
          TKET_ASSERT(false);
      }
    }
  }
  return;
}

std::pair<ctrl_op_map_t, Eigen::VectorXcd>
disentangle_final_qubit_from_diagonal(
    const Eigen::VectorXcd &full_diag, unsigned n_controls_) {
  // disentangle one qubit from the diagonal
  // results in a multiplexed-Rz targeting the target j
  Eigen::VectorXcd diag_vec =
      Eigen::VectorXcd::Constant(1ULL << n_controls_, 1);
  ctrl_op_map_t multip_rz;
  for (unsigned long long j = 0; j < (1ULL << n_controls_); j++) {
    // As full_diag as produced by MultiplexedU2Box::decompose
    // adds an extra qubit for the compensating diagonal, we know
    // that the control bitstring corresponding to j should match
    // up with the rotated bitstrings
    Complex a = full_diag[2 * j];
    Complex b = full_diag[2 * j + 1];
    // convert diag[a,b] into a p*Rz(alpha)
    double a_phase = std::arg(a);
    double b_phase = std::arg(b);
    double alpha = (b_phase - a_phase) / PI;
    Complex p = std::exp((b_phase + a_phase) * 0.5 * i_);
    if (std::abs(alpha) > EPS) {
      // bitstr should innately correspond to the rotated controls
      // when being constructed from diagonal provided
      // by MultiplexedU2Box::decompose
      multip_rz.insert(
          {dec_to_bin(j, n_controls_), get_op_ptr(OpType::Rz, alpha)});
    }
    diag_vec[j] *= p;
  }
  return std::make_pair(multip_rz, diag_vec);
}

void add_multi_rz(
    Circuit &circ, const std::vector<ctrl_op_map_t> &all_multiplexed_rz,
    unsigned n_controls_, unsigned n_targets_) {
  TKET_ASSERT(all_multiplexed_rz.size() == n_targets_);
  std::vector<std::vector<GateSpec>> all_decomps;
  // First get all GateSpec by constructing and decomposing
  // MultiplexedRotationBox
  for (unsigned target = 0; target < n_targets_; target++) {
    ctrl_op_map_t map = all_multiplexed_rz[target];
    if (!map.empty()) {
      all_decomps.push_back(
          MultiplexedRotationBox(all_multiplexed_rz[target]).decompose());
    } else {
      all_decomps.push_back({});
    }
  }
  TKET_ASSERT(!all_decomps.empty());
  unsigned reference_size = 0;
  for (unsigned i = 0; i < all_decomps.size(); i++) {
    if (all_decomps[i].empty()) continue;
    if (!reference_size) reference_size = all_decomps[i].size();
    TKET_ASSERT(reference_size == all_decomps[i].size());
  }

  // => no MultiplexedRz so we can carry on
  if (!reference_size) return;

  // Then iterate through all the commands, adding them to the circuit
  // in an interleaved manner
  for (unsigned i = 0; i < reference_size; i++) {
    for (unsigned target = 0; target < all_decomps.size(); target++) {
      if (!all_decomps[target].empty()) {
        GateSpec gate = all_decomps[target][i];
        unsigned rotated_index;
        switch (gate.type) {
          case OpType::CX:
            // we also need to map gate.qubit to the correct qubit
            // we know that the bitstrings for the "target"th target have been
            // left rotated by "target", so:
            rotated_index = (*gate.qubit + (target % n_targets_)) % n_controls_;
            circ.add_op<unsigned>(
                OpType::CX, {rotated_index, n_controls_ + target});
            break;
          case OpType::Rz:
            circ.add_op<unsigned>(
                OpType::Rz, *gate.angle, {n_controls_ + target});
            break;
          default:
            // this should never be hit
            TKET_ASSERT(false);
        }
      }
    }
  }
  return;
}

Eigen::VectorXcd combine_diagonals(
    const std::vector<Eigen::VectorXcd> &all_diags, unsigned n_controls_,
    unsigned n_targets_) {
  Eigen::VectorXcd combined_diag_vec =
      Eigen::VectorXcd::Constant(1ULL << n_controls_, 1);
  TKET_ASSERT(all_diags.size() == n_targets_);
  for (unsigned rotate = 0; rotate < n_targets_; rotate++) {
    Eigen::VectorXcd diag_vec = all_diags[rotate];
    TKET_ASSERT(diag_vec.size() == combined_diag_vec.size());
    // the "rotate" indexed diagonal vector in all_diags
    // will have indexing corresponding to a left rotation
    // of the input bit strings of "rotate"
    for (unsigned index = 0; index < diag_vec.size(); index++) {
      // to construct the diagonal vector correctly, we take
      // the value "index", convert it to a bitstring, right
      // rotate it by "rotate" and convert it back an integer.
      unsigned rotate_value = rotate % n_controls_;
      std::vector<bool> as_bits = dec_to_bin(index, n_controls_);
      // right rotate
      std::rotate(
          as_bits.begin(), as_bits.begin() + (as_bits.size() - rotate_value),
          as_bits.end());
      unsigned rotated_index = bin_to_dec(as_bits);
      TKET_ASSERT(rotated_index <= combined_diag_vec.size());
      // This gives a new index for updating the correct element of the
      // diagonal vector
      combined_diag_vec[rotated_index] *= diag_vec[index];
    }
  }
  return combined_diag_vec;
}

void MultiplexedTensoredU2Box::generate_circuit() const {
  /**
   * We will generate the circuit for MultiplexedTensoredU2Box by
   * separating it into Individual Multiplexed-U2 that are synthesised
   * separately. We will also reduce on depth by interleaving the {CX, U1} and
   * {CX, Rz} subcircuits of each Multiplexed-U2. We will save on gate count by
   * merging parts of the Diagonal gate produced by the Multiplexed-U2.
   */

  // Break the input into separate Multiplexors
  std::vector<MultiplexedU2Commands> m_u2_decomps;
  for (unsigned target = 0; target < n_targets_; target++) {
    ctrl_op_map_t u2_op_map;
    for (auto it = op_map_.begin(); it != op_map_.end(); it++) {
      // by rotating the control condition we change the order of CX gates for
      // each decomposition we can later use this to interleave the multiplexor
      // decompositions and so reduce depth
      std::vector<bool> control_condition = it->first;
      std::rotate(
          control_condition.begin(),
          control_condition.begin() + (target % n_controls_),
          control_condition.end());
      u2_op_map.insert({control_condition, it->second[target]});
    }
    m_u2_decomps.push_back(MultiplexedU2Box(u2_op_map).decompose());
  }

  // Next we split each diagonal vector for each MultiplexedU2 into a
  // Multiplexed-Rz on the target qubit and a Diagonal gate over the control
  // register
  std::vector<ctrl_op_map_t> all_multiplexed_rz;
  std::vector<Eigen::VectorXcd> all_diags;
  for (unsigned i = 0; i < m_u2_decomps.size(); i++) {
    std::pair<ctrl_op_map_t, Eigen::VectorXcd> disentangled =
        disentangle_final_qubit_from_diagonal(
            m_u2_decomps[i].diag, n_controls_);
    all_multiplexed_rz.push_back(disentangled.first);
    all_diags.push_back(disentangled.second);
  }

  // Final we merge the diagonals over the same qubits into a combined operator
  Eigen::VectorXcd combined_diag_vec =
      combine_diagonals(all_diags, n_controls_, n_targets_);

  // Now we can construct the circuit - first we add the U1 + CX segment of the
  // circuit construction with interleaving
  Circuit circ(n_controls_ + n_targets_);
  add_cx_u1(circ, m_u2_decomps, n_controls_, n_targets_);

  // Then add Multiplexed-Rz gates to circ
  add_multi_rz(circ, all_multiplexed_rz, n_controls_, n_targets_);

  // Finally add the combined diagonal vector to the circuit
  if ((combined_diag_vec - Eigen::VectorXcd::Constant(1ULL << n_controls_, 1))
          .cwiseAbs()
          .sum() > EPS) {
    std::vector<unsigned> control_qubits(n_controls_);
    std::iota(std::begin(control_qubits), std::end(control_qubits), 0);
    circ.add_box(DiagonalBox(combined_diag_vec), control_qubits);
  }
  circ_ = std::make_shared<Circuit>(circ);
}

REGISTER_OPFACTORY(MultiplexorBox, MultiplexorBox)
REGISTER_OPFACTORY(MultiplexedRotationBox, MultiplexedRotationBox)
REGISTER_OPFACTORY(MultiplexedU2Box, MultiplexedU2Box)
REGISTER_OPFACTORY(MultiplexedTensoredU2Box, MultiplexedTensoredU2Box)

}  // namespace tket