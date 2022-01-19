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

#pragma once

#include "OpType/OpType.hpp"
#include "Utils/Expression.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket {

/** A faithful representation of SU(2). */
class Rotation {
 public:
  /** Identity */
  Rotation()
      : rep_(Rep::id), s_(1), i_(0), j_(0), k_(0), optype_(OpType::noop) {}

  /**
   * Represent an X, Y or Z rotation
   *
   * @param optype one of @ref OpType::Rx, @ref OpType::Ry or @ref OpType::Rz
   * @param a angle in half-turns
   */
  Rotation(OpType optype, Expr a);

  /** Apply a second rotation */
  void apply(const Rotation& other);

  /** Is it the identity? */
  bool is_id() const { return rep_ == Rep::id; }

  /** Is it minus the identity? */
  bool is_minus_id() const { return rep_ == Rep::minus_id; }

  /**
   * Return the angle given the axis
   *
   * @param optype axis of rotation
   *
   * @return angle of rotation in half-turns, if the axis matches
   * @pre \p optype is @ref OpType::Rx, @ref OpType::Ry or @ref OpType::Rz
   */
  std::optional<Expr> angle(OpType optype) const;

  /**
   * Convert to a sequence of angles in a PQP representation
   *
   * @param p one of @ref OpType::Rx, @ref OpType::Ry or @ref OpType::Rz
   * @param q one of @ref OpType::Rx, @ref OpType::Ry or @ref OpType::Rz
   *
   * @return angles (p1,q,p2) (in half-turns) such that the rotation is
   *         equivalent to P(p1) followed by Q(q) followed by P(p2)
   *
   * @pre p != q
   */
  std::tuple<Expr, Expr, Expr> to_pqp(OpType p, OpType q) const;

  // Default copy and assignment is fine.

  friend std::ostream& operator<<(std::ostream& os, const Rotation& q);

 private:
  enum class Rep {
    id,       /**< identity rotation */
    minus_id, /**< minus the identity  */
    orth_rot, /**< rotation about X, Y or Z */
    quat      /**< general rotation */
  };
  Rep rep_;

  // We represent every rotation as a quaternion of unit norm:
  Expr s_;  // scalar
  Expr i_;  // i coordinate
  Expr j_;  // j coordinate
  Expr k_;  // k coordinate

  // If rep_ == Rep::orth_rot, we represent the rotation as an axis and angle:
  OpType optype_;
  Expr a_;
};

/**
 * Construct TK1 angles and phase from matrix
 *
 * @param U 2x2 unitary matrix
 *
 * @return [a,b,c,t] where a,b,c are the TK1 angles and t is the phase
 */
std::vector<double> tk1_angles_from_unitary(const Eigen::Matrix2cd& U);

/**
 * Construct matrix from TK1 angles and phase
 *
 * @param params [a,b,c,t] where a,b,c are the TK1 angles and t is the phase
 *
 * @return 2x2 unitary matrix
 * @pre no symbolic parameters
 */
Eigen::Matrix2cd get_matrix_from_tk1_angles(std::vector<Expr> params);

}  // namespace tket
