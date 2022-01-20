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

#include <optional>

#include "Circuit/Circuit.hpp"
#include "Gate/GatePtr.hpp"

namespace tket {

/**
 * @brief Abstract Squasher interface
 *
 * Subclasses must implement these methods to be used as a squasher
 * in SingleQubitSquash
 */
class AbstractSquasher {
 public:
  /**
   * @brief Whether the OpType can be added to current squash.
   *
   * @param ot OpType to be accepted.
   * @retval true `ot` can be squashed.
   * @retval false `ot` cannot be squashed.
   */
  virtual bool accepts(OpType ot) const = 0;

  /**
   * @brief Add gate to current squash.
   *
   * @pre `accepts(gp->get_type())` must return true
   * @param gp Gate to add to squash.
   */
  virtual void append(Gate_ptr gp) = 0;

  /**
   * @brief obtain the current squash as circuit and gate to be commuted
   *
   * Optionally use the commutation colour of the next gate to return an
   * additional Gate_ptr to be commuted through.
   * If no gate should be commuted, return nullptr.
   * If `commutation_colour==std::nullopt`, then the returned Gate_ptr
   * is expected to be `nullptr`
   *
   * @param commutation_colour
   * @return std::pair<Circuit, Gate_ptr>
   */
  virtual std::pair<Circuit, Gate_ptr> flush(
      std::optional<Pauli> commutation_colour) const = 0;

  /**
   * @brief Reset the current squash.
   */
  virtual void clear() = 0;

  virtual ~AbstractSquasher() {}
};

/**
 * @brief Squashes single qubit gates using given Squasher.
 */
class SingleQubitSquash {
 private:
  using Condition = std::optional<std::pair<std::list<VertPort>, unsigned>>;

 public:
  /**
   * @brief Construct a new Single Qubit Squash object.
   *
   * @param squasher The Squasher instance.
   * @param reversed whether squashing is made back-to-front or front to back.
   */
  SingleQubitSquash(
      std::unique_ptr<AbstractSquasher> squasher, bool reversed = false)
      : squasher_(std::move(squasher)),
        reversed_(reversed),
        circ_ptr_(nullptr) {}

  /**
   * @brief Squash an entire circuit, one qubit at a time.
   *
   * @param circ The circuit to be squashed.
   * @return true The squash succeeded.
   * @return false The circuit was not changed.
   */
  bool squash(Circuit &circ);

  /**
   * @brief Squash everything between in-edge and out-edge
   *
   * If `reversed=true`, then the in-edge should come after the out-edge
   * in the circuit.
   *
   * @param in Starting edge of the squash.
   * @param out Last edge of the squash.
   * @return true The squash succeeded.
   * @return false The circuit was not changed.
   */
  bool squash_between(const Edge &in, const Edge &out);

 private:
  std::unique_ptr<AbstractSquasher> squasher_;
  bool reversed_;
  // points to the current circuit during squashing
  Circuit *circ_ptr_;

  // substitute chain by a sub circuit, handling conditions
  // and backing up + restoring current edge
  void substitute(
      const Circuit &sub, const VertexVec &single_chain, Edge &e,
      const Condition &condition);

  // insert a gate at the given edge, respecting condition
  void insert_left_over_gate(
      Op_ptr left_over, const Edge &e, const Condition &condition);

  // whether a vertex can be squashed with the previous vertices
  bool is_squashable(Vertex v, OpType v_type) const;

  // whether the sub circuit is shorter than chain
  bool sub_is_better(
      const Circuit &sub, const std::vector<Gate_ptr> chain) const;

  // returns a description of the condition of current vertex
  Condition get_condition(Vertex v) const;

  // simple utils respecting reversed boolean
  Vertex next_vertex(const Edge &e) const;

  port_t next_port(const Edge &e) const;

  Edge prev_edge(const VertPort &pair) const;

  Edge next_edge(const Vertex &v, const Edge &e) const;

  bool is_last_optype(OpType type) const;

  static bool is_equal(
      const Circuit &circ, const std::vector<Gate_ptr> &gates,
      bool reversed = false);
};

}  // namespace tket
