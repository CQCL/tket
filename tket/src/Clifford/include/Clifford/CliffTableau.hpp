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
#include "Utils/BiMapHeaders.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "Utils/PauliStrings.hpp"

namespace tket {

class Circuit;

class CliffTableau {
  /*
   * Based off of the stabiliser-destabiliser tableau in Aaronson & Gottesman,
   * "Improved Simulation of Stabilizer Circuits",
   * https://arxiv.org/pdf/quant-ph/0406196.pdf
   *
   * This paper uses the rows of the tableau to describe the generators for the
   * set of stabilizers/destabilizers of a state. They then start with the state
   * |0>^n and perform column operations for each gate to show how all
   * generators are affected by the interacting qubits (columns).
   *
   * We flip this and instead have rows representing how the output of each wire
   * is affected by the inputs. The ZPauli on a given wire is the Pauli operator
   * that would be applied if an Rz gate were applied there and commuted through
   * to the inputs, similarly for the XPauli and Rx. This is the same as the
   * operator one would measure on the inputs if a standard/Hadamard basis
   * measurement were applied on the wire. Hence applying gates performs row
   * operations on the interacting wires.
   *
   * We can also consider the effect of adding a gate at the start of our
   * circuit (bringing it from the environment into the box) by column
   * operations.
   */
 public:
  /**
   * Construct the tableau for the identity over n qubits (given default qubit
   * names).
   */
  explicit CliffTableau(unsigned n);

  /**
   * Construct the tableau for the identity over specific qubits.
   */
  explicit CliffTableau(const qubit_vector_t &qbs);

  /**
   * Copy constructor
   */
  CliffTableau(const CliffTableau &other)
      : size_(other.size_),
        xpauli_x(other.xpauli_x),
        xpauli_z(other.xpauli_z),
        xpauli_phase(other.xpauli_phase),
        zpauli_x(other.zpauli_x),
        zpauli_z(other.zpauli_z),
        zpauli_phase(other.zpauli_phase),
        qubits_(other.qubits_){};

  /**
   * Access the Pauli string on the Z-channel of a given qubit with respect to
   * all inputs.
   */
  QubitPauliTensor get_zpauli(const Qubit &qb) const;

  /**
   * Access the Pauli string on the X-channel of a given qubit with respect to
   * all inputs.
   */
  QubitPauliTensor get_xpauli(const Qubit &qb) const;

  /**
   * Access all IDs for the qubits captured by the tableau.
   */
  std::set<Qubit> get_qubits() const;

  /**
   * Transform the tableau according to consuming a Clifford gate at either
   * end of the circuit.
   * Args are indices in the tableau's qubit ordering, so may not necessarily
   * match up to register indices.
   */
  void apply_S_at_front(unsigned qb);
  void apply_S_at_end(unsigned qb);
  void apply_V_at_front(unsigned qb);
  void apply_V_at_end(unsigned qb);
  void apply_CX_at_front(unsigned control, unsigned target);
  void apply_CX_at_end(unsigned control, unsigned target);

  /**
   * Transform the tableau according to consuming a Clifford gate at either
   * end of the circuit.
   */
  void apply_gate_at_front(OpType type, const std::vector<unsigned> &qbs);
  void apply_gate_at_front(OpType type, const qubit_vector_t &qbs);
  void apply_gate_at_end(OpType type, const std::vector<unsigned> &qbs);
  void apply_gate_at_end(OpType type, const qubit_vector_t &qbs);

  /**
   * Transform the tableau according to consuming a Clifford-phase pauli
   * gadget at either end of the circuit.
   *
   * @param pauli The string of the pauli gadget
   * @param half_pis The Clifford angle: {0, 1, 2, 3} represents {0, pi/2, pi,
   * -pi/2}
   */
  void apply_pauli_at_front(const QubitPauliTensor &pauli, unsigned half_pis);
  void apply_pauli_at_end(const QubitPauliTensor &pauli, unsigned half_pis);

  /**
   * Combine two tableaus in sequence.
   * Will throw an exception if the tableaus are not over the same set of
   * qubits.
   *
   * @param first first circuit
   * @param second second circuit
   *
   * @return The tableau corresponding to applying \p first, followed by \p
   * second
   */
  static CliffTableau compose(
      const CliffTableau &first, const CliffTableau &second);

  friend CliffTableau circuit_to_tableau(const Circuit &circ);
  friend Circuit tableau_to_circuit(const CliffTableau &tab);

  friend std::ostream &operator<<(std::ostream &os, CliffTableau const &tab);
  bool operator==(const CliffTableau &other) const;

 private:
  /** Number of qubits */
  const unsigned size_;

  /**
   * (A)pauli_(B): matrix representing the dependence of the A-channel of each
   * output on the (B) channel of each input
   *
   * (A)pauli_phase: whether or not there is a phase-flip on the A-channel of
   * each output
   */
  MatrixXb xpauli_x;
  MatrixXb xpauli_z;
  VectorXb xpauli_phase;
  MatrixXb zpauli_x;
  MatrixXb zpauli_z;
  VectorXb zpauli_phase;

  /** Map from qubit IDs to their row/column index in tableau */
  boost::bimap<Qubit, unsigned> qubits_;

  /**
   * Boolean encoding of Pauli
   * <x, z> = <false, false> ==> I
   * <x, z> = <false, true>  ==> Z
   * <x, z> = <true,  false> ==> X
   * <x, z> = <true,  true>  ==> Y
   */
  struct BoolPauli {
    bool x;
    bool z;

    /**
     * Lexicographic ordering by <x, z>
     */
    bool operator<(const BoolPauli &other) const {
      if (x == other.x) {
        return z < other.z;
      }
      return x < other.x;
    }
  };

  /**
   * Look-up table for Pauli multiplication with boolean encoding
   */
  static const std::map<
      std::pair<BoolPauli, BoolPauli>, std::pair<BoolPauli, Complex>>
      mult_lut;

  /**
   * Helper methods for manipulating the tableau when applying gates
   */
  void row_mult(
      const MatrixXb::RowXpr &xa, const MatrixXb::RowXpr &za, const bool &pa,
      const MatrixXb::RowXpr &xb, const MatrixXb::RowXpr &zb, const bool &pb,
      Complex phase, MatrixXb::RowXpr &xw, MatrixXb::RowXpr &zw, bool &pw);
  void col_mult(
      const MatrixXb::ColXpr &a, const MatrixXb::ColXpr &b, bool flip,
      MatrixXb::ColXpr &w, VectorXb &pw);
};

std::ostream &operator<<(std::ostream &os, CliffTableau const &tab);

}  // namespace tket
