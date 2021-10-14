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

#ifndef _TKET_Transform_H_
#define _TKET_Transform_H_

#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Architecture/Architectures.hpp"
#include "Characterisation/DeviceCharacterisation.hpp"
#include "Circuit/Circuit.hpp"
#include "Utils/Json.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "Utils/PauliStrings.hpp"

namespace tket {

class ControlDecompError : public std::logic_error {
 public:
  explicit ControlDecompError(const std::string& message)
      : std::logic_error(message) {}
};

/* Dictates whether synthesis of a PauliGraph should
    be done on the Paulis individually, making use of the pairwise
    interactions or collecting into mutually commuting sets. */
enum class PauliSynthStrat { Individual, Pairwise, Sets };

NLOHMANN_JSON_SERIALIZE_ENUM(
    PauliSynthStrat, {{PauliSynthStrat::Individual, "Individual"},
                      {PauliSynthStrat::Pairwise, "Pairwise"},
                      {PauliSynthStrat::Sets, "Sets"}});

class Transform {
 public:
  typedef std::function<bool(Circuit&)> Transformation;
  typedef std::function<unsigned(const Circuit&)> Metric;

  // the actual transformation to be applied
  // performs transformation in place and returns true iff made some change
  Transformation apply;  // this would ideally be `const`, but that deletes the
                         // copy assignment operator for Transform.

  explicit Transform(const Transformation& trans) : apply(trans) {}

  static const Transform id;  // identity Transform (does nothing to Circuit)

  ///////////////
  // Combinators//
  ///////////////

  // compose transforms in sequence
  // sequences return true if any transform made a change, even if it was
  // overwritten later
  friend Transform operator>>(const Transform& lhs, const Transform& rhs);
  static Transform sequence(std::vector<Transform>& tvec);

  // repeats a transform until it makes no changes (returns false)
  static Transform repeat(const Transform& trans);

  // repeats a transform and stops when the metric stops decreasing
  static Transform repeat_with_metric(
      const Transform& trans, const Metric& eval);

  static Transform repeat_while(const Transform& cond, const Transform& body);

  ///////////////
  // MeasurePass//
  ///////////////

  /**
   * Commute all measurement gates to the end of the circuit.
   * Throws a CircuitInvalidity exception if it is not possible to delay.
   */
  static Transform delay_measures();

  //////////////////
  // Decompositions//
  //////////////////

  /**
   * Decomposes all multi-qubit unitary gates into CX and single-qubit gates.
   *
   * Ignores boxes.
   *
   * Expects: any gates
   * Produces: CX and any single-qubit gates
   */
  static Transform decompose_multi_qubits_CX();

  /**
   * Decomposes all single-qubit unitary gates into TK1 gates. Ignores boxes.
   */
  static Transform decompose_single_qubits_TK1();

  /**
   * Starting with Rz, Ry and multi-qubit gates, replace all singles with TK1.
   */
  static Transform decompose_ZYZ_to_TK1();

  // converts all single-qubit gates into Rz and Rx gates
  // Expects: any gates
  // Produces: Rz, Rx and any multi-qubit gates
  static Transform decompose_ZX();

  // converts all single-qubit gates into Rz and Ry gates
  // Expects: any gates
  // Produces: Rz, Ry and any multi-qubit gates
  static Transform decompose_ZY();

  // converts all single-qubit gates into Rx and Ry gates
  // Expects: any gates
  // Produces: Rx, Ry and any multi-qubit gates
  static Transform decompose_XY();

  // converts all tk1-gates into rzrx gates
  // Expects: tk1-gates and any multi-qubit gates
  // Produces: Rz, Rx and any multi-qubit gates
  static Transform decompose_tk1_to_rzrx();

  // replaces CXs with ECR, SX, Rz
  // Expects: CX and any single-qubit gates
  // Produces: ECR and any single-qubit gates
  static Transform decompose_CX_to_ECR();

  // replaces CXs with ZZMax
  // Expects: CX and any single-qubit gates
  // Produces: ZZMax and any single-qubit gates
  static Transform decompose_CX_to_HQS2();

  // replaces Rz-Rx-Rz triples with Rz-PhasedX pairs
  // Expects: Rz, Rx, and any multi-qubit gates
  // Produces: Rz, PhasedX, and any multi-qubit gates
  static Transform decompose_ZX_to_HQS1();

  // converts all CX gates into Molmer-Sorensen gates by recognising exp(-i XX *
  // angle * pi/2) and converting rest to exp(-i XX * pi/4) Expects: CX, Rx, and
  // other single-qubit gates Produces: Molmer-Sorensen, U2, and other
  // single-qubit gates
  static Transform decompose_MolmerSorensen();

  // finds all ZZ gates which are in the form of CX and Rz gates and converts
  // into ZZ gates Expects: CX, Rz Produces: ZZ, CX, Rz
  static Transform decompose_ZZPhase();

  // identifies any tk1, Rz,Ry,Rx,u3,u2,u1 gate sequences that can be naively
  // decomposed into S,Z,X,V returns clifford sequences in a standard form
  // (Z)(X)(S)(V)(S)
  // Expects: tk1, U3, U2, U1, Rx, Ry, Rz, and any multi-qubit gates
  // Produces: Z, X, S, V, and any remaining non-clifford single-qubit gates,
  // and any multi-qubit gates
  static Transform decompose_cliffords_std();
  static Transform decompose_ZX_to_cliffords();

  // converts all valid CX-Rz-CX strings and CX-Rx-CX strings to 2qb
  // PhaseGadgets Expects: CX, Rz, Rx and any other single-qubit gates Produces:
  // PhaseGadgets, and any other gates
  static Transform decompose_PhaseGadgets();

  // converts all SWAP gates to given replacement circuit (not checked to
  // preserve unitary) Expects: SWAP gates, replacement circuit Produces:
  // Instances of the replacement circuit,
  static Transform decompose_SWAP(const Circuit& replacement_circuit);

  // converts all SWAP gates to 3 CX gates
  // providing an Architecture will prefer an orientation that reduces H
  // redirection cost -X-    -C-X-C-   -X-C-X-
  //  |   =  | | |  =  | | |
  // -X-    -X-C-X-   -C-X-C-
  // Expects: SWAP gates
  // Produces: CX gates
  static Transform decompose_SWAP_to_CX(
      const Architecture& arc = Architecture());

  // converts all BRIDGE (distance 2 CX) gates to 4 CX gates
  // -B-   -C-   -C---C---   ---C---C-
  //  R     |     |   |         |   |
  // -I- = --- = -X-C-X-C- = -C-X-C-X-
  //  D     |       |   |     |   |
  // -G-   -X-   ---X---X-   -X---X---
  //  E
  // Expects: BRIDGE gates
  // Produces: CX gates
  static Transform decompose_BRIDGE_to_CX();

  // converts -C- gates to -H-X-H- depending on direction of Architecture edges
  //           |              |
  //          -X-          -H-C-H-
  // Expects: CX gates
  // Produces CX and H gates
  static Transform decompose_CX_directed(const Architecture& arc);

  // converts arbitrarily controlled Ry gates. It does not use ancillae, so it
  // is not very depth-efficient Expects: CRys and any other gates returns Ry,
  // CX, H, T, Tdg + whatever other gates were there before
  static Transform decomp_controlled_Rys();

  // does not use ancillae
  // Expects: CCX + any other gates
  // returns CX, H, T, Tdg + any previous gates
  static Transform decomp_CCX();

  // does not use ancillae
  // Expects: any CnRys + CnXs + any other gates
  // returns Ry, CX, H, T, Tdg + any previous gates
  static Transform decomp_arbitrary_controlled_gates();

  /**
   * Replaces all boxes by their decomposition using Box::to_circuit
   * Expects: any gateset
   * returns potentially all gates
   */
  static Transform decomp_boxes();

  /**
   * Replaces all CX+Rz sub circuits by PhasePolyBox
   * Expects: only CX + Rz + H (and measure + reset + collapse + barrier)
   * returns: H + PhasePolyBox
   */
  static Transform compose_phase_poly_boxes();

  ///////////////
  // Rebase Pass//
  ///////////////

  // decomposes multiq gates not in the gate set to CXs, then replaces CXs with
  // the replacement (if CX is not allowed) then converts singleq gates no in
  // the gate set to U3 and replaces them using provided function Expects: any
  // gates Produces: gates in multiqs and singleqs
  static Transform rebase_factory(
      const OpTypeSet& multiqs, const Circuit& cx_replacement,
      const OpTypeSet& singleqs,
      const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
          tk1_replacement);

  // Multiqs: CX
  // Singleqs: tk1
  static Transform rebase_tket();

  // Multiqs: CZ
  // Singleqs: PhasedX, Rz
  static Transform rebase_cirq();

  // Multiqs: CZ
  // Singleqs: Rx, Rz
  static Transform rebase_quil();

  // Multiqs: SWAP, CX, CZ
  // Singleqs: H, X, Z, S, T, Rx, Rz
  static Transform rebase_pyzx();

  // Multiqs: SWAP, CRz, CX, CZ
  // Singleqs: H, X, Y, Z, S, T, V, Rx, Ry, Rz
  static Transform rebase_projectq();

  // Multiqs: ZZMax
  // Singleqs: PhasedX, Rz
  static Transform rebase_HQS();

  // Multiqs: XXPhase
  // Singleqs: PhasedX, Rz
  static Transform rebase_UMD();

  // Used for UniversalFrameRandomisation
  // Multiqs: CX
  // Singleqs: Rz, H
  static Transform rebase_UFR();

  // Multiqs: ECR
  // Singleqs: Rz, SX
  static Transform rebase_OQC();

  //////////////////
  // Synthesis Pass//
  //////////////////

  /*Synthesis passes preserve connectivity */

  /**
   * Synthesise a circuit consisting of CX and TK1 gates only.
   */
  static Transform synthesise_tket();

  // converts a circuit into the HQS primitives (Rz, PhasedX, ZZMax) whilst
  // optimising Expects: CX and any single-qubit gates Produces: ZZMax, PhasedX,
  // Rz
  static Transform synthesise_HQS();

  // converts a circuit into the OQC primitives (Rz, SX, ECR gates)
  // Expects: any gates
  // Produces: Rz, SX, ECR
  static Transform synthesise_OQC();

  // converts a circuit into the UMD primitives (Rz, PhasedX, XXPhase) whilst
  // optimising Expects: Any gate set Produces: XXPhase, PhasedX, Rz
  static Transform synthesise_UMD();

  //////////////////////////
  // Full Optimisation Pass//
  //////////////////////////

  /*these Transform passes do not preserve connectivity*/

  // simplifies a circuit using Clifford rules
  // Expects: CX and any single-qubit gates
  // Produces: CX, tk1
  static Transform clifford_simp(bool allow_swaps = true);

  // runs clifford_simp
  // Expects: Any gates
  // Produces: CX, tk1
  static Transform hyper_clifford_squash();

  // kitchen sink optimisation - phase gadget resynthesis, two-qubit Cartan
  // forms, Clifford Expects: Any gates Produces: CX, tk1
  static Transform canonical_hyper_clifford_squash();

  // only peephole optimisation, so no higher structure abstraction.
  // Two qubit Cartan, Clifford, synthesis
  // Expects: Any gates
  // Produces: CX, tk1
  static Transform peephole_optimise_2q();

  /**
   * Peephole optimisation including resynthesis of three-qubit gate sequences.
   *
   * @param allow_swaps whether to allow introduction of implicit wire swaps
   *
   * Produces: CX, tk1.
   */
  static Transform full_peephole_optimise(bool allow_swaps = true);

  //////////////////////
  // Phase Optimisation//
  //////////////////////

  // extends PhaseGadgets by a qubit on identifying a pair of CXs around it
  // Expects: CX, PhaseGadget, and any single-qubit gates
  // Produces: CX, PhaseGadget, and any single-qubit gates
  static Transform smash_CX_PhaseGadgets();

  // tries to match up ports on adjacent PhaseGadgets to enable maximal
  // annihilation after synthesis (ignoring any intervening gates) Expects: Any
  // gates Produces: The same gate set
  static Transform align_PhaseGadgets();

  /////////////////////////
  // Clifford Optimisation//
  /////////////////////////

  // These apply some of the Clifford rules in the paper "Optimising Clifford
  // Circuits with Quantomatic" All of these expect and produce CX, Z, X, S, V,
  // and other single qubit gates

  static Transform singleq_clifford_sweep();

  // mutli-qubit patterns that decrease the CX count
  // inserting swaps can sometimes cause errors elsewhere (e.g. routing), so
  // they can be turned off
  static Transform multiq_clifford_replacement(bool allow_swaps = false);

  static Transform clifford_reduction(bool allow_swaps = false);

  // copies Z through the target of a CX and X through the control
  static Transform copy_pi_through_CX();

  /////////////////////////////
  // Pauli Gadget Optimisation//
  /////////////////////////////

  /**
   * Depth-saving resynthesis of phase gadgets with alignment.
   *
   * Produces CX and TK1 gates.
   */
  static Transform optimise_via_PhaseGadget(
      CXConfigType cx_config = CXConfigType::Snake);

  static Transform pairwise_pauli_gadgets(
      CXConfigType cx_config = CXConfigType::Snake);

  // always returns true, as it leaves Circuit data structure
  static Transform synthesise_pauli_graph(
      PauliSynthStrat strat = PauliSynthStrat::Sets,
      CXConfigType cx_config = CXConfigType::Snake);

  // Assumes incoming circuit is composed of `CircBox`es with
  // `PauliExpBox`es inside
  static Transform special_UCC_synthesis(
      PauliSynthStrat strat = PauliSynthStrat::Sets,
      CXConfigType cx_config = CXConfigType::Snake);

  //////////////////////
  // Basic Optimisation//
  //////////////////////

  // removes gate-inverse pairs, merges rotations, removes identity rotations,
  // and removes redundant gates before measure Expects: Any gates Produces: The
  // same gate set
  static Transform remove_redundancies();

  /**
   * general u_squash by converting any chains of p, q gates (p, q in
   * {Rx,Ry,Rz}) to triples -p-q-p- or -q-p-q-
   *
   * if strict = false (default), then the chains will be squashed to either
   * -p-q-p- or -q-p-q-, and the third rotation will be commuted through the
   * next gate whenever possible if strict = true, then all chains will be
   * squashed to -p-q-p-
   *
   * Expects: p, q, and any multi-qubit gates
   * Produces: p, q, and any multi-qubit gates
   */
  static Transform squash_1qb_to_pqp(
      const OpType& q, const OpType& p, bool strict = false);

  // 1qb squashing into -Rz-Rx-Rz- or -Rx-Rz-Rx- form
  // Expects: Rx, Rz, and any multi-qubit gates
  // Produces: Rx, Rz, and any multi-qubit gates
  static Transform reduce_XZ_chains();

  /**
   * Squash all single-qubit gates to TK1.
   */
  static Transform squash_1qb_to_tk1();

  // identifies single-qubit chains and squashes them in the target gate set
  // Expects: any gates
  // Produces: singleqs and any multi-qubit gates
  static Transform squash_factory(
      const OpTypeSet& singleqs,
      const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
          tk1_replacement);

  // moves single qubit operations past multiqubit operations they commute with,
  // towards front of circuit (hardcoded)
  // Expects: Any gates
  // Produces: Any gates
  static Transform commute_through_multis();

  // commutes Rz gates through ZZMax, and combines adjacent ZZMax gates
  // Expects: ZZMax, Rz, Rx
  // Produces: ZZMax, Rz, Rx
  static Transform commute_and_combine_HQS2();

  // commutes single qubit gates through SWAP gates, leaving them on the
  // PhysicalQubit with best fidelity for the given OP Expects: any single qubit
  // gates, SWAP gates Produces: any single qubit gates, SWAP gates
  static Transform commute_SQ_gates_through_SWAPS(
      const avg_node_errors_t& node_errors);
  static Transform commute_SQ_gates_through_SWAPS(
      const op_node_errors_t& node_errors);

  // squashes each sequence of two-qubit instructions into their canonical 3-CX
  // form (Cartan/KAK decomposition) The optional CX gate fidelity cx_fidelity
  // is used to produce circuit decompositions that maximise expected circuit
  // fidelity, assuming perfect local operations (i.e. with fidelity 1.).
  // Expects: CX and any single qubit gates
  // Produces: CX and single tk1 qubit gates
  static Transform two_qubit_squash(double cx_fidelity = 1.);

  /**
   * Squash sequences of 3-qubit instructions into their canonical 20-CX form.
   *
   * The circuit should comprise only CX and single-qubit gates. The transform
   * may perform a combination of 2-qubit (KAK) and 3-qubit decompositions of
   * subcircuits, but only does so if this reduces the CX count.
   *
   * @return Transform implementing the squash
   */
  static Transform three_qubit_squash();

  ////////////////////////
  // Contextual Reduction//
  ////////////////////////

  /**
   * Remove all operations that have no @ref OpType::Output or
   * @ref OpType::ClOutput in their causal future.
   */
  static Transform remove_discarded_ops();

  /** Allow insertion of classical operations when simplifying? */
  enum class AllowClassical { Yes, No };

  /** Automatically create all qubits in zero state when simplifying? */
  enum class CreateAllQubits { Yes, No };

  /**
   * Simplify the circuit where it acts on known basis states.
   *
   * Whenever a gate transforms a known basis state to another known basis
   * state, remove it, inserting X gates where necessary to achieve the same
   * state. Beginning with Create and Reset vertices, move forward through the
   * circuit as far as we can in this way.
   *
   * Global phase is ignored.
   *
   * @param allow_classical if allowed, insert SetBits operations in place of
   *  Measure operations that act on known basis states
   * @param create_all_qubits if enabled, annotate all qubits as initialized
   *  to zero as part of the transform, before applying simplification
   * @param xcirc 1-qubit circuit implementing an X gate (if null, an X gate
   *  is used)
   */
  static Transform simplify_initial(
      AllowClassical allow_classical = AllowClassical::Yes,
      CreateAllQubits create_all_qubits = CreateAllQubits::No,
      std::shared_ptr<const Circuit> = 0);

  /**
   * Commute classical maps through measurements.
   *
   * A "classical map" is a pure quantum operation that acts as a permutation
   * of the computational basis states composed with a diagonal operator.
   *
   * A "measured classical map" is a classical map whose succeeding vertices
   * in the DAG are all @ref OpType::Measure.
   *
   * This transform replaces all measured classical maps that are followed by
   * Measure operations whose quantum output is discarded with classical
   * operations following the Measure.
   *
   * The process is repeated until no such replacements are possible.
   *
   * Global phase is not preserved.
   */
  static Transform simplify_measured();

  ///////////////////////////////////////////////
  // Minor components for other Transform passes//
  ///////////////////////////////////////////////

  // converts a tk1 gate to a PhasedXRz gate
  static Circuit tk1_to_PhasedXRz(
      const Expr& alpha, const Expr& beta, const Expr& gamma);

  static Circuit tk1_to_rzrx(
      const Expr& alpha, const Expr& beta, const Expr& gamma);

  static Circuit tk1_to_rzh(
      const Expr& alpha, const Expr& beta, const Expr& gamma);

  static Circuit tk1_to_tk1(
      const Expr& alpha, const Expr& beta, const Expr& gamma);

  static Circuit tk1_to_rzsx(
      const Expr& alpha, const Expr& beta, const Expr& gamma);

  // helper class subcircuits representing 2qb interactions
  struct Interaction {
    Interaction(const Qubit& _q0, const Qubit& _q1) : q0(_q0), q1(_q1) {}
    Qubit q0;  // Qubit numbers
    Qubit q1;
    Edge e0;  // In edges starting interaction
    Edge e1;
    unsigned count;      // Number of two qubit gates in interaction
    VertexSet vertices;  // Vertices in interaction subcircuit
  };

  static Circuit cnx_normal_decomp(unsigned n);
  static Circuit decomposed_CnRy(const Op_ptr op, unsigned arity);
  static Circuit incrementer_borrow_n_qubits(unsigned n);
  static Circuit incrementer_borrow_1_qubit(unsigned n);
};

inline const Transform Transform::id = Transform([](const Circuit&) {
  return false;
});  // returns `false` as it does not change the Circuit in any way

}  // namespace tket

#endif
