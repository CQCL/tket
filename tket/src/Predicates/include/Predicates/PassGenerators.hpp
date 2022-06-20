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

#include "ArchAwareSynth/SteinerForest.hpp"
#include "CompilerPass.hpp"
#include "Mapping/LexiRoute.hpp"
#include "Mapping/RoutingMethod.hpp"
#include "Transformations/ContextualReduction.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/PauliOptimisation.hpp"

namespace tket {

/* a wrapper method for the rebase_factory in Transforms */
PassPtr gen_rebase_pass(
    const OpTypeSet& allowed_gates, const Circuit& cx_replacement,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement);

/* a wrapper method for the squash_factory in Transforms */
PassPtr gen_squash_pass(
    const OpTypeSet& singleqs,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement);
PassPtr gen_euler_pass(const OpType& q, const OpType& p, bool strict = false);
PassPtr gen_clifford_simp_pass(bool allow_swaps = true);

/**
 * Pass to rename some or all qubits according to the given map.
 *
 * Qubits in the map that do not occur in the circuit are ignored.
 */
PassPtr gen_rename_qubits_pass(const std::map<Qubit, Qubit>& qm);

PassPtr gen_placement_pass(const PlacementPtr& placement_ptr);

PassPtr gen_naive_placement_pass(const Architecture& arc);
/* This higher order function generates a Routing pass using the
std::vector<RoutingMethodPtr> object */
PassPtr gen_full_mapping_pass(
    const Architecture& arc, const PlacementPtr& placement_ptr,
    const std::vector<RoutingMethodPtr>& config);
PassPtr gen_default_mapping_pass(
    const Architecture& arc, bool delay_measures = true);
PassPtr gen_cx_mapping_pass(
    const Architecture& arc, const PlacementPtr& placement_ptr,
    const std::vector<RoutingMethodPtr>& config, bool directed_cx,
    bool delay_measures);
PassPtr gen_routing_pass(
    const Architecture& arc, const std::vector<RoutingMethodPtr>& config);
PassPtr gen_directed_cx_routing_pass(
    const Architecture& arc, const std::vector<RoutingMethodPtr>& config);

/**
 * execute architecture aware synthesis on a given architecture for an allready
 * place circuit, only for circuit which contains Cx+Rz+H gates
 * this pass is not able to handle implicit wire swaps
 * @param arc architecture to route on
 * @param lookahead parameter for the recursion depth in the algorithm, the
 * value should be > 0
 * @param cnotsynthtype parameter for the type of cnot synth
 * @return passpointer to perform architecture aware synthesis
 */
PassPtr aas_routing_pass(
    const Architecture& arc, const unsigned lookahead = 1,
    const aas::CNotSynthType cnotsynthtype = aas::CNotSynthType::Rec);

/**
 * execute architecture aware synthesis on a given architecture for any circuit.
 * all unplaced qubits will be placed in this pass
 * @param arc architecture to route on
 * @param lookahead parameter for the recursion depth in the algorithm, the
 * value should be > 0
 * @param cnotsynthtype parameter for the type of cnot synth
 * @return passpointer to perform architecture aware synthesis
 */
PassPtr gen_full_mapping_pass_phase_poly(
    const Architecture& arc, const unsigned lookahead = 1,
    const aas::CNotSynthType cnotsynthtype = aas::CNotSynthType::Rec);

/**
 * pass to place all not yet placed qubits of the circuit to the given
 * architecture for the architecture aware synthesis.
 * @param arc achitecture to place the circuit on
 * @return passpointer to perfomr the mapping
 */
PassPtr gen_placement_pass_phase_poly(const Architecture& arc);

PassPtr gen_decompose_routing_gates_to_cxs_pass(
    const Architecture& arc = Architecture(), bool directed = false);

/* generates a decomposition pass that converts all SWAPs into a chosen
 * replacement circuit */
PassPtr gen_user_defined_swap_decomp_pass(const Circuit& replacement_circ);

/**
 * Pass to decompose sequences of two-qubit operations into optimal gate
 * sequence
 *
 * Generate a KAK decomposition pass that can decompose any two-qubit subcircuit
 * into <= 3 CX + local operations.
 * If an estimate for the CX gate fidelity cx_fidelity is given, the
 * decomposition will trade off the accuracy of the circuit approximation
 * against the loss in circuit fidelity induced by additional noisy CX gates.
 * This will result in an output circuit that is not necessarily equivalent to
 * the input, but that maximises expected circuit fidelity
 */
PassPtr KAKDecomposition(double cx_fidelity = 1.);

/**
 * @brief Decomposes each TK2 gate into two-qubit gates.
 *
 * We currently support CX, ZZMax and ZZPhase.
 *
 * If one or more gate fidelities are provided, the two-qubit gate
 * type achieving the highest fidelity will be chosen for the
 * decomposition, as measured using squared trace fidelity.
 * If no fidelities are provided, the TK2 gates will be decomposed
 * exactly using CX gates.
 *
 * All TK2(α, β, γ) gates must be normalised to the Weyl chamber, i.e.
 * 0.5 ≥ 𝛼 ≥ 𝛽 ≥ |𝛾|.
 *
 * Gate fidelities are passed as keyword arguments to perform noise-aware
 * decompositions. We currently support `CX_fidelity`, `ZZMax_fidelity` and
 * `ZZPhase_fidelity`. If provided, the `CX` and `ZZMax` fidelities must be
 * given by a single floating point fidelity. The `ZZPhase` fidelity is given as
 * a lambda float -> float, mapping a ZZPhase angle parameter to its fidelity.
 * These parameters will be used to return the optimal decomposition of each TK2
 * gate, taking noise into consideration.
 *
 * If the TK2 angles are symbolic values, the decomposition will be exact
 * (i.e. not noise-aware). It is not possible in general to obtain optimal
 * decompositions for arbitrary symbolic parameters, so consider substituting
 * for concrete values if possible.
 *
 * @param fid The two-qubit gate fidelities (optional).
 * @return PassPtr
 */
PassPtr DecomposeTK2(const Transforms::TwoQbFidelities& fid);
PassPtr DecomposeTK2();

/**
 * Resynthesize and squash three-qubit interactions.
 *
 * Steps through the circuit accumulating sequences of 2- and 3-qubit
 * interactions, where possible squashing them into subcircuits having lower CX
 * count, then applies Clifford simplification.
 *
 * @param allow_swaps whether to allow introduction of implicit swaps
 */
PassPtr ThreeQubitSquash(bool allow_swaps = true);

/**
 * Performs peephole optimisation including resynthesis of 2- and 3-qubit gate
 * sequences, and converts to a circuit containing only CX and TK1 gates.
 *
 * @param allow_swaps whether to allow introduction of implicit swaps
 */
PassPtr FullPeepholeOptimise(bool allow_swaps = true);

/* generates an optimisation pass that converts a circuit into phase
gadgets and optimises them using techniques from
https://arxiv.org/abs/1906.01734 */
PassPtr gen_optimise_phase_gadgets(
    CXConfigType cx_config = CXConfigType::Snake);

/* generates an optimisation pass that converts a circuit into Pauli
gadgets and optimises them using techniques from
https://arxiv.org/abs/1906.01734 */
PassPtr gen_pairwise_pauli_gadgets(
    CXConfigType cx_config = CXConfigType::Snake);

/* generates an optimisation pass that converts a circuit into a graph
of Pauli gadgets and optimises them using strategies from <paper to come> */
PassPtr gen_synthesise_pauli_graph(
    Transforms::PauliSynthStrat strat = Transforms::PauliSynthStrat::Sets,
    CXConfigType cx_config = CXConfigType::Snake);

/* generates an optimisation pass that converts a circuit built using
term sequencing techniques from <paper to come> into a graph of Pauli
gadgets and optimises them. */
PassPtr gen_special_UCC_synthesis(
    Transforms::PauliSynthStrat strat = Transforms::PauliSynthStrat::Sets,
    CXConfigType cx_config = CXConfigType::Snake);

/**
 * Generate a pass to simplify the circuit where it acts on known basis states.
 *
 * @param allow_classical allow replacement of measures by pure classical set-
 *      bit operations when the measure acts on a known state
 * @param create_all_qubits if enabled, annotate all qubits as initialized
 *      to zero as part of the transform, before applying simplification
 * @param xcirc 1-qubit circuit implementing an X gate (if null, an X gate is
 *      used)
 */
PassPtr gen_simplify_initial(
    Transforms::AllowClassical allow_classical =
        Transforms::AllowClassical::Yes,
    Transforms::CreateAllQubits create_all_qubits =
        Transforms::CreateAllQubits::No,
    std::shared_ptr<const Circuit> xcirc = 0);

/**
 * Generate a pass to perform simplifications dependent on qubit state.
 *
 * @param allow_classical allow insertion of classical operations
 * @param xcirc 1-qubit circuit implementing an X gate (if null, an X gate is
 *      used)
 */
PassPtr gen_contextual_pass(
    Transforms::AllowClassical allow_classical =
        Transforms::AllowClassical::Yes,
    std::shared_ptr<const Circuit> xcirc = 0);

/**
 * Builds a sequence of PauliSimp (gen_synthesise_pauli_graph) and
 * FullPeepholeOptimise
 */
PassPtr PauliSquash(Transforms::PauliSynthStrat strat, CXConfigType cx_config);

/**
 * Turns all PhasedX and NPhasedX gates into global gates
 *
 * Replaces any PhasedX gates with global NPhasedX gates.
 * By default, this transform will squash all single-qubit gates to PhasedX and
 * Rz gates before proceeding further. Existing non-global NPhasedX will not
 * be preserved. This is the recommended setting for best performance.
 *
 * If squashing is disabled, each non-global PhasedX gate will be replaced with
 * two global NPhasedX, but any other gates will be left untouched.
 *
 * @param squash Whether to squash the circuit in pre-processing
 *      (default: true).
 *
 * If squash=true (default), the `GlobalisePhasedX().apply` method will always
 * return true. For squash=false, `apply()` will return true if the circuit was
 * changed and false otherwise.
 *
 * It is not recommended to use this pass with symbolic expressions, as
 * in certain cases a blow-up in symbolic expression sizes may occur.
 */
PassPtr GlobalisePhasedX(bool squash = true);

}  // namespace tket
