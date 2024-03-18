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

#pragma once

#include "PGBox.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/PauliGraphRefactor/PauliGraph.hpp"

namespace tket {

// Copied from Converters/Converters.hpp to allow synthesis and testing of reset
// and boxes
ChoiMixTableau circuit_to_cm_tableau(const Circuit& circ);
std::pair<Circuit, qubit_map_t> cm_tableau_to_exact_circuit(
    const ChoiMixTableau& tab, CXConfigType cx_config = CXConfigType::Snake);
std::pair<Circuit, qubit_map_t> cm_tableau_to_unitary_extension_circuit(
    const ChoiMixTableau& tab, const std::vector<Qubit>& init_names = {},
    const std::vector<Qubit>& post_names = {},
    CXConfigType cx_config = CXConfigType::Snake);

/**
 * Converts the Circuit to a PauliGraph representing the same circuit by
 * iterating through the circuit, accumulating a UnitaryRevTableau of Clifford
 * operations and yielding PGOps for any other Op. Qubit initialisations are
 * recorded as part of the PGInputTableau, and discards are incorporated into
 * the accumulated tableau to give the PGOutputTableau.
 */
pg::PauliGraph circuit_to_pauli_graph3(
    const Circuit& circ, bool collect_cliffords = true);

/**
 * THIS IS A LEGACY METHOD DESIGNED TO REPLICATE THE SYNTHESIS METHODS FOR THE
 * EXISTING PAULIGRAPH
 *
 * Synthesises a circuit equivalent to the PauliGraph by adding each
 * pauli gadget to the circuit as a PauliExpBox individually
 * in the order given by TopSortIterator.
 * The tableau is then synthesised at the end.
 */
Circuit pauli_graph3_to_pauli_exp_box_circuit_individually(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * THIS IS A LEGACY METHOD DESIGNED TO REPLICATE THE SYNTHESIS METHODS FOR THE
 * EXISTING PAULIGRAPH
 *
 * Synthesises a circuit equivalent to the PauliGraph by inserting pairs of
 * pauli gadgets as PauliExpPairBoxes into the circuit
 * The tableau is then synthesised at the end.
 */
Circuit pauli_graph3_to_pauli_exp_box_circuit_pairwise(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * THIS IS A LEGACY METHOD DESIGNED TO REPLICATE THE SYNTHESIS METHODS FOR THE
 * EXISTING PAULIGRAPH
 *
 * Synthesises a circuit equivalent to the PauliGraph by building
 * sets of mutually commuting pauli gadgets and
 * inserting them into the circuit as PauliExpCommutingSetBoxes
 * The tableau is then synthesised at the end.
 */

Circuit pauli_graph3_to_pauli_exp_box_circuit_sets(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * Synthesises a circuit equivalent to the PauliGraph by synthesising each
 * vertex individually in the order given by TopSortIterator. The tableaux are
 * synthesised at each end using the default synthesis method for
 * ChoiMixTableau.
 */
Circuit pauli_graph3_to_circuit_individual(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

Circuit pauli_graph3_to_circuit_sets(
    const pg::PauliGraph& pg, CXConfigType cx_config = CXConfigType::Snake);

}  // namespace tket
