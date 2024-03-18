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

#include "tket/Circuit/Circuit.hpp"
#include "tket/Clifford/ChoiMixTableau.hpp"
#include "tket/Clifford/UnitaryTableau.hpp"
#include "tket/PauliGraph/PauliGraph.hpp"
#include "tket/ZX/ZXDiagram.hpp"

namespace tket {

/**
 * Construct the tableau for a given circuit.
 * Will throw an exception if it contains non-Clifford gates.
 */
UnitaryTableau circuit_to_unitary_tableau(const Circuit &circ);
UnitaryRevTableau circuit_to_unitary_rev_tableau(const Circuit &circ);

/**
 * Constructs a circuit producing the same effect as the tableau.
 * Uses the method from Aaronson-Gottesman: Improved Simulation of
 * Stabilizer Circuits, Theorem 8.
 * CAUTION: GATE COUNT IS ATROCIOUS IN PRACTICE
 */
Circuit unitary_tableau_to_circuit(const UnitaryTableau &tab);
Circuit unitary_rev_tableau_to_circuit(const UnitaryRevTableau &tab);

/**
 * Construct a ChoiMixTableau for a given circuit.
 * Will incorporate qubit initialisations and discarding into the circuit.
 * Will throw an exception if it contains non-Clifford gates.
 */
ChoiMixTableau circuit_to_cm_tableau(const Circuit &circ);

/**
 * Constructs a circuit producing the same effect as a ChoiMixTableau.
 * Uses a naive synthesis method until we develop a good heuristic.
 * Since Circuit does not support distinct qubit addresses for inputs and
 * outputs, also returns a map from the output qubit IDs in the tableau to their
 * corresponding outputs in the circuit
 */
std::pair<Circuit, unit_map_t> cm_tableau_to_circuit(
    const ChoiMixTableau &circ);

/**
 * Convert a tableau for a unitary to the equivalent ChoiMixTableau. This
 * enables composition with non-unitary stabiliser operations.
 */
ChoiMixTableau unitary_tableau_to_cm_tableau(const UnitaryTableau &tab);
ChoiMixTableau unitary_rev_tableau_to_cm_tableau(const UnitaryRevTableau &tab);

/**
 * Convert a ChoiMixTableau representing a unitary to a more specialised
 * tableau. This enables simpler and faster calculations of Pauli conjugations
 * (i.e. pushing a given Pauli string from one side of the tableau to the
 * other). If this is provided with a non-unitary process (requiring the names
 * of input and output qubits to be identical, and have 2n rows for n qubits),
 * an exception is thrown.
 */
UnitaryTableau cm_tableau_to_unitary_tableau(const ChoiMixTableau &tab);
UnitaryRevTableau cm_tableau_to_unitary_rev_tableau(const ChoiMixTableau &tab);

PauliGraph circuit_to_pauli_graph(const Circuit &circ);

/**
 * Synthesises a circuit equivalent to the PauliGraph by adding each
 * pauli gadget to the circuit as a PauliExpBox individually
 * in the order given by TopSortIterator.
 * The tableau is then synthesised at the end.
 */
Circuit pauli_graph_to_pauli_exp_box_circuit_individually(
    const PauliGraph &pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * Synthesises a circuit equivalent to the PauliGraph by inserting pairs of
 * pauli gadgets as PauliExpPairBoxes into the circuit
 * The tableau is then synthesised at the end.
 */
Circuit pauli_graph_to_pauli_exp_box_circuit_pairwise(
    const PauliGraph &pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * Synthesises a circuit equivalent to the PauliGraph by building
 * sets of mutually commuting pauli gadgets and
 * inserting them into the circuit as PauliExpCommutingSetBoxes
 * The tableau is then synthesised at the end.
 */

Circuit pauli_graph_to_pauli_exp_box_circuit_sets(
    const PauliGraph &pg, CXConfigType cx_config = CXConfigType::Snake);

/**
 * Construct a zx diagram from a given circuit.
 * Return the zx diagram and a map between the zx boundary vertices and the
 * circuit boundary vertices.
 */
std::pair<zx::ZXDiagram, boost::bimap<zx::ZXVert, Vertex>> circuit_to_zx(
    const Circuit &circuit);

/**
 * Takes a unitary ZX diagram in MBQC form with the promise that a gflow exists.
 * Produces an equivalent circuit using the gate extraction method from
 * Backens et al., "There and Back Again: A Circuit Extraction Tale".
 */
Circuit zx_to_circuit(const zx::ZXDiagram &diag);

}  // namespace tket
