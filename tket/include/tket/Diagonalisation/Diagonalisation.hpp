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

#include "DiagUtils.hpp"
#include "tket/Circuit/Boxes.hpp"
#include "tket/Circuit/Circuit.hpp"

namespace tket {

enum class DiagonalisationMethod {
  Greedy,
  JGM,
  vdBT_PE,
  vdBT_CX,
  vdBT_greedy1,
  vdBT_greedy2,
  CSW_CZ,
  CSW_CX
};

NLOHMANN_JSON_SERIALIZE_ENUM(
    DiagonalisationMethod,
    {
        {DiagonalisationMethod::Greedy, "Greedy"},
        {DiagonalisationMethod::JGM, "JGM"},
        {DiagonalisationMethod::vdBT_PE, "vdBT_PE"},
        {DiagonalisationMethod::vdBT_CX, "vdBT_CX"},
        {DiagonalisationMethod::vdBT_greedy1, "vdBT_greedy1"},
        {DiagonalisationMethod::vdBT_greedy2, "vdBT_greedy2"},
        {DiagonalisationMethod::CSW_CZ, "CSW_CZ"},
        {DiagonalisationMethod::CSW_CX, "CSW_CX"},
    });

/**
 * Diagonalise a mutually commuting set of Pauli strings. Modifies the
 * list of Pauli strings in place, and returns the Clifford circuit
 * required to generate the initial set.
 */
Circuit mutual_diagonalise(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType cx_config,
    DiagonalisationMethod diag_meth = DiagonalisationMethod::Greedy);

/**
 * Check whether there are any qubits which only requires a single
 * qubit Clifford to make all Paulis I or Z
 */
void check_easy_diagonalise(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> &qubits,
    Circuit &circ);

/**
 * Given two qubits, attempt to find a basis in which a single CX will
 * make the Paulis on one of qubits fully diagonal
 */
std::optional<std::pair<Pauli, Pauli>> check_pair_compatibility(
    const Qubit &qb1, const Qubit &qb2,
    const std::list<SpSymPauliTensor> &gadgets);

/**
 * Diagonalise a qubit greedily by finding the Pauli Gadget with
 * the lowest residual support over the non-diagonal qubits, and apply
 * single qubit Cliffords and CXs to make it a `ZIII...I` string
 */
void greedy_diagonalise(
    const std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> &qubits,
    Conjugations &conjugations, Circuit &circ, CXConfigType cx_config);

/**
 * Implements mutual_diagonalise for the Greedy method.
 */
Circuit mutual_diagonalise_greedy(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType cx_config);
/**
 * Applies Clifford conjugations to a SpSymPauliTensor
 */
void apply_conjugations(
    SpSymPauliTensor &qps, const Conjugations &conjugations);

/**
 *  Given two qubits on which to conjugate a CX gate, try to conjugate with a
 * XXPhase3 instead. If successful, undoes conjugations that must be undone and
 * replaces it with XXPhase3 conjugation. Returns true if successful and false
 * otherwise.
 */
bool conjugate_with_xxphase3(
    const Qubit &qb_a, const Qubit &qb_b, Conjugations &conjugations,
    Circuit &cliff_circ);

ChoiMixTableau tab_from_gadgets(const std::list<SpSymPauliTensor> &gadgets);

/**
 * Implements mutual_diagonalise for the method in Appendix A of Jena, Genin,
 * Mosca, "Pauli Partitioning with Respect to Gate Sets",
 * https://arxiv.org/pdf/1907.07859.pdf
 */
Circuit mutual_diagonalise_JGM(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType cx_config);

/**
 * Implements mutual_diagonalise for the methods in Section 4 of van den Berg &
 * Temme, "Circuit optimisation of Hamiltonian simulation by simultaneous
 * diagonalization of Pauli clusters",
 * https://quantum-journal.org/papers/q-2020-09-12-322/
 */
Circuit mutual_diagonalise_vdBT_PE(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType cx_config);
Circuit mutual_diagonalise_vdBT_CX(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType cx_config);
Circuit mutual_diagonalise_vdBT_greedy(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType cx_config, bool subsort_by_singles);

/**
 * Implements mutual_diagonalise for the methods in Section 3 of Crawford, van
 * Straaten, Wang, Parks, Campbell, Brierley, "Efficient quantum measurement of
 * Pauli operators in the presence of finite sampling error",
 * https://arxiv.org/pdf/1908.06942.pdf
 */
Circuit mutual_diagonalise_CSW_CZ(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType cx_config);
Circuit mutual_diagonalise_CSW_CX(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType cx_config);

}  // namespace tket
