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

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Utils/Constants.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/Exceptions.hpp"
#include "Utils/Json.hpp"
#include "Utils/UnitID.hpp"
namespace tket {

/** Pauli not supported */
class UnknownPauli : public std::logic_error {
 public:
  UnknownPauli()
      : std::logic_error("Unknown Pauli. This code should be unreachable!") {}
};
/** OpType not supported */
class UnknownOpType : public std::logic_error {
 public:
  UnknownOpType()
      : std::logic_error(
            "Unknown OpType received when applying conjugations.") {}
};

class UnknownCXConfigType : public std::logic_error {
 public:
  UnknownCXConfigType()
      : std::logic_error(
            "Unknown CXConfigType received when decomposing gadget.") {}
};

class NonDefaultQubit : public std::logic_error {
 public:
  NonDefaultQubit()
      : std::logic_error("Only default register Qubits are supported.") {}
};

typedef Eigen::SparseMatrix<Complex, Eigen::ColMajor> CmplxSpMat;

/** Symbols for the Pauli operators (and identity) */
enum Pauli { I, X, Y, Z };

NLOHMANN_JSON_SERIALIZE_ENUM(
    Pauli, {
               {Pauli::I, "I"},
               {Pauli::X, "X"},
               {Pauli::Y, "Y"},
               {Pauli::Z, "Z"},
           });

/**
 * Whenever a decomposition choice of Pauli gadgets is presented,
 * users may use either Snake (a.k.a. cascade, ladder), Tree (i.e. CX
 * balanced tree) or Star (i.e. CXs target a common qubit).
 */
enum class CXConfigType { Snake, Tree, Star, MultiQGate };

NLOHMANN_JSON_SERIALIZE_ENUM(
    CXConfigType, {{CXConfigType::Snake, "Snake"},
                   {CXConfigType::Tree, "Tree"},
                   {CXConfigType::Star, "Star"},
                   {CXConfigType::MultiQGate, "MultiQGate"}});

typedef std::map<Qubit, Pauli> QubitPauliMap;

/**
 * A string of Pauli letters from the alphabet {I, X, Y, Z},
 * implemented as a sparse list, indexed by qubit
 */
class QubitPauliString {
 public:
  QubitPauliMap map;

  /**
   * Construct an identity map
   */
  QubitPauliString() : map() {}

  /**
   * Construct a single Pauli term
   */
  QubitPauliString(const Qubit &qubit, Pauli p) : map({{qubit, p}}) {}

  /**
   * Construct a string of many Pauli terms. Shortcut to make a
   * string for a default qubit register, without explicitly
   * constructing the QubitPauliMap
   *
   * @param _paulis list of Pauli letters
   */
  explicit QubitPauliString(const std::list<Pauli> &_paulis);

  /**
   * Construct several terms from lists
   * (useful for python hackery, and therefore not extended to QubitPauliTensor,
   * which is not exposed to python)
   *
   * @param qubits list of Qubits with corresponding Paulis
   * @param paulis Pauli letters
   */
  QubitPauliString(
      const std::list<Qubit> &qubits, const std::list<Pauli> &paulis);

  /**
   * Construct a string of many Pauli terms
   *
   * @param _map sparse representation of the full tensor
   */
  explicit QubitPauliString(const QubitPauliMap &_map) : map(_map) {}

  /**
   * Determine whether two strings are equivalent (ignoring I terms)
   *
   * @param other second string
   *
   * @return true iff *this and other represent the same numerical tensor
   */
  bool operator==(const QubitPauliString &other) const;

  /**
   * Determine whether two strings are not equivalent (ignoring I terms)
   *
   * @param other second string
   *
   * @return true iff *this and other represent different numerical tensors
   */
  bool operator!=(const QubitPauliString &other) const;

  /**
   * Determine the lexicographic ordering of two strings.
   * Ignores I terms.
   * Orders individual terms by the ordering of Qubits.
   *
   * @param other other string
   *
   * @return true iff *this orders strictly before other
   */
  bool operator<(const QubitPauliString &other) const;

  /**
   * Compares the lexicographical ordering of two strings.
   * Ignores I terms.
   * Ordered individual terms by the ordering of Qubits.
   *
   * @param other other string
   *
   * @returns -1 if *this orders strictly before other
   * @returns 0 if *this and other represent the same numerical tensor
   * @returns 1 if *this orders strictly after other
   */
  int compare(const QubitPauliString &other) const;

  /**
   * Removes I terms to compress the sparse representation.
   */
  void compress();

  /**
   * Determines whether or not two strings commute.
   *
   * @param other String to compare this to
   *
   * @return true if \p other commutes with this
   * @return false if \p other anti-commutes with this
   */
  bool commutes_with(const QubitPauliString &other) const;

  /**
   * Finds qubits with the same (non-trivial) Pauli term.
   *
   * @param other String to compare this to
   *
   * @return All qubits q where this->map[q] == other.map[q] != I
   */
  std::set<Qubit> common_qubits(const QubitPauliString &other) const;

  /**
   * Finds qubits that only occur in this string.
   *
   * @param other String to compare this to
   *
   * @return All qubits q where this->map[q] != I and other.map[q] == I
   */
  std::set<Qubit> own_qubits(const QubitPauliString &other) const;

  /**
   * Finds qubits with different (non-trivial) Pauli terms.
   *
   * @param other String to compare this to
   *
   * @return All qubits q where I != this->map[q] != other.map[q] != I
   */
  std::set<Qubit> conflicting_qubits(const QubitPauliString &other) const;

  /**
   * Readable string for sparse operator
   */
  std::string to_str() const;

  /**
   * Gets Pauli for a given Qubit
   *
   * @param q Qubit to lookup
   *
   * @return this->map[q] if defined, Pauli::I otherwise
   */
  Pauli get(const Qubit &q) const;

  /**
   * Updates this->map[q] to p
   *
   * @param q Qubit to update
   * @param p Pauli to set this->map[q] = p
   */
  void set(const Qubit &q, Pauli p);

  /**
   * Hash method, required for top-level python. Does not
   * distinguish between equivalent Strings, i.e. ignores I terms.
   */
  friend std::size_t hash_value(const QubitPauliString &qps);

  /**
   * Calculate a sparse matrix corresponding to this string.
   *
   * @param n_qubits Number of qubits the string covers
   *
   * @return A complex sparse matrix corresponding to the n_qubit operator
   */

  CmplxSpMat to_sparse_matrix() const;
  CmplxSpMat to_sparse_matrix(const unsigned n_qubits) const;
  CmplxSpMat to_sparse_matrix(const qubit_vector_t &qubits) const;

  Eigen::VectorXcd dot_state(const Eigen::VectorXcd &state) const;
  Eigen::VectorXcd dot_state(
      const Eigen::VectorXcd &state, const qubit_vector_t &qubits) const;

  Complex state_expectation(const Eigen::VectorXcd &state) const;
  Complex state_expectation(
      const Eigen::VectorXcd &state, const qubit_vector_t &qubits) const;
};

JSON_DECL(QubitPauliString);

// a sum of QubitPauliString with complex coefficients
typedef std::vector<std::pair<QubitPauliString, Complex>> OperatorSum;

/**
 * A simple struct for Pauli strings with +/- phase,
 * used to represent Pauli strings in a stabiliser subgroup
 */
struct PauliStabiliser {
  std::vector<Pauli> string;

  /** true -> 1
   *  false -> -1
   */
  bool coeff;

  PauliStabiliser() {}

  PauliStabiliser(const std::vector<Pauli> string, const bool coeff);

  /**
   * Determine whether two PauliStabilisers are equivalent
   *
   * @param other second PauliStabiliser
   *
   * @return true iff *this and other represent PauliStabiliser
   */
  bool operator==(const PauliStabiliser &other) const;

  /**
   * Determine whether two PauliStabilisers are not equivalent
   *
   * @param other second PauliStabiliser
   *
   * @return true iff *this and other represent different PauliStabiliser
   */
  bool operator!=(const PauliStabiliser &other) const;
};

typedef std::vector<PauliStabiliser> PauliStabiliserList;

JSON_DECL(PauliStabiliser)

/**
 * Calculate a sparse matrix corresponding to a sum of QubitPauliString.
 *
 * @param total_operator Operator specified by sum of pauli strings
 * @param n_qubits Number of qubits the string covers
 *
 * @return A complex sparse matrix corresponding to the n_qubit operator
 */
CmplxSpMat operator_tensor(
    const OperatorSum &total_operator, unsigned n_qubits);
CmplxSpMat operator_tensor(
    const OperatorSum &total_operator, const qubit_vector_t &qubits);
/**
 * Calculate expectation value of QubitPauliString with respect to a state.
 * <state|Operator|state>
 *
 * @param total_operator Operator specified by sum of pauli strings
 * @param state state, encoded by a complex vector
 *
 * @return Complex expectation value
 */
Complex operator_expectation(
    const OperatorSum &total_operator, const Eigen::VectorXcd &state);
Complex operator_expectation(
    const OperatorSum &total_operator, const Eigen::VectorXcd &state,
    const qubit_vector_t &qubits);

/**
 * A tensor of Pauli terms
 * \f$ P = i^k \sigma_1 \otimes \sigma_2 \otimes \cdots \otimes \sigma_n \f$.
 * This is a QubitPauliString but with an additional complex coefficient.
 */
class QubitPauliTensor {
 public:
  typedef std::map<std::pair<Pauli, Pauli>, std::pair<Complex, Pauli>>
      Mult_Matrix;
  static const Mult_Matrix &get_mult_matrix();

  QubitPauliString string;
  Complex coeff;

  /**
   * Construct an identity operator
   */
  QubitPauliTensor() : string(), coeff(1.) {}

  /**
   * Construct a constant multiple of the identity
   */
  explicit QubitPauliTensor(Complex _coeff) : string(), coeff(_coeff) {}

  /**
   * Construct a single Pauli term
   */
  QubitPauliTensor(const Qubit &qubit, Pauli p)
      : string({{qubit, p}}), coeff(1.) {}

  /**
   * Construct a tensor of many Pauli terms. Shortcut to make a
   * tensor for a default qubit register, without explicitly
   * constructing the QubitPauliMap
   *
   * @param _paulis list of Pauli letters
   */
  explicit QubitPauliTensor(const std::list<Pauli> &_paulis)
      : string({_paulis}), coeff(1.) {}

  /**
   * Construct a constant multiple of a single Pauli term
   */
  QubitPauliTensor(const Qubit &qubit, Pauli p, Complex _coeff)
      : string({{qubit, p}}), coeff(_coeff) {}

  /**
   * Construct a tensor product of many Pauli terms
   *
   * @param _string sparse representation of the full tensor
   */
  explicit QubitPauliTensor(const QubitPauliString &_string)
      : string(_string), coeff(1.) {}

  /**
   * Construct a tensor product of many Pauli terms
   *
   * @param _map sparse representation of the full tensor
   */
  explicit QubitPauliTensor(const QubitPauliMap &_map)
      : string(_map), coeff(1.) {}

  /**
   * Construct an arbitrary QubitPauliTensor
   *
   * @param _string sparse representation of the full tensor
   * @param _coeff complex coefficient
   */
  QubitPauliTensor(const QubitPauliString &_string, Complex _coeff)
      : string(_string), coeff(_coeff) {}

  /**
   * Construct an arbitrary QubitPauliTensor
   *
   * @param _map sparse representation of the full tensor
   * @param _coeff complex coefficient
   */
  QubitPauliTensor(const QubitPauliMap &_map, Complex _coeff)
      : string(_map), coeff(_coeff) {}

  /**
   * Calculate the product of two tensors
   *
   * @param other second tensor
   *
   * @return *this x other
   */
  QubitPauliTensor operator*(const QubitPauliTensor &other) const;

  /**
   * Determine whether two tensors are equivalent (ignoring I terms)
   *
   * @param other second tensor
   *
   * @return true iff *this and other represent the same numerical tensor
   */
  bool operator==(const QubitPauliTensor &other) const;

  /**
   * Determine whether two tensors are not equivalent (ignoring I terms)
   *
   * @param other second tensor
   *
   * @return true iff *this and other represent different numerical tensors
   */
  bool operator!=(const QubitPauliTensor &other) const;

  /**
   * Determine the lexicographic ordering of two tensors.
   * Ignores I terms.
   * Orders individual terms by the ordering of Qubits.
   * If two terms have equivalent tensors, they are ordered by coefficient
   * (lexicographically according to the pair <real, imag>).
   *
   * @param other second tensor
   *
   * @return true iff *this orders strictly before other
   */
  bool operator<(const QubitPauliTensor &other) const;

  /**
   * Removes I terms to compress the sparse representation.
   * Wrapper method for `QubitPauliString` method.
   */
  void compress();

  /**
   * Determines whether or not two tensor commute.
   * Wrapper method for `QubitPauliString` method.
   *
   * @param other Tensor to compare this to
   *
   * @return true if \p other commutes with this
   * @return false if \p other anti-commutes with this
   */
  bool commutes_with(const QubitPauliTensor &other) const;

  /**
   * Finds qubits with the same (non-trivial) Pauli term.
   * Wrapper method for `QubitPauliString` method.
   *
   * @param other Tensor to compare this to
   *
   * @return All qubits q where this->string.map[q] == other.string.map[q] != I
   */
  std::set<Qubit> common_qubits(const QubitPauliTensor &other) const;

  /**
   * Finds qubits that only occur in this tensor.
   * Wrapper method for `QubitPauliString` method.
   *
   * @param other Tensor to compare this to
   *
   * @return All qubits q where this->string.map[q] != I and other.string.map[q]
   * == I
   */
  std::set<Qubit> own_qubits(const QubitPauliTensor &other) const;

  /**
   * Finds qubits with different (non-trivial) Pauli terms.
   * Wrapper method for `QubitPauliString` method.
   *
   * @param other Tensor to compare this to
   *
   * @return All qubits q where I != this->string.map[q] != other.string.map[q]
   * != I
   */
  std::set<Qubit> conflicting_qubits(const QubitPauliTensor &other) const;

  /**
   * Readable string for sparse operator
   */
  std::string to_str() const;

  /**
   * Hash method, required for top-level python. Does not
   * distinguish between equivalent Strings, i.e. ignores I terms.
   */
  friend std::size_t hash_value(const QubitPauliTensor &qpt);
};

/**
 * Multiply coefficient by a scalar
 *
 * @param a scalar multiplier
 * @param qpt tensor
 *
 * @return result of the scalar multiplication
 */
QubitPauliTensor operator*(Complex a, const QubitPauliTensor &qpt);

}  // namespace tket
