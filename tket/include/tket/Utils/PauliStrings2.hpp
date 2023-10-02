// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tket/Utils/Constants.hpp"
#include "tket/Utils/EigenConfig.hpp"
#include "tket/Utils/Expression.hpp"
#include "tket/Utils/Json.hpp"
#include "tket/Utils/UnitID.hpp"

namespace tket {

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

typedef Eigen::SparseMatrix<Complex, Eigen::ColMajor> CmplxSpMat;

typedef std::map<Qubit, Pauli> QubitPauliMap;
typedef std::vector<Pauli> DensePauliMap;

struct no_coeff_t {};
typedef unsigned quarter_turns_t;

JSON_DECL(no_coeff_t)

template <typename CoeffType>
CoeffType default_coeff() = delete;

template <>
no_coeff_t default_coeff<no_coeff_t>();
template <>
quarter_turns_t default_coeff<quarter_turns_t>();
template <>
Complex default_coeff<Complex>();
template <>
Expr default_coeff<Expr>();

template <typename OriginalContainer, typename NewContainer>
NewContainer cast_container(const OriginalContainer &) = delete;

template <>
QubitPauliMap cast_container<QubitPauliMap, QubitPauliMap>(
    const QubitPauliMap &cont);
template <>
QubitPauliMap cast_container<DensePauliMap, QubitPauliMap>(
    const DensePauliMap &cont);
template <>
DensePauliMap cast_container<QubitPauliMap, DensePauliMap>(
    const QubitPauliMap &cont);
template <>
DensePauliMap cast_container<DensePauliMap, DensePauliMap>(
    const DensePauliMap &cont);

template <typename OriginalCoeff, typename NewCoeff>
NewCoeff cast_coeff(const OriginalCoeff &coeff) = delete;

template <>
no_coeff_t cast_coeff<no_coeff_t, no_coeff_t>(const no_coeff_t &);
template <>
quarter_turns_t cast_coeff<no_coeff_t, quarter_turns_t>(const no_coeff_t &);
template <>
Complex cast_coeff<no_coeff_t, Complex>(const no_coeff_t &);
template <>
Expr cast_coeff<no_coeff_t, Expr>(const no_coeff_t &);

template <>
no_coeff_t cast_coeff<quarter_turns_t, no_coeff_t>(const quarter_turns_t &);
template <>
quarter_turns_t cast_coeff<quarter_turns_t, quarter_turns_t>(
    const quarter_turns_t &coeff);
template <>
Complex cast_coeff<quarter_turns_t, Complex>(const quarter_turns_t &coeff);
template <>
Expr cast_coeff<quarter_turns_t, Expr>(const quarter_turns_t &coeff);

template <>
no_coeff_t cast_coeff<Complex, no_coeff_t>(const Complex &);
template <>
quarter_turns_t cast_coeff<Complex, quarter_turns_t>(const Complex &coeff);
template <>
Complex cast_coeff<Complex, Complex>(const Complex &coeff);
template <>
Expr cast_coeff<Complex, Expr>(const Complex &coeff);

template <>
no_coeff_t cast_coeff<Expr, no_coeff_t>(const Expr &);
template <>
quarter_turns_t cast_coeff<Expr, quarter_turns_t>(const Expr &coeff);
template <>
Complex cast_coeff<Expr, Complex>(const Expr &coeff);
template <>
Expr cast_coeff<Expr, Expr>(const Expr &coeff);

template <typename PauliContainer>
int compare_containers(
    const PauliContainer &first, const PauliContainer &second) = delete;

template <>
int compare_containers<QubitPauliMap>(
    const QubitPauliMap &first, const QubitPauliMap &second);
template <>
int compare_containers<DensePauliMap>(
    const DensePauliMap &first, const DensePauliMap &second);

template <typename CoeffType>
int compare_coeffs(const CoeffType &first, const CoeffType &second) = delete;

template <>
int compare_coeffs<no_coeff_t>(const no_coeff_t &, const no_coeff_t &);
template <>
int compare_coeffs<quarter_turns_t>(
    const quarter_turns_t &first, const quarter_turns_t &second);
template <>
int compare_coeffs<Complex>(const Complex &first, const Complex &second);
template <>
int compare_coeffs<Expr>(const Expr &first, const Expr &second);

std::set<Qubit> common_qubits(
    const QubitPauliMap &first, const QubitPauliMap &second);

std::set<unsigned> common_indices(
    const DensePauliMap &first, const DensePauliMap &second);

std::set<Qubit> own_qubits(
    const QubitPauliMap &first, const QubitPauliMap &second);

std::set<unsigned> own_indices(
    const DensePauliMap &first, const DensePauliMap &second);

std::set<Qubit> conflicting_qubits(
    const QubitPauliMap &first, const QubitPauliMap &second);

std::set<unsigned> conflicting_indices(
    const DensePauliMap &first, const DensePauliMap &second);

template <typename PauliContainer>
bool commuting_containers(
    const PauliContainer &first, const PauliContainer &second) = delete;

template <>
bool commuting_containers<QubitPauliMap>(
    const QubitPauliMap &first, const QubitPauliMap &second);
template <>
bool commuting_containers<DensePauliMap>(
    const DensePauliMap &first, const DensePauliMap &second);

template <typename PauliContainer>
void print_paulis(std::ostream &os, const PauliContainer &paulis) = delete;

template <>
void print_paulis<QubitPauliMap>(std::ostream &os, const QubitPauliMap &paulis);
template <>
void print_paulis<DensePauliMap>(std::ostream &os, const DensePauliMap &paulis);

template <typename CoeffType>
void print_coeff(std::ostream &os, const CoeffType &coeff) = delete;

template <>
void print_coeff<no_coeff_t>(std::ostream &, const no_coeff_t &);
template <>
void print_coeff<quarter_turns_t>(
    std::ostream &os, const quarter_turns_t &coeff);
template <>
void print_coeff<Complex>(std::ostream &os, const Complex &coeff);
template <>
void print_coeff<Expr>(std::ostream &os, const Expr &coeff);

template <typename PauliContainer>
void hash_combine_paulis(std::size_t &seed, const PauliContainer &paulis) =
    delete;

template <>
void hash_combine_paulis<QubitPauliMap>(
    std::size_t &seed, const QubitPauliMap &paulis);
template <>
void hash_combine_paulis<DensePauliMap>(
    std::size_t &seed, const DensePauliMap &paulis);

template <typename CoeffType>
void hash_combine_coeff(std::size_t &seed, const CoeffType &coeff) = delete;

template <>
void hash_combine_coeff<no_coeff_t>(std::size_t &, const no_coeff_t &);
template <>
void hash_combine_coeff<quarter_turns_t>(
    std::size_t &seed, const quarter_turns_t &coeff);
template <>
void hash_combine_coeff<Complex>(std::size_t &seed, const Complex &coeff);
template <>
void hash_combine_coeff<Expr>(std::size_t &seed, const Expr &coeff);

template <typename PauliContainer>
unsigned n_ys(const PauliContainer &paulis) = delete;

template <>
unsigned n_ys<QubitPauliMap>(const QubitPauliMap &paulis);
template <>
unsigned n_ys<DensePauliMap>(const DensePauliMap &paulis);

const std::map<std::pair<Pauli, Pauli>, std::pair<quarter_turns_t, Pauli>> &
get_mult_matrix();

template <typename PauliContainer>
std::pair<quarter_turns_t, PauliContainer> multiply_strings(
    const PauliContainer &first, const PauliContainer &second) = delete;

template <>
std::pair<quarter_turns_t, QubitPauliMap> multiply_strings<QubitPauliMap>(
    const QubitPauliMap &first, const QubitPauliMap &second);
template <>
std::pair<quarter_turns_t, DensePauliMap> multiply_strings<DensePauliMap>(
    const DensePauliMap &first, const DensePauliMap &second);

template <typename CoeffType>
CoeffType multiply_coeffs(const CoeffType &first, const CoeffType &second) =
    delete;

template <>
no_coeff_t multiply_coeffs<no_coeff_t>(const no_coeff_t &, const no_coeff_t &);
template <>
quarter_turns_t multiply_coeffs<quarter_turns_t>(
    const quarter_turns_t &first, const quarter_turns_t &second);
template <>
Complex multiply_coeffs<Complex>(const Complex &first, const Complex &second);
template <>
Expr multiply_coeffs<Expr>(const Expr &first, const Expr &second);

template <typename PauliContainer>
CmplxSpMat to_sparse_matrix(const PauliContainer &paulis) = delete;

template <>
CmplxSpMat to_sparse_matrix<QubitPauliMap>(const QubitPauliMap &paulis);
template <>
CmplxSpMat to_sparse_matrix<DensePauliMap>(const DensePauliMap &paulis);

template <typename PauliContainer>
CmplxSpMat to_sparse_matrix(const PauliContainer &paulis, unsigned n_qubits) =
    delete;

template <>
CmplxSpMat to_sparse_matrix<QubitPauliMap>(
    const QubitPauliMap &paulis, unsigned n_qubits);
template <>
CmplxSpMat to_sparse_matrix<DensePauliMap>(
    const DensePauliMap &paulis, unsigned n_qubits);

template <typename PauliContainer>
CmplxSpMat to_sparse_matrix(
    const PauliContainer &paulis, const qubit_vector_t &qubits) = delete;

template <>
CmplxSpMat to_sparse_matrix<QubitPauliMap>(
    const QubitPauliMap &paulis, const qubit_vector_t &qubits);
template <>
CmplxSpMat to_sparse_matrix<DensePauliMap>(
    const DensePauliMap &paulis, const qubit_vector_t &qubits);

template <typename PauliContainer, typename CoeffType>
class PauliTensor {
  static_assert(
      std::is_same<PauliContainer, QubitPauliMap>::value ||
          std::is_same<PauliContainer, DensePauliMap>::value,
      "PauliTensor must be either dense or qubit-indexed.");
  static_assert(
      std::is_same<CoeffType, no_coeff_t>::value ||
          std::is_same<CoeffType, quarter_turns_t>::value ||
          std::is_same<CoeffType, Complex>::value ||
          std::is_same<CoeffType, Expr>::value,
      "PauliTensor must either support no coefficient, quarter turns, or "
      "arbitrary complex (floating point or symbolic) coefficients.");

  static const CoeffType default_coeff;

 public:
  PauliContainer string;
  CoeffType coeff;

  PauliTensor() : string(), coeff(default_coeff) {}

  PauliTensor(
      const PauliContainer &_string, const CoeffType &_coeff = default_coeff)
      : string(_string), coeff(_coeff) {}

  template <typename PC = PauliContainer>
  PauliTensor(
      const Qubit &q, Pauli p, const CoeffType &_coeff = default_coeff,
      typename std::enable_if<std::is_same<PC, QubitPauliMap>::value>::type * =
          0)
      : string({{q, p}}), coeff(_coeff) {}

  template <typename PC = PauliContainer>
  PauliTensor(
      const DensePauliMap &_string, const CoeffType &_coeff = default_coeff,
      typename std::enable_if<std::is_same<PC, QubitPauliMap>::value>::type * =
          0)
      : string(cast_container<DensePauliMap, QubitPauliMap>(_string)),
        coeff(_coeff) {}

  template <typename CastContainer, typename CastCoeffType>
  operator PauliTensor<CastContainer, CastCoeffType>() const {
    return PauliTensor<CastContainer, CastCoeffType>(
        cast_container<PauliContainer, CastContainer>(string),
        cast_coeff<CoeffType, CastCoeffType>(coeff));
  }

  int compare(const PauliTensor<PauliContainer, CoeffType> &other) const {
    int cont_comp =
        compare_containers<PauliContainer>(this->string, other.string);
    if (cont_comp != 0) return cont_comp;
    return compare_coeffs<CoeffType>(this->coeff, other.coeff);
  }

  bool operator==(const PauliTensor<PauliContainer, CoeffType> &other) const {
    return (compare(other) == 0);
  }

  bool operator!=(const PauliTensor<PauliContainer, CoeffType> &other) const {
    return !(*this == other);
  }

  bool operator<(const PauliTensor<PauliContainer, CoeffType> &other) const {
    return compare(other) < 0;
  }

  template <typename PC = PauliContainer>
  typename std::enable_if<std::is_same<PC, QubitPauliMap>::value>::type
  compress() {
    QubitPauliMap::iterator it = string.begin();
    while (it != string.end()) {
      if (it->second == Pauli::I)
        it = string.erase(it);
      else
        ++it;
    }
  }

  template <typename OtherCoeffType>
  bool commutes_with(
      const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return commuting_containers<PauliContainer>(string, other.string);
  }

  template <typename OtherCoeffType, typename PC = PauliContainer>
  typename std::enable_if<
      std::is_same<PC, QubitPauliMap>::value, std::set<Qubit>>::type
  common_qubits(
      const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return tket::common_qubits(string, other.string);
  }

  template <typename OtherCoeffType, typename PC = PauliContainer>
  typename std::enable_if<
      std::is_same<PC, QubitPauliMap>::value, std::set<Qubit>>::type
  own_qubits(const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return tket::own_qubits(string, other.string);
  }

  template <typename OtherCoeffType, typename PC = PauliContainer>
  typename std::enable_if<
      std::is_same<PC, QubitPauliMap>::value, std::set<Qubit>>::type
  conflicting_qubits(
      const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return tket::conflicting_qubits(string, other.string);
  }

  template <typename OtherCoeffType, typename PC = PauliContainer>
  typename std::enable_if<
      std::is_same<PC, DensePauliMap>::value, std::set<unsigned>>::type
  common_indices(
      const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return tket::common_indices(string, other.string);
  }

  template <typename OtherCoeffType, typename PC = PauliContainer>
  typename std::enable_if<
      std::is_same<PC, DensePauliMap>::value, std::set<unsigned>>::type
  own_indices(const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return tket::own_indices(string, other.string);
  }

  template <typename OtherCoeffType, typename PC = PauliContainer>
  typename std::enable_if<
      std::is_same<PC, DensePauliMap>::value, std::set<unsigned>>::type
  conflicting_indices(
      const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return tket::conflicting_indices(string, other.string);
  }

  template <typename PC = PauliContainer>
  typename std::enable_if<std::is_same<PC, QubitPauliMap>::value, Pauli>::type
  get(const Qubit &qb) const {
    QubitPauliMap::const_iterator i = string.find(qb);
    if (i == string.end())
      return Pauli::I;
    else
      return i->second;
  }
  template <typename PC = PauliContainer>
  typename std::enable_if<std::is_same<PC, DensePauliMap>::value, Pauli>::type
  get(unsigned qb) const {
    if (qb >= string.size())
      return Pauli::I;
    else
      return string.at(qb);
  }

  template <typename PC = PauliContainer>
  typename std::enable_if<std::is_same<PC, QubitPauliMap>::value>::type set(
      const Qubit &qb, Pauli p) {
    QubitPauliMap::iterator i = string.find(qb);
    if (i == string.end()) {
      if (p != Pauli::I) string.insert({qb, p});
    } else {
      if (p == Pauli::I)
        string.erase(i);
      else
        i->second = p;
    }
  }
  template <typename PC = PauliContainer>
  typename std::enable_if<std::is_same<PC, DensePauliMap>::value>::type set(
      unsigned qb, Pauli p) {
    if (qb >= string.size()) string.resize(qb + 1, Pauli::I);
    string.at(qb) = p;
  }

  std::string to_str() const {
    std::stringstream ss;
    print_coeff<CoeffType>(ss, coeff);
    print_paulis<PauliContainer>(ss, string);
    return ss.str();
  }

  std::size_t hash_value() const {
    std::size_t seed = 0;
    hash_combine_paulis<PauliContainer>(seed, string);
    hash_combine_coeff<CoeffType>(seed, coeff);
    return seed;
  }

  void transpose() {
    if (n_ys<PauliContainer>(string) % 2 == 1)
      coeff = multiply_coeffs<CoeffType>(
          coeff, cast_coeff<quarter_turns_t, CoeffType>(2));
  }

  PauliTensor<PauliContainer, CoeffType> operator*(
      const PauliTensor<PauliContainer, CoeffType> &other) const {
    std::pair<quarter_turns_t, PauliContainer> prod =
        multiply_strings<PauliContainer>(this->string, other.string);
    CoeffType new_coeff = multiply_coeffs<CoeffType>(this->coeff, other.coeff);
    new_coeff = multiply_coeffs<CoeffType>(
        new_coeff, cast_coeff<quarter_turns_t, CoeffType>(prod.first));
    return PauliTensor<PauliContainer, CoeffType>(prod.second, new_coeff);
  }

  template <typename CT = CoeffType>
  typename std::enable_if<std::is_same<CT, Expr>::value, SymSet>::type
  free_symbols() const {
    return expr_free_symbols(coeff);
  }
  template <typename CT = CoeffType>
  typename std::enable_if<
      std::is_same<CT, Expr>::value,
      PauliTensor<PauliContainer, CoeffType>>::type
  symbol_substitution(const SymEngine::map_basic_basic &sub_map) const {
    return PauliTensor<PauliContainer, CoeffType>(string, coeff.subs(sub_map));
  }

  unsigned size() const { return string.size(); }

  CmplxSpMat to_sparse_matrix() const {
    return cast_coeff<CoeffType, Complex>(coeff) *
           tket::to_sparse_matrix<PauliContainer>(string);
  }
  CmplxSpMat to_sparse_matrix(const unsigned n_qubits) const {
    return cast_coeff<CoeffType, Complex>(coeff) *
           tket::to_sparse_matrix<PauliContainer>(string, n_qubits);
  }
  CmplxSpMat to_sparse_matrix(const qubit_vector_t &qubits) const {
    return cast_coeff<CoeffType, Complex>(coeff) *
           tket::to_sparse_matrix<PauliContainer>(string, qubits);
  }
};

template <typename PauliContainer, typename CoeffType>
const CoeffType PauliTensor<PauliContainer, CoeffType>::default_coeff =
    tket::default_coeff<CoeffType>();

template <typename PauliContainer, typename CoeffType>
void to_json(
    nlohmann::json &j, const PauliTensor<PauliContainer, CoeffType> &tensor) {
  j["string"] = tensor.string;
  j["coeff"] = tensor.coeff;
}

template <typename PauliContainer, typename CoeffType>
void from_json(
    const nlohmann::json &j, PauliTensor<PauliContainer, CoeffType> &tensor) {
  tensor = PauliTensor<PauliContainer, CoeffType>(
      j.at("string").get<PauliContainer>(), j.at("coeff").get<CoeffType>());
}

typedef PauliTensor<QubitPauliMap, no_coeff_t> SpPauliString;
typedef PauliTensor<DensePauliMap, no_coeff_t> PauliString;
typedef PauliTensor<QubitPauliMap, quarter_turns_t> SpPauliStabiliser;
typedef PauliTensor<DensePauliMap, quarter_turns_t> PauliStabiliser;
typedef PauliTensor<QubitPauliMap, Complex> SpCxPauliTensor;
typedef PauliTensor<DensePauliMap, Complex> CxPauliTensor;
typedef PauliTensor<QubitPauliMap, Expr> SpSymPauliTensor;
typedef PauliTensor<DensePauliMap, Expr> SymPauliTensor;

typedef std::vector<PauliStabiliser> PauliStabiliserVec;

}  // namespace tket
