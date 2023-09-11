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

typedef std::map<Qubit, Pauli> QubitPauliMap;
typedef std::vector<Pauli> DensePauliMap;

struct no_coeff_t {};
typedef unsigned quarter_turns_t;

template <typename CoeffType>
CoeffType default_coeff() = delete;

template <>
no_coeff_t default_coeff<no_coeff_t>() {
  return {};
}
template <>
quarter_turns_t default_coeff<quarter_turns_t>() {
  return 0;
}
template <>
Complex default_coeff<Complex>() {
  return 1.;
}
template <>
Expr default_coeff<Expr>() {
  return 1;
}

template <typename OriginalContainer, typename NewContainer>
NewContainer cast_container(const OriginalContainer &) = delete;

template <>
QubitPauliMap cast_container<QubitPauliMap, QubitPauliMap>(
    const QubitPauliMap &cont) {
  return cont;
}

template <>
QubitPauliMap cast_container<DensePauliMap, QubitPauliMap>(
    const DensePauliMap &cont) {
  QubitPauliMap res;
  for (unsigned i = 0; i < cont.size(); ++i) {
    Pauli pi = cont.at(i);
    if (pi != Pauli::I) res.insert({Qubit(i), cont.at(i)});
  }
  return res;
}

template <>
DensePauliMap cast_container<QubitPauliMap, DensePauliMap>(
    const QubitPauliMap &cont) {
  if (cont.empty()) return {};
  unsigned max_index = 0;
  for (const std::pair<const Qubit, Pauli> &pair : cont) {
    if (pair.first.reg_info() != register_info_t{UnitType::Qubit, 1} ||
        pair.first.reg_name() != q_default_reg())
      throw std::logic_error(
          "Cannot cast a QubitPauliMap with non-default register qubits to a "
          "DensePauliMap");
    unsigned i = pair.first.index().front();
    if (i > max_index) max_index = i;
  }
  DensePauliMap res(max_index + 1, Pauli::I);
  for (const std::pair<const Qubit, Pauli> &pair : cont) {
    res.at(pair.first.index().front()) = pair.second;
  }
  return res;
}

template <>
DensePauliMap cast_container<DensePauliMap, DensePauliMap>(
    const DensePauliMap &cont) {
  return cont;
}

template <typename OriginalCoeff, typename NewCoeff>
NewCoeff cast_coeff(const OriginalCoeff &coeff) = delete;

template <>
no_coeff_t cast_coeff<no_coeff_t, no_coeff_t>(const no_coeff_t &) {
  return {};
}
template <>
quarter_turns_t cast_coeff<no_coeff_t, quarter_turns_t>(const no_coeff_t &) {
  return 0;
}
template <>
Complex cast_coeff<no_coeff_t, Complex>(const no_coeff_t &) {
  return 1.;
}
template <>
Expr cast_coeff<no_coeff_t, Expr>(const no_coeff_t &) {
  return 1.;
}

template <>
no_coeff_t cast_coeff<quarter_turns_t, no_coeff_t>(const quarter_turns_t &) {
  return {};
}
template <>
quarter_turns_t cast_coeff<quarter_turns_t, quarter_turns_t>(
    const quarter_turns_t &coeff) {
  return coeff;
}
template <>
Complex cast_coeff<quarter_turns_t, Complex>(const quarter_turns_t &coeff) {
  switch (coeff % 4) {
    case 0: {
      return 1.;
    }
    case 1: {
      return i_;
    }
    case 2: {
      return -1.;
    }
    default: {
      return -i_;
    }
  }
}
template <>
Expr cast_coeff<quarter_turns_t, Expr>(const quarter_turns_t &coeff) {
  switch (coeff % 4) {
    case 0: {
      return 1.;
    }
    case 1: {
      return Expr(SymEngine::I);
    }
    case 2: {
      return -1.;
    }
    default: {
      return -Expr(SymEngine::I);
    }
  }
}

template <>
no_coeff_t cast_coeff<Complex, no_coeff_t>(const Complex &) {
  return {};
}
template <>
quarter_turns_t cast_coeff<Complex, quarter_turns_t>(const Complex &coeff) {
  if (std::abs(coeff - 1.) < EPS)
    return 0;
  else if (std::abs(coeff - i_) < EPS)
    return 1;
  else if (std::abs(coeff + 1.) < EPS)
    return 2;
  else if (std::abs(coeff + i_) < EPS)
    return 3;
  else
    throw std::logic_error(
        "Could not cast PauliTensor coefficient to quarter turns: not a power "
        "of i.");
}
template <>
Complex cast_coeff<Complex, Complex>(const Complex &coeff) {
  return coeff;
}
template <>
Expr cast_coeff<Complex, Expr>(const Complex &coeff) {
  return Expr(SymEngine::make_rcp<const SymEngine::ComplexDouble>(coeff));
}

template <>
no_coeff_t cast_coeff<Expr, no_coeff_t>(const Expr &) {
  return {};
}
template <>
quarter_turns_t cast_coeff<Expr, quarter_turns_t>(const Expr &coeff) {
  std::optional<Complex> ev = eval_expr_c(coeff);
  if (ev)
    return cast_coeff<Complex, quarter_turns_t>(*ev);
  else
    throw std::logic_error(
        "Could not cast symbolic PauliTensor to quarter turns.");
}
template <>
Complex cast_coeff<Expr, Complex>(const Expr &coeff) {
  std::optional<Complex> ev = eval_expr_c(coeff);
  if (ev)
    return *ev;
  else
    throw std::logic_error(
        "Could not cast symbolic PauliTensor to complex coefficient.");
}
template <>
Expr cast_coeff<Expr, Expr>(const Expr &coeff) {
  return coeff;
}

template <typename PauliContainer>
int compare_containers(
    const PauliContainer &first, const PauliContainer &second) = delete;

template <>
int compare_containers<QubitPauliMap>(
    const QubitPauliMap &first, const QubitPauliMap &second) {
  QubitPauliMap::const_iterator p1_it = first.begin();
  QubitPauliMap::const_iterator p2_it = second.begin();
  while (p1_it != first.end()) {
    if (p1_it->second == Pauli::I) {
      ++p1_it;
      continue;
    }
    while (p2_it != second.end() && p2_it->second == Pauli::I) {
      ++p2_it;
    }
    if (p2_it == second.end()) return 1;
    // QubitPauliString order should reflect ILO
    // i.e. IZ < ZI (Zq1 < Zq0)
    // Hence we first order by reverse of leading qubit
    if (p1_it->first < p2_it->first) return 1;
    if (p2_it->first < p1_it->first) return -1;
    // and then by increasing order of Pauli letter on the same qubit
    if (p1_it->second < p2_it->second) return -1;
    if (p1_it->second > p2_it->second) return 1;
    ++p1_it;
    ++p2_it;
  }
  while (p2_it != second.end() && p2_it->second == Pauli::I) {
    ++p2_it;
  }
  return (p2_it == second.end()) ? 0 : -1;
}

template <>
int compare_containers<DensePauliMap>(
    const DensePauliMap &first, const DensePauliMap &second) {
  DensePauliMap::const_iterator p1_it = first.begin();
  DensePauliMap::const_iterator p2_it = second.begin();
  while (p1_it != first.end() && p2_it != second.end()) {
    if (*p1_it < *p2_it) return -1;
    if (*p1_it > *p2_it) return 1;
    ++p1_it;
    ++p2_it;
  }
  while (p1_it != first.end() && *p1_it == Pauli::I) ++p1_it;
  if (p1_it != first.end()) return 1;
  while (p2_it != second.end() && *p2_it == Pauli::I) ++p2_it;
  return (p2_it == second.end()) ? 0 : -1;
}

template <typename CoeffType>
int compare_coeffs(const CoeffType &first, const CoeffType &second) = delete;

template <>
int compare_coeffs<no_coeff_t>(const no_coeff_t &, const no_coeff_t &) {
  return 0;
}
template <>
int compare_coeffs<quarter_turns_t>(
    const quarter_turns_t &first, const quarter_turns_t &second) {
  if (first % 4 < second % 4) return -1;
  return (first % 4 == second % 4) ? 0 : 1;
}
template <>
int compare_coeffs<Complex>(const Complex &first, const Complex &second) {
  if (first.real() < second.real()) return -1;
  if (first.real() > second.real()) return 1;
  if (first.imag() < second.imag()) return -1;
  return (first.imag() == second.imag()) ? 0 : 1;
}
template <>
int compare_coeffs<Expr>(const Expr &first, const Expr &second) {
  return first.get_basic()->compare(second);
}

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
    const QubitPauliMap &first, const QubitPauliMap &second) {
  return (conflicting_qubits(first, second).size() % 2) == 0;
}

template <>
bool commuting_containers<DensePauliMap>(
    const DensePauliMap &first, const DensePauliMap &second) {
  return (conflicting_indices(first, second).size() % 2) == 0;
}

template <typename PauliContainer>
void print_paulis(std::ostream &os, const PauliContainer &paulis) = delete;

template <>
void print_paulis<QubitPauliMap>(
    std::ostream &os, const QubitPauliMap &paulis) {
  os << "(";
  QubitPauliMap::const_iterator i = paulis.begin();
  while (i != paulis.end()) {
    switch (i->second) {
      case Pauli::I: {
        os << "I";
        break;
      }
      case Pauli::X: {
        os << "X";
        break;
      }
      case Pauli::Y: {
        os << "Y";
        break;
      }
      case Pauli::Z: {
        os << "Z";
        break;
      }
    }
    os << i->first.repr();
    ++i;
    if (i != paulis.end()) os << ", ";
  }
  os << ")";
}

template <>
void print_paulis<DensePauliMap>(
    std::ostream &os, const DensePauliMap &paulis) {
  for (const Pauli &p : paulis) {
    switch (p) {
      case Pauli::I: {
        os << "I";
        break;
      }
      case Pauli::X: {
        os << "X";
        break;
      }
      case Pauli::Y: {
        os << "Y";
        break;
      }
      case Pauli::Z: {
        os << "Z";
        break;
      }
    }
  }
}

template <typename CoeffType>
void print_coeff(std::ostream &os, const CoeffType &coeff) = delete;

template <>
void print_coeff<no_coeff_t>(std::ostream &, const no_coeff_t &) {}

template <>
void print_coeff<quarter_turns_t>(
    std::ostream &os, const quarter_turns_t &coeff) {
  switch (coeff % 4) {
    case 1: {
      os << "i*";
      break;
    }
    case 2: {
      os << "-";
      break;
    }
    case 3: {
      os << "-i*";
      break;
    }
    default: {
      break;
    }
  }
}

template <>
void print_coeff<Complex>(std::ostream &os, const Complex &coeff) {
  if (coeff == -1.) {
    os << "-";
  } else if (coeff != 1.) {
    os << coeff << "*";
  }
}

template <>
void print_coeff<Expr>(std::ostream &os, const Expr &coeff) {
  if (coeff == -1.) {
    os << "-";
  } else if (coeff != -1.) {
    os << coeff << "*";
  }
}

template <typename PauliContainer>
void hash_combine_paulis(std::size_t &seed, const PauliContainer &paulis) =
    delete;

template <>
void hash_combine_paulis<QubitPauliMap>(
    std::size_t &seed, const QubitPauliMap &paulis) {
  for (const std::pair<const Qubit, Pauli> &qp : paulis) {
    if (qp.second != Pauli::I) {
      boost::hash_combine(seed, qp.first);
      boost::hash_combine(seed, qp.second);
    }
  }
}

template <>
void hash_combine_paulis<DensePauliMap>(
    std::size_t &seed, const DensePauliMap &paulis) {
  DensePauliMap::const_reverse_iterator i = paulis.rbegin();
  while (i != paulis.rend() && *i == Pauli::I) {
    ++i;
  }
  while (i != paulis.rend()) {
    boost::hash_combine(seed, *i);
    ++i;
  }
}

template <typename CoeffType>
void hash_combine_coeff(std::size_t &seed, const CoeffType &coeff) = delete;

template <>
void hash_combine_coeff<no_coeff_t>(std::size_t &, const no_coeff_t &) {}

template <>
void hash_combine_coeff<quarter_turns_t>(
    std::size_t &seed, const quarter_turns_t &coeff) {
  boost::hash_combine(seed, coeff % 4);
}

template <>
void hash_combine_coeff<Complex>(std::size_t &seed, const Complex &coeff) {
  boost::hash_combine(seed, coeff);
}

template <>
void hash_combine_coeff<Expr>(std::size_t &seed, const Expr &coeff) {
  boost::hash_combine(seed, coeff.get_basic()->hash());
}

template <typename PauliContainer>
unsigned n_ys(const PauliContainer &paulis) = delete;

template <>
unsigned n_ys<QubitPauliMap>(const QubitPauliMap &paulis) {
  unsigned n = 0;
  for (const std::pair<const Qubit, Pauli> &qp : paulis) {
    if (qp.second == Pauli::Y) ++n;
  }
  return n;
}

template <>
unsigned n_ys<DensePauliMap>(const DensePauliMap &paulis) {
  unsigned n = 0;
  for (const Pauli &p : paulis) {
    if (p == Pauli::Y) ++n;
  }
  return n;
}

const std::map<std::pair<Pauli, Pauli>, std::pair<quarter_turns_t, Pauli>> &
get_mult_matrix() {
  static const std::map<
      std::pair<Pauli, Pauli>, std::pair<quarter_turns_t, Pauli>>
      mult_matrix{
          {{Pauli::I, Pauli::I}, {0, Pauli::I}},
          {{Pauli::I, Pauli::X}, {0, Pauli::X}},
          {{Pauli::I, Pauli::Y}, {0, Pauli::Y}},
          {{Pauli::I, Pauli::Z}, {0, Pauli::Z}},
          {{Pauli::X, Pauli::I}, {0, Pauli::X}},
          {{Pauli::X, Pauli::X}, {0, Pauli::I}},
          {{Pauli::X, Pauli::Y}, {1, Pauli::Z}},
          {{Pauli::X, Pauli::Z}, {3, Pauli::Y}},
          {{Pauli::Y, Pauli::I}, {0, Pauli::Y}},
          {{Pauli::Y, Pauli::X}, {3, Pauli::Z}},
          {{Pauli::Y, Pauli::Y}, {0, Pauli::I}},
          {{Pauli::Y, Pauli::Z}, {1, Pauli::X}},
          {{Pauli::Z, Pauli::I}, {0, Pauli::Z}},
          {{Pauli::Z, Pauli::X}, {1, Pauli::Y}},
          {{Pauli::Z, Pauli::Y}, {3, Pauli::X}},
          {{Pauli::Z, Pauli::Z}, {0, Pauli::I}},
      };
  return mult_matrix;
}

template <typename PauliContainer>
std::pair<quarter_turns_t, PauliContainer> multiply_strings(
    const PauliContainer &first, const PauliContainer &second) = delete;

template <>
std::pair<quarter_turns_t, QubitPauliMap> multiply_strings<QubitPauliMap>(
    const QubitPauliMap &first, const QubitPauliMap &second) {
  quarter_turns_t total_turns = 0;
  QubitPauliMap result;
  QubitPauliMap::const_iterator fi = first.begin();
  QubitPauliMap::const_iterator si = second.begin();
  while (fi != first.end()) {
    while (si != second.end() && si->first < fi->first) {
      result.insert(*si);
      ++si;
    }
    if (si != second.end() && si->first == fi->first) {
      // Pauli in the same position, so need to multiply
      const std::pair<quarter_turns_t, Pauli> &prod =
          get_mult_matrix().at({fi->second, si->second});
      total_turns += prod.first;
      if (prod.second != Pauli::I) {
        result.insert({fi->first, prod.second});
      }
      ++si;
    } else {
      result.insert(*fi);
    }
    ++fi;
  }
  while (si != second.end()) {
    result.insert(*si);
    ++si;
  }
  return {total_turns, result};
}

template <>
std::pair<quarter_turns_t, DensePauliMap> multiply_strings<DensePauliMap>(
    const DensePauliMap &first, const DensePauliMap &second) {
  quarter_turns_t total_turns = 0;
  DensePauliMap result;
  DensePauliMap::const_iterator fi = first.begin();
  DensePauliMap::const_iterator si = second.begin();
  while (fi != first.end() && si != second.end()) {
    const std::pair<quarter_turns_t, Pauli> &prod =
        get_mult_matrix().at({*fi, *si});
    total_turns += prod.first;
    result.push_back(prod.second);
    ++fi;
    ++si;
  }
  while (fi != first.end()) {
    result.push_back(*fi);
    ++fi;
  }
  while (si != second.end()) {
    result.push_back(*si);
    ++si;
  }
  return {total_turns, result};
}

template <typename CoeffType>
CoeffType multiply_coeffs(const CoeffType &first, const CoeffType &second) =
    delete;

template <>
no_coeff_t multiply_coeffs<no_coeff_t>(const no_coeff_t &, const no_coeff_t &) {
  return {};
}

template <>
quarter_turns_t multiply_coeffs<quarter_turns_t>(
    const quarter_turns_t &first, const quarter_turns_t &second) {
  return (first + second) % 4;
}

template <>
Complex multiply_coeffs<Complex>(const Complex &first, const Complex &second) {
  return first * second;
}

template <>
Expr multiply_coeffs<Expr>(const Expr &first, const Expr &second) {
  return first * second;
}

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

  static const CoeffType default_coeff = tket::default_coeff<CoeffType>();

  PauliContainer string;
  CoeffType coeff;

  PauliTensor(const CoeffType &_coeff = default_coeff)
      : string(), coeff(_coeff) {}

  PauliTensor(const std::initializer_list<Pauli> &paulis);

  PauliTensor(
      const std::list<Pauli> &paulis, const CoeffType &_coeff = default_coeff);

  PauliTensor(
      const PauliContainer &_string, const CoeffType &_coeff = default_coeff)
      : string(_string), coeff(_coeff) {}

  template <typename CastContainer, typename CastCoeffType>
  operator PauliTensor<CastContainer, CastCoeffType>() const {
    return PauliTensor(
        cast_container<PauliContainer, CastContainer>(string),
        cast_coeff<CoeffType, CastCoeffType>(coeff));
  }

  int compare(const PauliTensor<PauliContainer, CoeffType> &other) const {
    int coeff_comp = compare_coeffs<CoeffType>(this->coeff, other.coeff);
    if (coeff_comp != 0) return coeff_comp;
    return compare_containers<PauliContainer>(this->string, other.string);
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

  std::enable_if<std::is_same<PauliContainer, QubitPauliMap>::value, void>
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

  template <typename OtherCoeffType>
  std::enable_if<
      std::is_same<PauliContainer, QubitPauliMap>::value, std::set<Qubit>>
  common_qubits(
      const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return common_qubits(string, other.string);
  }

  template <typename OtherCoeffType>
  std::enable_if<
      std::is_same<PauliContainer, QubitPauliMap>::value, std::set<Qubit>>
  own_qubits(const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return own_qubits(string, other.string);
  }

  template <typename OtherCoeffType>
  std::enable_if<
      std::is_same<PauliContainer, QubitPauliMap>::value, std::set<Qubit>>
  conflicting_qubits(
      const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return conflicting_qubits(string, other.string);
  }

  template <typename OtherCoeffType>
  std::enable_if<
      std::is_same<PauliContainer, DensePauliMap>::value, std::set<unsigned>>
  common_indices(
      const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return common_indices(string, other.string);
  }

  template <typename OtherCoeffType>
  std::enable_if<
      std::is_same<PauliContainer, DensePauliMap>::value, std::set<unsigned>>
  own_indices(const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return own_indices(string, other.string);
  }

  template <typename OtherCoeffType>
  std::enable_if<
      std::is_same<PauliContainer, DensePauliMap>::value, std::set<unsigned>>
  conflicting_indices(
      const PauliTensor<PauliContainer, OtherCoeffType> &other) const {
    return conflicting_indices(string, other.string);
  }

  template <typename OtherCoeffType>
  std::enable_if<std::is_same<PauliContainer, QubitPauliMap>::value, Pauli> get(
      const Qubit &qb) const {
    QubitPauliMap::const_iterator i = string.find(qb);
    if (i == string.end())
      return Pauli::I;
    else
      return i->second;
  }
  template <typename OtherCoeffType>
  std::enable_if<std::is_same<PauliContainer, DensePauliMap>::value, Pauli> get(
      unsigned qb) const {
    if (qb >= string.size())
      return Pauli::I;
    else
      return string.at(qb);
  }

  template <typename OtherCoeffType>
  std::enable_if<std::is_same<PauliContainer, QubitPauliMap>::value, void> set(
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
  template <typename OtherCoeffType>
  std::enable_if<std::is_same<PauliContainer, DensePauliMap>::value, void> set(
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
};

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
typedef PauliTensor<QubitPauliMap, Complex> SpPauliTensor;
typedef PauliTensor<DensePauliMap, Complex> PauliTensor;
typedef PauliTensor<QubitPauliMap, Expr> SpSymPauliTensor;
typedef PauliTensor<DensePauliMap, Expr> SymPauliTensor;

}  // namespace tket
