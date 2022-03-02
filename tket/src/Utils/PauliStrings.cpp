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

#include "PauliStrings.hpp"

#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "Utils/Assert.hpp"
#include "Utils/Constants.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/Json.hpp"

namespace tket {

static const CmplxSpMat const_2x2_matrix(
    Complex tl, Complex tr, Complex bl, Complex br) {
  CmplxSpMat m(2, 2);
  if (tl != czero) {
    m.insert(0, 0) = tl;
  }
  if (tr != czero) {
    m.insert(0, 1) = tr;
  }
  if (bl != czero) {
    m.insert(1, 0) = bl;
  }
  if (br != czero) {
    m.insert(1, 1) = br;
  }
  return m;
}

static const CmplxSpMat &pauli_sparse_mat(Pauli p) {
  static const CmplxSpMat I_mat = const_2x2_matrix(1, 0, 0, 1);
  static const CmplxSpMat X_mat = const_2x2_matrix(0, 1, 1, 0);
  static const CmplxSpMat Y_mat = const_2x2_matrix(0, -i_, i_, 0);
  static const CmplxSpMat Z_mat = const_2x2_matrix(1, 0, 0, -1);
  switch (p) {
    case Pauli::X:
      return X_mat;
    case Pauli::Y:
      return Y_mat;
    case Pauli::Z:
      return Z_mat;
    default:
      TKET_ASSERT(p == Pauli::I);
      return I_mat;
  }
}

class StateNotPowerTwo : public std::logic_error {
 public:
  StateNotPowerTwo()
      : std::logic_error("Statevector size is not a power of two.") {}
};

unsigned get_n_qb_from_statevector(const Eigen::VectorXcd &state) {
  // allowing room for big states
  unsigned long long n = state.size();
  if (!(n && (!(n & (n - 1))))) {
    throw StateNotPowerTwo();
  }

  unsigned count = 0;
  while (n) {
    n >>= 1;
    ++count;
  }
  return count - 1;
}

CmplxSpMat QubitPauliString::to_sparse_matrix() const {
  qubit_vector_t qubits(map.size());
  unsigned i = 0;
  for (const std::pair<const Qubit, Pauli> &pair : map) {
    qubits[i] = pair.first;
    ++i;
  }
  return to_sparse_matrix(qubits);
}

CmplxSpMat QubitPauliString::to_sparse_matrix(const unsigned n_qubits) const {
  qubit_vector_t qubits(n_qubits);
  for (unsigned i = 0; i < n_qubits; ++i) {
    qubits[i] = Qubit(i);
  }
  return to_sparse_matrix(qubits);
}

CmplxSpMat QubitPauliString::to_sparse_matrix(
    const qubit_vector_t &qubits) const {
  std::vector<Pauli> paulis(qubits.size(), Pauli::I);
  std::map<Qubit, unsigned> index_map;
  unsigned index = 0;
  for (const Qubit &q : qubits) {
    index_map.insert({q, index});
    ++index;
  }
  if (index_map.size() != qubits.size())
    throw std::logic_error(
        "Qubit list given to to_sparse_matrix contains repeats");
  for (const std::pair<const Qubit, Pauli> &pair : map) {
    std::map<Qubit, unsigned>::iterator found = index_map.find(pair.first);
    if (found == index_map.end())
      throw std::logic_error(
          "Qubit list given to to_sparse_matrix doesn't contain " +
          pair.first.repr());
    paulis.at(found->second) = pair.second;
  }
  CmplxSpMat result(1, 1);
  result.insert(0, 0) = 1;
  for (Pauli p : paulis) {
    const CmplxSpMat pauli_mat = pauli_sparse_mat(p);
    result = Eigen::KroneckerProductSparse(result, pauli_mat).eval();
  }
  return result;
}

CmplxSpMat operator_tensor(
    const OperatorSum &total_operator, unsigned n_qubits) {
  qubit_vector_t qubits(n_qubits);
  for (unsigned i = 0; i < n_qubits; ++i) {
    qubits[i] = Qubit(i);
  }
  return operator_tensor(total_operator, qubits);
}

CmplxSpMat operator_tensor(
    const OperatorSum &total_operator, const qubit_vector_t &qubits) {
  CmplxSpMat sum = total_operator[0].second *
                   total_operator[0].first.to_sparse_matrix(qubits);
  for (unsigned j = 1; j < total_operator.size(); j++) {
    sum += total_operator[j].second *
           total_operator[j].first.to_sparse_matrix(qubits);
  }
  return sum;
}

Eigen::VectorXcd QubitPauliString::dot_state(
    const Eigen::VectorXcd &state) const {
  const unsigned n_qubits = get_n_qb_from_statevector(state);
  return (to_sparse_matrix(n_qubits) * state);
}

Eigen::VectorXcd QubitPauliString::dot_state(
    const Eigen::VectorXcd &state, const qubit_vector_t &qubits) const {
  if (state.size() != 1 << qubits.size())
    throw std::logic_error(
        "Size of statevector does not match number of qubits passed to "
        "dot_state");
  return (to_sparse_matrix(qubits) * state);
}

Complex QubitPauliString::state_expectation(
    const Eigen::VectorXcd &state) const {
  return state.dot(dot_state(state));
}

Complex QubitPauliString::state_expectation(
    const Eigen::VectorXcd &state, const qubit_vector_t &qubits) const {
  return state.dot(dot_state(state, qubits));
}

Complex operator_expectation(
    const OperatorSum &total_operator, const Eigen::VectorXcd &state) {
  Complex exp(0, 0);

  for (unsigned j = 0; j < total_operator.size(); j++) {
    exp += total_operator[j].second *
           total_operator[j].first.state_expectation(state);
  }

  return exp;
}

Complex operator_expectation(
    const OperatorSum &total_operator, const Eigen::VectorXcd &state,
    const qubit_vector_t &qubits) {
  Complex exp(0, 0);

  for (unsigned j = 0; j < total_operator.size(); j++) {
    exp += total_operator[j].second *
           total_operator[j].first.state_expectation(state, qubits);
  }

  return exp;
}

QubitPauliString::QubitPauliString(const std::list<Pauli> &_paulis) {
  unsigned qb_i = 0;
  for (Pauli p : _paulis) {
    map[Qubit(qb_i)] = p;
    ++qb_i;
  }
}

QubitPauliString::QubitPauliString(
    const std::list<Qubit> &qubits, const std::list<Pauli> &paulis) {
  if (qubits.size() != paulis.size()) {
    throw std::logic_error(
        "Mismatch of Qubits and Paulis upon QubitPauliString "
        "construction");
  }
  std::list<Pauli>::const_iterator p_it = paulis.begin();
  for (const Qubit &qb : qubits) {
    Pauli p = *p_it;
    if (map.find(qb) != map.end()) {
      throw std::logic_error(
          "Non-unique Qubit inserted into QubitPauliString map");
    }
    map[qb] = p;
    ++p_it;
  }
}

bool QubitPauliString::operator==(const QubitPauliString &other) const {
  return compare(other) == 0;
}

bool QubitPauliString::operator!=(const QubitPauliString &other) const {
  return !(*this == other);
}

bool QubitPauliString::operator<(const QubitPauliString &other) const {
  return compare(other) < 0;
}

int QubitPauliString::compare(const QubitPauliString &other) const {
  QubitPauliMap::const_iterator p1_it = this->map.begin();
  QubitPauliMap::const_iterator p2_it = other.map.begin();
  while (p1_it != this->map.end()) {
    if (p1_it->second == Pauli::I) {
      ++p1_it;
      continue;
    }
    while (p2_it != other.map.end() && p2_it->second == Pauli::I) {
      ++p2_it;
    }
    if (p2_it == other.map.end()) return 1;
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
  while (p2_it != other.map.end() && p2_it->second == Pauli::I) {
    ++p2_it;
  }
  return (p2_it == other.map.end()) ? 0 : -1;
}

void QubitPauliString::compress() {
  QubitPauliMap::iterator i = map.begin();
  while (i != map.end()) {
    if (i->second == Pauli::I) {
      i = map.erase(i);
    } else {
      ++i;
    }
  }
}

bool QubitPauliString::commutes_with(const QubitPauliString &other) const {
  return (conflicting_qubits(other).size() % 2) == 0;
}

std::set<Qubit> QubitPauliString::common_qubits(
    const QubitPauliString &other) const {
  std::set<Qubit> common;
  for (const std::pair<const Qubit, Pauli> &p : map) {
    QubitPauliMap::const_iterator found = other.map.find(p.first);
    if (p.second == Pauli::I) continue;
    if (found != other.map.end() && found->second == p.second) {
      common.insert(p.first);
    }
  }
  return common;
}

std::set<Qubit> QubitPauliString::own_qubits(
    const QubitPauliString &other) const {
  std::set<Qubit> own;
  for (const std::pair<const Qubit, Pauli> &p : map) {
    if (p.second == Pauli::I) continue;
    QubitPauliMap::const_iterator found = other.map.find(p.first);
    if (found == other.map.end() || found->second == Pauli::I) {
      own.insert(p.first);
    }
  }
  return own;
}

std::set<Qubit> QubitPauliString::conflicting_qubits(
    const QubitPauliString &other) const {
  std::set<Qubit> conflicts;
  for (const std::pair<const Qubit, Pauli> &p : map) {
    if (p.second == Pauli::I) continue;
    QubitPauliMap::const_iterator found = other.map.find(p.first);
    if (found != other.map.end() && found->second != Pauli::I &&
        found->second != p.second) {
      conflicts.insert(p.first);
    }
  }
  return conflicts;
}

std::string QubitPauliString::to_str() const {
  std::stringstream d;
  d << "(";
  QubitPauliMap::const_iterator i = map.begin();
  while (i != map.end()) {
    switch (i->second) {
      case Pauli::I: {
        d << "I";
        break;
      }
      case Pauli::X: {
        d << "X";
        break;
      }
      case Pauli::Y: {
        d << "Y";
        break;
      }
      case Pauli::Z: {
        d << "Z";
        break;
      }
    }
    d << i->first.repr();
    i++;
    if (i != map.end()) d << ", ";
  }
  d << ")";
  return d.str();
}

Pauli QubitPauliString::get(const Qubit &q) const {
  QubitPauliMap::const_iterator i = map.find(q);
  if (i == map.end())
    return Pauli::I;
  else
    return i->second;
}

void QubitPauliString::set(const Qubit &q, Pauli p) {
  QubitPauliMap::iterator i = map.find(q);
  if (i == map.end()) {
    if (p != Pauli::I) map.insert({q, p});
  } else {
    if (p == Pauli::I)
      map.erase(i);
    else
      i->second = p;
  }
}

std::size_t hash_value(const QubitPauliString &qps) {
  size_t seed = 0;
  for (const std::pair<const Qubit, Pauli> &qb_p : qps.map) {
    if (qb_p.second != Pauli::I) {
      boost::hash_combine(seed, qb_p.first);
      boost::hash_combine(seed, qb_p.second);
    }
  }
  return seed;
}

void to_json(nlohmann::json &j, const QubitPauliString &paulistr) {
  j = nlohmann::json::array();
  for (const auto &[qb, pauli] : paulistr.map) {
    j.push_back({qb, pauli});
  }
}

void from_json(const nlohmann::json &j, QubitPauliString &paulistr) {
  for (const auto &qb_pauli : j) {
    paulistr.set(qb_pauli[0].get<Qubit>(), qb_pauli[1].get<Pauli>());
  }
}

const QubitPauliTensor::Mult_Matrix &QubitPauliTensor::get_mult_matrix() {
  static const Mult_Matrix mult_matrix{
      {{Pauli::I, Pauli::I}, {1., Pauli::I}},
      {{Pauli::I, Pauli::X}, {1., Pauli::X}},
      {{Pauli::I, Pauli::Y}, {1., Pauli::Y}},
      {{Pauli::I, Pauli::Z}, {1., Pauli::Z}},
      {{Pauli::X, Pauli::I}, {1., Pauli::X}},
      {{Pauli::X, Pauli::X}, {1., Pauli::I}},
      {{Pauli::X, Pauli::Y}, {i_, Pauli::Z}},
      {{Pauli::X, Pauli::Z}, {-i_, Pauli::Y}},
      {{Pauli::Y, Pauli::I}, {1., Pauli::Y}},
      {{Pauli::Y, Pauli::X}, {-i_, Pauli::Z}},
      {{Pauli::Y, Pauli::Y}, {1., Pauli::I}},
      {{Pauli::Y, Pauli::Z}, {i_, Pauli::X}},
      {{Pauli::Z, Pauli::I}, {1., Pauli::Z}},
      {{Pauli::Z, Pauli::X}, {i_, Pauli::Y}},
      {{Pauli::Z, Pauli::Y}, {-i_, Pauli::X}},
      {{Pauli::Z, Pauli::Z}, {1., Pauli::I}}};
  return mult_matrix;
}

QubitPauliTensor QubitPauliTensor::operator*(
    const QubitPauliTensor &other) const {
  QubitPauliTensor result(this->coeff * other.coeff);
  QubitPauliMap::const_iterator p1i = this->string.map.begin();
  QubitPauliMap::const_iterator p2i = other.string.map.begin();
  while (p1i != this->string.map.end()) {
    while (p2i != other.string.map.end() && p2i->first < p1i->first) {
      result.string.map.insert(*p2i);
      p2i++;
    }
    if (p2i != other.string.map.end() && p2i->first == p1i->first) {
      // Pauli in the same position, so need to multiply
      const std::pair<Complex, Pauli> &prod =
          QubitPauliTensor::get_mult_matrix().at({p1i->second, p2i->second});
      result.coeff *= prod.first;
      if (prod.second != Pauli::I) {
        result.string.map.insert({p1i->first, prod.second});
      }
      p2i++;
    } else {
      result.string.map.insert(*p1i);
    }
    p1i++;
  }
  while (p2i != other.string.map.end()) {
    result.string.map.insert(*p2i);
    p2i++;
  }
  return result;
}

bool QubitPauliTensor::operator==(const QubitPauliTensor &other) const {
  if (this->coeff != other.coeff) return false;
  return (this->string == other.string);
}

bool QubitPauliTensor::operator!=(const QubitPauliTensor &other) const {
  return !(*this == other);
}

bool QubitPauliTensor::operator<(const QubitPauliTensor &other) const {
  int comp = this->string.compare(other.string);
  if (comp < 0) return true;
  if (comp > 0) return false;
  if (this->coeff.real() < other.coeff.real()) return true;
  if (this->coeff.real() > other.coeff.real()) return false;
  return (this->coeff.imag() < other.coeff.imag());
}

void QubitPauliTensor::compress() { string.compress(); }

bool QubitPauliTensor::commutes_with(const QubitPauliTensor &other) const {
  return (string.commutes_with(other.string));
}

std::set<Qubit> QubitPauliTensor::common_qubits(
    const QubitPauliTensor &other) const {
  return string.common_qubits(other.string);
}

std::set<Qubit> QubitPauliTensor::own_qubits(
    const QubitPauliTensor &other) const {
  return string.own_qubits(other.string);
}

std::set<Qubit> QubitPauliTensor::conflicting_qubits(
    const QubitPauliTensor &other) const {
  return (string.conflicting_qubits(other.string));
}

std::string QubitPauliTensor::to_str() const {
  std::stringstream d;
  if (coeff == -1.) {
    d << "-";
  } else if (coeff != 1.) {
    d << coeff << "*";
  }
  d << string.to_str();
  return d.str();
}

std::size_t hash_value(const QubitPauliTensor &qpt) {
  size_t seed = hash_value(qpt.string);
  boost::hash_combine(seed, qpt.coeff);
  return seed;
}

QubitPauliTensor operator*(Complex a, const QubitPauliTensor &qpt) {
  QubitPauliTensor result = qpt;
  result.coeff *= a;
  return result;
}

PauliStabiliser::PauliStabiliser(
    const std::vector<Pauli> string, const bool coeff)
    : string(string), coeff(coeff) {
  if (string.size() == 0) {
    throw NotValid("Pauli stabiliser cannot be empty.");
  }
  if (std::adjacent_find(string.begin(), string.end(), std::not_equal_to<>()) ==
          string.end() &&
      string[0] == Pauli::I) {
    throw NotValid("Pauli stabiliser cannot be identity.");
  }
}

bool PauliStabiliser::operator==(const PauliStabiliser &other) const {
  return coeff == other.coeff && string == other.string;
}

bool PauliStabiliser::operator!=(const PauliStabiliser &other) const {
  return coeff != other.coeff || string != other.string;
}

void to_json(nlohmann::json &j, const PauliStabiliser &pauli_stabiliser) {
  j["string"] = pauli_stabiliser.string;
  j["coeff"] = pauli_stabiliser.coeff;
}

void from_json(const nlohmann::json &j, PauliStabiliser &pauli_stabiliser) {
  pauli_stabiliser = PauliStabiliser(
      j.at("string").get<std::vector<Pauli>>(), j.at("coeff").get<bool>());
}
}  // namespace tket
