// Copyright Quantinuum
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

#include "GreedyPauliOptimisation.hpp"

namespace tket {

namespace Transforms {

namespace GreedyPauliSimp {

struct hash_pauli_pauli {
  size_t operator()(const std::pair<Pauli, Pauli>& pair) const {
    return pair.first * 10 + pair.second;
  }
};
struct hash_optype_pauli {
  size_t operator()(const std::pair<OpType, Pauli>& pair) const {
    return static_cast<unsigned>(pair.first) * 10 + pair.second;
  }
};
// Output in range [0, 144)
constexpr size_t hash_triple(const std::tuple<TQEType, Pauli, Pauli>& t) {
  return (static_cast<size_t>(std::get<0>(t)) << 4) |
         (static_cast<size_t>(std::get<1>(t)) << 2) |
         static_cast<size_t>(std::get<2>(t));
}
struct hash_quadruple {
  size_t operator()(const std::tuple<Pauli, Pauli, Pauli, Pauli>& t) const {
    return static_cast<unsigned>(std::get<0>(t)) * 1000 +
           (std::get<1>(t) + 1) * 100 + (std::get<2>(t) + 1) * 10 +
           std::get<3>(t);
  }
};

/**
 * @brief Transform a pair of anti-commuting pauli letters at the
 * right-hand-side to Z/X For example, Sdg; H; X/Y = Z/X; Sdg; H
 */
const static std::unordered_map<
    const std::pair<Pauli, Pauli>, std::vector<OpType>, hash_pauli_pauli>
    AA_TO_ZX = {
        {{Pauli::X, Pauli::Y}, {OpType::Sdg, OpType::H}},
        {{Pauli::X, Pauli::Z}, {OpType::H}},
        {{Pauli::Y, Pauli::X}, {OpType::Vdg}},
        {{Pauli::Y, Pauli::Z}, {OpType::H, OpType::S}},
        {{Pauli::Z, Pauli::X}, {}},
        {{Pauli::Z, Pauli::Y}, {OpType::S}}};

const static std::unordered_map<OpType, OpType> SQ_CLIFF_DAGGER = {
    {OpType::H, OpType::H},
    {OpType::S, OpType::Sdg},
    {OpType::Sdg, OpType::S},
    {OpType::V, OpType::Vdg},
    {OpType::Vdg, OpType::V}};

/**
 * @brief Given a SQ Clifford gate g and a Pauli operator P, return Pauli
 * P', and sign k such that g;P; = k* P';g
 *
 */
const static std::unordered_map<
    std::pair<OpType, Pauli>, std::pair<Pauli, bool>, hash_optype_pauli>
    SQ_CLIFF_MAP = {
        {{OpType::H, Pauli::X}, {Pauli::Z, true}},
        {{OpType::S, Pauli::X}, {Pauli::Y, false}},
        {{OpType::Sdg, Pauli::X}, {Pauli::Y, true}},
        {{OpType::V, Pauli::X}, {Pauli::X, true}},
        {{OpType::Vdg, Pauli::X}, {Pauli::X, true}},
        {{OpType::X, Pauli::X}, {Pauli::X, true}},
        {{OpType::Y, Pauli::X}, {Pauli::X, false}},
        {{OpType::Z, Pauli::X}, {Pauli::X, false}},
        {{OpType::H, Pauli::Y}, {Pauli::Y, false}},
        {{OpType::S, Pauli::Y}, {Pauli::X, true}},
        {{OpType::Sdg, Pauli::Y}, {Pauli::X, false}},
        {{OpType::V, Pauli::Y}, {Pauli::Z, false}},
        {{OpType::Vdg, Pauli::Y}, {Pauli::Z, true}},
        {{OpType::X, Pauli::Y}, {Pauli::Y, false}},
        {{OpType::Y, Pauli::Y}, {Pauli::Y, true}},
        {{OpType::Z, Pauli::Y}, {Pauli::Y, false}},
        {{OpType::H, Pauli::Z}, {Pauli::X, true}},
        {{OpType::S, Pauli::Z}, {Pauli::Z, true}},
        {{OpType::Sdg, Pauli::Z}, {Pauli::Z, true}},
        {{OpType::V, Pauli::Z}, {Pauli::Y, true}},
        {{OpType::Vdg, Pauli::Z}, {Pauli::Y, false}},
        {{OpType::X, Pauli::Z}, {Pauli::Z, false}},
        {{OpType::Y, Pauli::Z}, {Pauli::Z, false}},
        {{OpType::Z, Pauli::Z}, {Pauli::Z, true}}};

/**
 * @brief Given TQE;P(0);Q(1), return P'(0), Q'(0), and sign k such that
 * TQE;P(0);Q(1) = k* P'(0);Q'(1);TQE
 */
struct TQE_PAULI_MAP {
  using Key = std::tuple<TQEType, Pauli, Pauli>;
  using Value = std::tuple<Pauli, Pauli, bool>;

 private:
  // Readable array of key-value pairs
  constexpr static std::array<std::pair<Key, Value>, 144> TQEPairs = {
      {{{TQEType::XX, Pauli::X, Pauli::X}, {Pauli::X, Pauli::X, true}},
       {{TQEType::XY, Pauli::X, Pauli::X}, {Pauli::I, Pauli::X, true}},
       {{TQEType::XZ, Pauli::X, Pauli::X}, {Pauli::I, Pauli::X, true}},
       {{TQEType::YX, Pauli::X, Pauli::X}, {Pauli::X, Pauli::I, true}},
       {{TQEType::YY, Pauli::X, Pauli::X}, {Pauli::Z, Pauli::Z, true}},
       {{TQEType::YZ, Pauli::X, Pauli::X}, {Pauli::Z, Pauli::Y, false}},
       {{TQEType::ZX, Pauli::X, Pauli::X}, {Pauli::X, Pauli::I, true}},
       {{TQEType::ZY, Pauli::X, Pauli::X}, {Pauli::Y, Pauli::Z, false}},
       {{TQEType::ZZ, Pauli::X, Pauli::X}, {Pauli::Y, Pauli::Y, true}},
       {{TQEType::XX, Pauli::X, Pauli::Y}, {Pauli::I, Pauli::Y, true}},
       {{TQEType::XY, Pauli::X, Pauli::Y}, {Pauli::X, Pauli::Y, true}},
       {{TQEType::XZ, Pauli::X, Pauli::Y}, {Pauli::I, Pauli::Y, true}},
       {{TQEType::YX, Pauli::X, Pauli::Y}, {Pauli::Z, Pauli::Z, false}},
       {{TQEType::YY, Pauli::X, Pauli::Y}, {Pauli::X, Pauli::I, true}},
       {{TQEType::YZ, Pauli::X, Pauli::Y}, {Pauli::Z, Pauli::X, true}},
       {{TQEType::ZX, Pauli::X, Pauli::Y}, {Pauli::Y, Pauli::Z, true}},
       {{TQEType::ZY, Pauli::X, Pauli::Y}, {Pauli::X, Pauli::I, true}},
       {{TQEType::ZZ, Pauli::X, Pauli::Y}, {Pauli::Y, Pauli::X, false}},
       {{TQEType::XX, Pauli::X, Pauli::Z}, {Pauli::I, Pauli::Z, true}},
       {{TQEType::XY, Pauli::X, Pauli::Z}, {Pauli::I, Pauli::Z, true}},
       {{TQEType::XZ, Pauli::X, Pauli::Z}, {Pauli::X, Pauli::Z, true}},
       {{TQEType::YX, Pauli::X, Pauli::Z}, {Pauli::Z, Pauli::Y, true}},
       {{TQEType::YY, Pauli::X, Pauli::Z}, {Pauli::Z, Pauli::X, false}},
       {{TQEType::YZ, Pauli::X, Pauli::Z}, {Pauli::X, Pauli::I, true}},
       {{TQEType::ZX, Pauli::X, Pauli::Z}, {Pauli::Y, Pauli::Y, false}},
       {{TQEType::ZY, Pauli::X, Pauli::Z}, {Pauli::Y, Pauli::X, true}},
       {{TQEType::ZZ, Pauli::X, Pauli::Z}, {Pauli::X, Pauli::I, true}},
       {{TQEType::XX, Pauli::X, Pauli::I}, {Pauli::X, Pauli::I, true}},
       {{TQEType::XY, Pauli::X, Pauli::I}, {Pauli::X, Pauli::I, true}},
       {{TQEType::XZ, Pauli::X, Pauli::I}, {Pauli::X, Pauli::I, true}},
       {{TQEType::YX, Pauli::X, Pauli::I}, {Pauli::X, Pauli::X, true}},
       {{TQEType::YY, Pauli::X, Pauli::I}, {Pauli::X, Pauli::Y, true}},
       {{TQEType::YZ, Pauli::X, Pauli::I}, {Pauli::X, Pauli::Z, true}},
       {{TQEType::ZX, Pauli::X, Pauli::I}, {Pauli::X, Pauli::X, true}},
       {{TQEType::ZY, Pauli::X, Pauli::I}, {Pauli::X, Pauli::Y, true}},
       {{TQEType::ZZ, Pauli::X, Pauli::I}, {Pauli::X, Pauli::Z, true}},
       {{TQEType::XX, Pauli::Y, Pauli::X}, {Pauli::Y, Pauli::I, true}},
       {{TQEType::XY, Pauli::Y, Pauli::X}, {Pauli::Z, Pauli::Z, false}},
       {{TQEType::XZ, Pauli::Y, Pauli::X}, {Pauli::Z, Pauli::Y, true}},
       {{TQEType::YX, Pauli::Y, Pauli::X}, {Pauli::Y, Pauli::X, true}},
       {{TQEType::YY, Pauli::Y, Pauli::X}, {Pauli::I, Pauli::X, true}},
       {{TQEType::YZ, Pauli::Y, Pauli::X}, {Pauli::I, Pauli::X, true}},
       {{TQEType::ZX, Pauli::Y, Pauli::X}, {Pauli::Y, Pauli::I, true}},
       {{TQEType::ZY, Pauli::Y, Pauli::X}, {Pauli::X, Pauli::Z, true}},
       {{TQEType::ZZ, Pauli::Y, Pauli::X}, {Pauli::X, Pauli::Y, false}},
       {{TQEType::XX, Pauli::Y, Pauli::Y}, {Pauli::Z, Pauli::Z, true}},
       {{TQEType::XY, Pauli::Y, Pauli::Y}, {Pauli::Y, Pauli::I, true}},
       {{TQEType::XZ, Pauli::Y, Pauli::Y}, {Pauli::Z, Pauli::X, false}},
       {{TQEType::YX, Pauli::Y, Pauli::Y}, {Pauli::I, Pauli::Y, true}},
       {{TQEType::YY, Pauli::Y, Pauli::Y}, {Pauli::Y, Pauli::Y, true}},
       {{TQEType::YZ, Pauli::Y, Pauli::Y}, {Pauli::I, Pauli::Y, true}},
       {{TQEType::ZX, Pauli::Y, Pauli::Y}, {Pauli::X, Pauli::Z, false}},
       {{TQEType::ZY, Pauli::Y, Pauli::Y}, {Pauli::Y, Pauli::I, true}},
       {{TQEType::ZZ, Pauli::Y, Pauli::Y}, {Pauli::X, Pauli::X, true}},
       {{TQEType::XX, Pauli::Y, Pauli::Z}, {Pauli::Z, Pauli::Y, false}},
       {{TQEType::XY, Pauli::Y, Pauli::Z}, {Pauli::Z, Pauli::X, true}},
       {{TQEType::XZ, Pauli::Y, Pauli::Z}, {Pauli::Y, Pauli::I, true}},
       {{TQEType::YX, Pauli::Y, Pauli::Z}, {Pauli::I, Pauli::Z, true}},
       {{TQEType::YY, Pauli::Y, Pauli::Z}, {Pauli::I, Pauli::Z, true}},
       {{TQEType::YZ, Pauli::Y, Pauli::Z}, {Pauli::Y, Pauli::Z, true}},
       {{TQEType::ZX, Pauli::Y, Pauli::Z}, {Pauli::X, Pauli::Y, true}},
       {{TQEType::ZY, Pauli::Y, Pauli::Z}, {Pauli::X, Pauli::X, false}},
       {{TQEType::ZZ, Pauli::Y, Pauli::Z}, {Pauli::Y, Pauli::I, true}},
       {{TQEType::XX, Pauli::Y, Pauli::I}, {Pauli::Y, Pauli::X, true}},
       {{TQEType::XY, Pauli::Y, Pauli::I}, {Pauli::Y, Pauli::Y, true}},
       {{TQEType::XZ, Pauli::Y, Pauli::I}, {Pauli::Y, Pauli::Z, true}},
       {{TQEType::YX, Pauli::Y, Pauli::I}, {Pauli::Y, Pauli::I, true}},
       {{TQEType::YY, Pauli::Y, Pauli::I}, {Pauli::Y, Pauli::I, true}},
       {{TQEType::YZ, Pauli::Y, Pauli::I}, {Pauli::Y, Pauli::I, true}},
       {{TQEType::ZX, Pauli::Y, Pauli::I}, {Pauli::Y, Pauli::X, true}},
       {{TQEType::ZY, Pauli::Y, Pauli::I}, {Pauli::Y, Pauli::Y, true}},
       {{TQEType::ZZ, Pauli::Y, Pauli::I}, {Pauli::Y, Pauli::Z, true}},
       {{TQEType::XX, Pauli::Z, Pauli::X}, {Pauli::Z, Pauli::I, true}},
       {{TQEType::XY, Pauli::Z, Pauli::X}, {Pauli::Y, Pauli::Z, true}},
       {{TQEType::XZ, Pauli::Z, Pauli::X}, {Pauli::Y, Pauli::Y, false}},
       {{TQEType::YX, Pauli::Z, Pauli::X}, {Pauli::Z, Pauli::I, true}},
       {{TQEType::YY, Pauli::Z, Pauli::X}, {Pauli::X, Pauli::Z, false}},
       {{TQEType::YZ, Pauli::Z, Pauli::X}, {Pauli::X, Pauli::Y, true}},
       {{TQEType::ZX, Pauli::Z, Pauli::X}, {Pauli::Z, Pauli::X, true}},
       {{TQEType::ZY, Pauli::Z, Pauli::X}, {Pauli::I, Pauli::X, true}},
       {{TQEType::ZZ, Pauli::Z, Pauli::X}, {Pauli::I, Pauli::X, true}},
       {{TQEType::XX, Pauli::Z, Pauli::Y}, {Pauli::Y, Pauli::Z, false}},
       {{TQEType::XY, Pauli::Z, Pauli::Y}, {Pauli::Z, Pauli::I, true}},
       {{TQEType::XZ, Pauli::Z, Pauli::Y}, {Pauli::Y, Pauli::X, true}},
       {{TQEType::YX, Pauli::Z, Pauli::Y}, {Pauli::X, Pauli::Z, true}},
       {{TQEType::YY, Pauli::Z, Pauli::Y}, {Pauli::Z, Pauli::I, true}},
       {{TQEType::YZ, Pauli::Z, Pauli::Y}, {Pauli::X, Pauli::X, false}},
       {{TQEType::ZX, Pauli::Z, Pauli::Y}, {Pauli::I, Pauli::Y, true}},
       {{TQEType::ZY, Pauli::Z, Pauli::Y}, {Pauli::Z, Pauli::Y, true}},
       {{TQEType::ZZ, Pauli::Z, Pauli::Y}, {Pauli::I, Pauli::Y, true}},
       {{TQEType::XX, Pauli::Z, Pauli::Z}, {Pauli::Y, Pauli::Y, true}},
       {{TQEType::XY, Pauli::Z, Pauli::Z}, {Pauli::Y, Pauli::X, false}},
       {{TQEType::XZ, Pauli::Z, Pauli::Z}, {Pauli::Z, Pauli::I, true}},
       {{TQEType::YX, Pauli::Z, Pauli::Z}, {Pauli::X, Pauli::Y, false}},
       {{TQEType::YY, Pauli::Z, Pauli::Z}, {Pauli::X, Pauli::X, true}},
       {{TQEType::YZ, Pauli::Z, Pauli::Z}, {Pauli::Z, Pauli::I, true}},
       {{TQEType::ZX, Pauli::Z, Pauli::Z}, {Pauli::I, Pauli::Z, true}},
       {{TQEType::ZY, Pauli::Z, Pauli::Z}, {Pauli::I, Pauli::Z, true}},
       {{TQEType::ZZ, Pauli::Z, Pauli::Z}, {Pauli::Z, Pauli::Z, true}},
       {{TQEType::XX, Pauli::Z, Pauli::I}, {Pauli::Z, Pauli::X, true}},
       {{TQEType::XY, Pauli::Z, Pauli::I}, {Pauli::Z, Pauli::Y, true}},
       {{TQEType::XZ, Pauli::Z, Pauli::I}, {Pauli::Z, Pauli::Z, true}},
       {{TQEType::YX, Pauli::Z, Pauli::I}, {Pauli::Z, Pauli::X, true}},
       {{TQEType::YY, Pauli::Z, Pauli::I}, {Pauli::Z, Pauli::Y, true}},
       {{TQEType::YZ, Pauli::Z, Pauli::I}, {Pauli::Z, Pauli::Z, true}},
       {{TQEType::ZX, Pauli::Z, Pauli::I}, {Pauli::Z, Pauli::I, true}},
       {{TQEType::ZY, Pauli::Z, Pauli::I}, {Pauli::Z, Pauli::I, true}},
       {{TQEType::ZZ, Pauli::Z, Pauli::I}, {Pauli::Z, Pauli::I, true}},
       {{TQEType::XX, Pauli::I, Pauli::X}, {Pauli::I, Pauli::X, true}},
       {{TQEType::XY, Pauli::I, Pauli::X}, {Pauli::X, Pauli::X, true}},
       {{TQEType::XZ, Pauli::I, Pauli::X}, {Pauli::X, Pauli::X, true}},
       {{TQEType::YX, Pauli::I, Pauli::X}, {Pauli::I, Pauli::X, true}},
       {{TQEType::YY, Pauli::I, Pauli::X}, {Pauli::Y, Pauli::X, true}},
       {{TQEType::YZ, Pauli::I, Pauli::X}, {Pauli::Y, Pauli::X, true}},
       {{TQEType::ZX, Pauli::I, Pauli::X}, {Pauli::I, Pauli::X, true}},
       {{TQEType::ZY, Pauli::I, Pauli::X}, {Pauli::Z, Pauli::X, true}},
       {{TQEType::ZZ, Pauli::I, Pauli::X}, {Pauli::Z, Pauli::X, true}},
       {{TQEType::XX, Pauli::I, Pauli::Y}, {Pauli::X, Pauli::Y, true}},
       {{TQEType::XY, Pauli::I, Pauli::Y}, {Pauli::I, Pauli::Y, true}},
       {{TQEType::XZ, Pauli::I, Pauli::Y}, {Pauli::X, Pauli::Y, true}},
       {{TQEType::YX, Pauli::I, Pauli::Y}, {Pauli::Y, Pauli::Y, true}},
       {{TQEType::YY, Pauli::I, Pauli::Y}, {Pauli::I, Pauli::Y, true}},
       {{TQEType::YZ, Pauli::I, Pauli::Y}, {Pauli::Y, Pauli::Y, true}},
       {{TQEType::ZX, Pauli::I, Pauli::Y}, {Pauli::Z, Pauli::Y, true}},
       {{TQEType::ZY, Pauli::I, Pauli::Y}, {Pauli::I, Pauli::Y, true}},
       {{TQEType::ZZ, Pauli::I, Pauli::Y}, {Pauli::Z, Pauli::Y, true}},
       {{TQEType::XX, Pauli::I, Pauli::Z}, {Pauli::X, Pauli::Z, true}},
       {{TQEType::XY, Pauli::I, Pauli::Z}, {Pauli::X, Pauli::Z, true}},
       {{TQEType::XZ, Pauli::I, Pauli::Z}, {Pauli::I, Pauli::Z, true}},
       {{TQEType::YX, Pauli::I, Pauli::Z}, {Pauli::Y, Pauli::Z, true}},
       {{TQEType::YY, Pauli::I, Pauli::Z}, {Pauli::Y, Pauli::Z, true}},
       {{TQEType::YZ, Pauli::I, Pauli::Z}, {Pauli::I, Pauli::Z, true}},
       {{TQEType::ZX, Pauli::I, Pauli::Z}, {Pauli::Z, Pauli::Z, true}},
       {{TQEType::ZY, Pauli::I, Pauli::Z}, {Pauli::Z, Pauli::Z, true}},
       {{TQEType::ZZ, Pauli::I, Pauli::Z}, {Pauli::I, Pauli::Z, true}},
       {{TQEType::XX, Pauli::I, Pauli::I}, {Pauli::I, Pauli::I, true}},
       {{TQEType::XY, Pauli::I, Pauli::I}, {Pauli::I, Pauli::I, true}},
       {{TQEType::XZ, Pauli::I, Pauli::I}, {Pauli::I, Pauli::I, true}},
       {{TQEType::YX, Pauli::I, Pauli::I}, {Pauli::I, Pauli::I, true}},
       {{TQEType::YY, Pauli::I, Pauli::I}, {Pauli::I, Pauli::I, true}},
       {{TQEType::YZ, Pauli::I, Pauli::I}, {Pauli::I, Pauli::I, true}},
       {{TQEType::ZX, Pauli::I, Pauli::I}, {Pauli::I, Pauli::I, true}},
       {{TQEType::ZY, Pauli::I, Pauli::I}, {Pauli::I, Pauli::I, true}},
       {{TQEType::ZZ, Pauli::I, Pauli::I}, {Pauli::I, Pauli::I, true}}}};

  // Initialize the lookup table
  constexpr static std::array<Value, 144> LookupTable = []() {
    std::array<Value, 144> lookupTable;
    for (const auto& [key, val] : TQEPairs) {
      lookupTable[hash_triple(key)] = val;
    }
    return lookupTable;
  }();

  // Pre-compute cases where the TQE gate commutes with the two Pauli operators
  constexpr static std::array<bool, 144> CommuteTable = []() {
    std::array<bool, 144> lookupTable;
    for (const auto& [key, val] : TQEPairs) {
      const auto [_tqe, p0, p1] = key;
      const auto [new_p0, new_p1, sign] = val;
      lookupTable[hash_triple(key)] = (p0 == new_p0) && (p1 == new_p1) && sign;
    }
    return lookupTable;
  }();

  // Pre-compute cost of replacing TQE;P(0);Q(1) with k*P'(0);Q'(1);TQE
  constexpr static std::array<int, 144> CostTable = []() {
    std::array<int, 144> costTable;
    for (const auto& [key, val] : TQEPairs) {
      const auto [_tqe, p0, p1] = key;
      const auto [new_p0, new_p1, _sign] = val;
      costTable[hash_triple(key)] = (p0 == Pauli::I) + (p1 == Pauli::I) -
                                    (new_p0 == Pauli::I) - (new_p1 == Pauli::I);
    }
    return costTable;
  }();

 public:
  static const Value at(const Key& key) {
    return LookupTable[hash_triple(key)];
  }

  static bool tqe_commutes(const Key& key) {
    return CommuteTable[hash_triple(key)];
  }

  static int cost_increase(const Key& key) {
    return CostTable[hash_triple(key)];
  }
};

/**
 * @brief Given non-identities P(0), Q(0),
 * return a list of TQEs, T, such that t;P(0);Q(1) = P'(0);Q'(1);t,
 * for all t in T, and one of P'(0), Q'(1) is identity.
 */
const static std::unordered_map<
    std::pair<Pauli, Pauli>, std::vector<TQEType>, hash_pauli_pauli>
    TQE_REDUCTION_MAP = {
        {{Pauli::X, Pauli::X},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::ZX}},
        {{Pauli::X, Pauli::Y},
         {TQEType::XX, TQEType::XZ, TQEType::YY, TQEType::ZY}},
        {{Pauli::X, Pauli::Z},
         {TQEType::XX, TQEType::XY, TQEType::YZ, TQEType::ZZ}},
        {{Pauli::Y, Pauli::X},
         {TQEType::XX, TQEType::YY, TQEType::YZ, TQEType::ZX}},
        {{Pauli::Y, Pauli::Y},
         {TQEType::XY, TQEType::YX, TQEType::YZ, TQEType::ZY}},
        {{Pauli::Y, Pauli::Z},
         {TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZZ}},
        {{Pauli::Z, Pauli::X},
         {TQEType::XX, TQEType::YX, TQEType::ZY, TQEType::ZZ}},
        {{Pauli::Z, Pauli::Y},
         {TQEType::XY, TQEType::YY, TQEType::ZX, TQEType::ZZ}},
        {{Pauli::Z, Pauli::Z},
         {TQEType::XZ, TQEType::YZ, TQEType::ZX, TQEType::ZY}}};

/**
 * @brief Given (P(0), P(1), Q(0), Q(1)),
 * where both pairs (P(0), Q(0)) and (P(1), Q(1)) non-trivially commute,
 * return the TQE gates that can map one pair to identities.
 */
const static std::unordered_map<
    std::tuple<Pauli, Pauli, Pauli, Pauli>, std::vector<TQEType>,
    hash_quadruple>
    CC_TO_IC_OR_CI_MAP = {
        {{Pauli::X, Pauli::X, Pauli::X, Pauli::X},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::ZX}},
        {{Pauli::X, Pauli::X, Pauli::X, Pauli::I}, {}},
        {{Pauli::X, Pauli::X, Pauli::I, Pauli::X}, {}},
        {{Pauli::X, Pauli::X, Pauli::I, Pauli::I},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::ZX}},
        {{Pauli::X, Pauli::Y, Pauli::X, Pauli::Y},
         {TQEType::XX, TQEType::XZ, TQEType::YY, TQEType::ZY}},
        {{Pauli::X, Pauli::Y, Pauli::X, Pauli::I}, {}},
        {{Pauli::X, Pauli::Y, Pauli::I, Pauli::Y}, {}},
        {{Pauli::X, Pauli::Y, Pauli::I, Pauli::I},
         {TQEType::XX, TQEType::XZ, TQEType::YY, TQEType::ZY}},
        {{Pauli::X, Pauli::Z, Pauli::X, Pauli::Z},
         {TQEType::XX, TQEType::XY, TQEType::YZ, TQEType::ZZ}},
        {{Pauli::X, Pauli::Z, Pauli::X, Pauli::I}, {}},
        {{Pauli::X, Pauli::Z, Pauli::I, Pauli::Z}, {}},
        {{Pauli::X, Pauli::Z, Pauli::I, Pauli::I},
         {TQEType::XX, TQEType::XY, TQEType::YZ, TQEType::ZZ}},
        {{Pauli::X, Pauli::I, Pauli::X, Pauli::X}, {}},
        {{Pauli::X, Pauli::I, Pauli::X, Pauli::Y}, {}},
        {{Pauli::X, Pauli::I, Pauli::X, Pauli::Z}, {}},
        {{Pauli::X, Pauli::I, Pauli::I, Pauli::X}, {}},
        {{Pauli::X, Pauli::I, Pauli::I, Pauli::Y}, {}},
        {{Pauli::X, Pauli::I, Pauli::I, Pauli::Z}, {}},
        {{Pauli::Y, Pauli::X, Pauli::Y, Pauli::X},
         {TQEType::XX, TQEType::YY, TQEType::YZ, TQEType::ZX}},
        {{Pauli::Y, Pauli::X, Pauli::Y, Pauli::I}, {}},
        {{Pauli::Y, Pauli::X, Pauli::I, Pauli::X}, {}},
        {{Pauli::Y, Pauli::X, Pauli::I, Pauli::I},
         {TQEType::XX, TQEType::YY, TQEType::YZ, TQEType::ZX}},
        {{Pauli::Y, Pauli::Y, Pauli::Y, Pauli::Y},
         {TQEType::XY, TQEType::YX, TQEType::YZ, TQEType::ZY}},
        {{Pauli::Y, Pauli::Y, Pauli::Y, Pauli::I}, {}},
        {{Pauli::Y, Pauli::Y, Pauli::I, Pauli::Y}, {}},
        {{Pauli::Y, Pauli::Y, Pauli::I, Pauli::I},
         {TQEType::XY, TQEType::YX, TQEType::YZ, TQEType::ZY}},
        {{Pauli::Y, Pauli::Z, Pauli::Y, Pauli::Z},
         {TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZZ}},
        {{Pauli::Y, Pauli::Z, Pauli::Y, Pauli::I}, {}},
        {{Pauli::Y, Pauli::Z, Pauli::I, Pauli::Z}, {}},
        {{Pauli::Y, Pauli::Z, Pauli::I, Pauli::I},
         {TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZZ}},
        {{Pauli::Y, Pauli::I, Pauli::Y, Pauli::X}, {}},
        {{Pauli::Y, Pauli::I, Pauli::Y, Pauli::Y}, {}},
        {{Pauli::Y, Pauli::I, Pauli::Y, Pauli::Z}, {}},
        {{Pauli::Y, Pauli::I, Pauli::I, Pauli::X}, {}},
        {{Pauli::Y, Pauli::I, Pauli::I, Pauli::Y}, {}},
        {{Pauli::Y, Pauli::I, Pauli::I, Pauli::Z}, {}},
        {{Pauli::Z, Pauli::X, Pauli::Z, Pauli::X},
         {TQEType::XX, TQEType::YX, TQEType::ZY, TQEType::ZZ}},
        {{Pauli::Z, Pauli::X, Pauli::Z, Pauli::I}, {}},
        {{Pauli::Z, Pauli::X, Pauli::I, Pauli::X}, {}},
        {{Pauli::Z, Pauli::X, Pauli::I, Pauli::I},
         {TQEType::XX, TQEType::YX, TQEType::ZY, TQEType::ZZ}},
        {{Pauli::Z, Pauli::Y, Pauli::Z, Pauli::Y},
         {TQEType::XY, TQEType::YY, TQEType::ZX, TQEType::ZZ}},
        {{Pauli::Z, Pauli::Y, Pauli::Z, Pauli::I}, {}},
        {{Pauli::Z, Pauli::Y, Pauli::I, Pauli::Y}, {}},
        {{Pauli::Z, Pauli::Y, Pauli::I, Pauli::I},
         {TQEType::XY, TQEType::YY, TQEType::ZX, TQEType::ZZ}},
        {{Pauli::Z, Pauli::Z, Pauli::Z, Pauli::Z},
         {TQEType::XZ, TQEType::YZ, TQEType::ZX, TQEType::ZY}},
        {{Pauli::Z, Pauli::Z, Pauli::Z, Pauli::I}, {}},
        {{Pauli::Z, Pauli::Z, Pauli::I, Pauli::Z}, {}},
        {{Pauli::Z, Pauli::Z, Pauli::I, Pauli::I},
         {TQEType::XZ, TQEType::YZ, TQEType::ZX, TQEType::ZY}},
        {{Pauli::Z, Pauli::I, Pauli::Z, Pauli::X}, {}},
        {{Pauli::Z, Pauli::I, Pauli::Z, Pauli::Y}, {}},
        {{Pauli::Z, Pauli::I, Pauli::Z, Pauli::Z}, {}},
        {{Pauli::Z, Pauli::I, Pauli::I, Pauli::X}, {}},
        {{Pauli::Z, Pauli::I, Pauli::I, Pauli::Y}, {}},
        {{Pauli::Z, Pauli::I, Pauli::I, Pauli::Z}, {}},
        {{Pauli::I, Pauli::X, Pauli::X, Pauli::X}, {}},
        {{Pauli::I, Pauli::X, Pauli::X, Pauli::I}, {}},
        {{Pauli::I, Pauli::X, Pauli::Y, Pauli::X}, {}},
        {{Pauli::I, Pauli::X, Pauli::Y, Pauli::I}, {}},
        {{Pauli::I, Pauli::X, Pauli::Z, Pauli::X}, {}},
        {{Pauli::I, Pauli::X, Pauli::Z, Pauli::I}, {}},
        {{Pauli::I, Pauli::Y, Pauli::X, Pauli::Y}, {}},
        {{Pauli::I, Pauli::Y, Pauli::X, Pauli::I}, {}},
        {{Pauli::I, Pauli::Y, Pauli::Y, Pauli::Y}, {}},
        {{Pauli::I, Pauli::Y, Pauli::Y, Pauli::I}, {}},
        {{Pauli::I, Pauli::Y, Pauli::Z, Pauli::Y}, {}},
        {{Pauli::I, Pauli::Y, Pauli::Z, Pauli::I}, {}},
        {{Pauli::I, Pauli::Z, Pauli::X, Pauli::Z}, {}},
        {{Pauli::I, Pauli::Z, Pauli::X, Pauli::I}, {}},
        {{Pauli::I, Pauli::Z, Pauli::Y, Pauli::Z}, {}},
        {{Pauli::I, Pauli::Z, Pauli::Y, Pauli::I}, {}},
        {{Pauli::I, Pauli::Z, Pauli::Z, Pauli::Z}, {}},
        {{Pauli::I, Pauli::Z, Pauli::Z, Pauli::I}, {}},
        {{Pauli::I, Pauli::I, Pauli::X, Pauli::X},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::ZX}},
        {{Pauli::I, Pauli::I, Pauli::X, Pauli::Y},
         {TQEType::XX, TQEType::XZ, TQEType::YY, TQEType::ZY}},
        {{Pauli::I, Pauli::I, Pauli::X, Pauli::Z},
         {TQEType::XX, TQEType::XY, TQEType::YZ, TQEType::ZZ}},
        {{Pauli::I, Pauli::I, Pauli::Y, Pauli::X},
         {TQEType::XX, TQEType::YY, TQEType::YZ, TQEType::ZX}},
        {{Pauli::I, Pauli::I, Pauli::Y, Pauli::Y},
         {TQEType::XY, TQEType::YX, TQEType::YZ, TQEType::ZY}},
        {{Pauli::I, Pauli::I, Pauli::Y, Pauli::Z},
         {TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZZ}},
        {{Pauli::I, Pauli::I, Pauli::Z, Pauli::X},
         {TQEType::XX, TQEType::YX, TQEType::ZY, TQEType::ZZ}},
        {{Pauli::I, Pauli::I, Pauli::Z, Pauli::Y},
         {TQEType::XY, TQEType::YY, TQEType::ZX, TQEType::ZZ}},
        {{Pauli::I, Pauli::I, Pauli::Z, Pauli::Z},
         {TQEType::XZ, TQEType::YZ, TQEType::ZX, TQEType::ZY}}};

/**
 * @brief Given (P(0), P(1), Q(0), Q(1)),
 * where both pairs (P(0), Q(0)) and (P(1), Q(1)) anti-commute,
 * return the TQE gates that can map both to non-trivial commuting pairs.
 */
const static std::unordered_map<
    std::tuple<Pauli, Pauli, Pauli, Pauli>, std::vector<TQEType>,
    hash_quadruple>
    AA_TO_CC_MAP = {
        {{Pauli::X, Pauli::X, Pauli::Y, Pauli::Y},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}},
        {{Pauli::X, Pauli::X, Pauli::Y, Pauli::Z},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::X, Pauli::X, Pauli::Z, Pauli::Y},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::X, Pauli::X, Pauli::Z, Pauli::Z},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}},
        {{Pauli::X, Pauli::Y, Pauli::Y, Pauli::X},
         {TQEType::XX, TQEType::XZ, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}},
        {{Pauli::X, Pauli::Y, Pauli::Y, Pauli::Z},
         {TQEType::XX, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::X, Pauli::Y, Pauli::Z, Pauli::X},
         {TQEType::XX, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::X, Pauli::Y, Pauli::Z, Pauli::Z},
         {TQEType::XX, TQEType::XZ, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}},
        {{Pauli::X, Pauli::Z, Pauli::Y, Pauli::X},
         {TQEType::XX, TQEType::XY, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::X, Pauli::Z, Pauli::Y, Pauli::Y},
         {TQEType::XX, TQEType::XY, TQEType::YX, TQEType::YZ, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::X, Pauli::Z, Pauli::Z, Pauli::X},
         {TQEType::XX, TQEType::XY, TQEType::YX, TQEType::YZ, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::X, Pauli::Z, Pauli::Z, Pauli::Y},
         {TQEType::XX, TQEType::XY, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::Y, Pauli::X, Pauli::X, Pauli::Y},
         {TQEType::XX, TQEType::XZ, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}},
        {{Pauli::Y, Pauli::X, Pauli::X, Pauli::Z},
         {TQEType::XX, TQEType::XY, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::Y, Pauli::X, Pauli::Z, Pauli::Y},
         {TQEType::XX, TQEType::XY, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::Y, Pauli::X, Pauli::Z, Pauli::Z},
         {TQEType::XX, TQEType::XZ, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}},
        {{Pauli::Y, Pauli::Y, Pauli::X, Pauli::X},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}},
        {{Pauli::Y, Pauli::Y, Pauli::X, Pauli::Z},
         {TQEType::XX, TQEType::XY, TQEType::YX, TQEType::YZ, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::Y, Pauli::Y, Pauli::Z, Pauli::X},
         {TQEType::XX, TQEType::XY, TQEType::YX, TQEType::YZ, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::Y, Pauli::Y, Pauli::Z, Pauli::Z},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}},
        {{Pauli::Y, Pauli::Z, Pauli::X, Pauli::X},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::Y, Pauli::Z, Pauli::X, Pauli::Y},
         {TQEType::XX, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::Y, Pauli::Z, Pauli::Z, Pauli::X},
         {TQEType::XX, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::Y, Pauli::Z, Pauli::Z, Pauli::Y},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::Z, Pauli::X, Pauli::X, Pauli::Y},
         {TQEType::XX, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::Z, Pauli::X, Pauli::X, Pauli::Z},
         {TQEType::XX, TQEType::XY, TQEType::YX, TQEType::YZ, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::Z, Pauli::X, Pauli::Y, Pauli::Y},
         {TQEType::XX, TQEType::XY, TQEType::YX, TQEType::YZ, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::Z, Pauli::X, Pauli::Y, Pauli::Z},
         {TQEType::XX, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZY,
          TQEType::ZZ}},
        {{Pauli::Z, Pauli::Y, Pauli::X, Pauli::X},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::Z, Pauli::Y, Pauli::X, Pauli::Z},
         {TQEType::XX, TQEType::XY, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::Z, Pauli::Y, Pauli::Y, Pauli::X},
         {TQEType::XX, TQEType::XY, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::Z, Pauli::Y, Pauli::Y, Pauli::Z},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YY, TQEType::ZX,
          TQEType::ZZ}},
        {{Pauli::Z, Pauli::Z, Pauli::X, Pauli::X},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}},
        {{Pauli::Z, Pauli::Z, Pauli::X, Pauli::Y},
         {TQEType::XX, TQEType::XZ, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}},
        {{Pauli::Z, Pauli::Z, Pauli::Y, Pauli::X},
         {TQEType::XX, TQEType::XZ, TQEType::YY, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}},
        {{Pauli::Z, Pauli::Z, Pauli::Y, Pauli::Y},
         {TQEType::XY, TQEType::XZ, TQEType::YX, TQEType::YZ, TQEType::ZX,
          TQEType::ZY}}};

/**
 * @brief Given (P(0), P(1), Q(0), Q(1)),
 * where P(0), Q(0) anti-commute and P(1), Q(1) non-trivially commute (not both
 * identity), return the TQE gate that maps P(1), Q(1) to identities.
 */
const static std::unordered_map<
    std::tuple<Pauli, Pauli, Pauli, Pauli>, std::vector<TQEType>,
    hash_quadruple>
    AC_TO_AI_MAP = {
        {{Pauli::X, Pauli::X, Pauli::Y, Pauli::X}, {TQEType::ZX}},
        {{Pauli::X, Pauli::X, Pauli::Y, Pauli::I}, {TQEType::YX}},
        {{Pauli::X, Pauli::X, Pauli::Z, Pauli::X}, {TQEType::YX}},
        {{Pauli::X, Pauli::X, Pauli::Z, Pauli::I}, {TQEType::ZX}},
        {{Pauli::X, Pauli::Y, Pauli::Y, Pauli::Y}, {TQEType::ZY}},
        {{Pauli::X, Pauli::Y, Pauli::Y, Pauli::I}, {TQEType::YY}},
        {{Pauli::X, Pauli::Y, Pauli::Z, Pauli::Y}, {TQEType::YY}},
        {{Pauli::X, Pauli::Y, Pauli::Z, Pauli::I}, {TQEType::ZY}},
        {{Pauli::X, Pauli::Z, Pauli::Y, Pauli::Z}, {TQEType::ZZ}},
        {{Pauli::X, Pauli::Z, Pauli::Y, Pauli::I}, {TQEType::YZ}},
        {{Pauli::X, Pauli::Z, Pauli::Z, Pauli::Z}, {TQEType::YZ}},
        {{Pauli::X, Pauli::Z, Pauli::Z, Pauli::I}, {TQEType::ZZ}},
        {{Pauli::X, Pauli::I, Pauli::Y, Pauli::X}, {TQEType::XX}},
        {{Pauli::X, Pauli::I, Pauli::Y, Pauli::Y}, {TQEType::XY}},
        {{Pauli::X, Pauli::I, Pauli::Y, Pauli::Z}, {TQEType::XZ}},
        {{Pauli::X, Pauli::I, Pauli::Z, Pauli::X}, {TQEType::XX}},
        {{Pauli::X, Pauli::I, Pauli::Z, Pauli::Y}, {TQEType::XY}},
        {{Pauli::X, Pauli::I, Pauli::Z, Pauli::Z}, {TQEType::XZ}},
        {{Pauli::Y, Pauli::X, Pauli::X, Pauli::X}, {TQEType::ZX}},
        {{Pauli::Y, Pauli::X, Pauli::X, Pauli::I}, {TQEType::XX}},
        {{Pauli::Y, Pauli::X, Pauli::Z, Pauli::X}, {TQEType::XX}},
        {{Pauli::Y, Pauli::X, Pauli::Z, Pauli::I}, {TQEType::ZX}},
        {{Pauli::Y, Pauli::Y, Pauli::X, Pauli::Y}, {TQEType::ZY}},
        {{Pauli::Y, Pauli::Y, Pauli::X, Pauli::I}, {TQEType::XY}},
        {{Pauli::Y, Pauli::Y, Pauli::Z, Pauli::Y}, {TQEType::XY}},
        {{Pauli::Y, Pauli::Y, Pauli::Z, Pauli::I}, {TQEType::ZY}},
        {{Pauli::Y, Pauli::Z, Pauli::X, Pauli::Z}, {TQEType::ZZ}},
        {{Pauli::Y, Pauli::Z, Pauli::X, Pauli::I}, {TQEType::XZ}},
        {{Pauli::Y, Pauli::Z, Pauli::Z, Pauli::Z}, {TQEType::XZ}},
        {{Pauli::Y, Pauli::Z, Pauli::Z, Pauli::I}, {TQEType::ZZ}},
        {{Pauli::Y, Pauli::I, Pauli::X, Pauli::X}, {TQEType::YX}},
        {{Pauli::Y, Pauli::I, Pauli::X, Pauli::Y}, {TQEType::YY}},
        {{Pauli::Y, Pauli::I, Pauli::X, Pauli::Z}, {TQEType::YZ}},
        {{Pauli::Y, Pauli::I, Pauli::Z, Pauli::X}, {TQEType::YX}},
        {{Pauli::Y, Pauli::I, Pauli::Z, Pauli::Y}, {TQEType::YY}},
        {{Pauli::Y, Pauli::I, Pauli::Z, Pauli::Z}, {TQEType::YZ}},
        {{Pauli::Z, Pauli::X, Pauli::X, Pauli::X}, {TQEType::YX}},
        {{Pauli::Z, Pauli::X, Pauli::X, Pauli::I}, {TQEType::XX}},
        {{Pauli::Z, Pauli::X, Pauli::Y, Pauli::X}, {TQEType::XX}},
        {{Pauli::Z, Pauli::X, Pauli::Y, Pauli::I}, {TQEType::YX}},
        {{Pauli::Z, Pauli::Y, Pauli::X, Pauli::Y}, {TQEType::YY}},
        {{Pauli::Z, Pauli::Y, Pauli::X, Pauli::I}, {TQEType::XY}},
        {{Pauli::Z, Pauli::Y, Pauli::Y, Pauli::Y}, {TQEType::XY}},
        {{Pauli::Z, Pauli::Y, Pauli::Y, Pauli::I}, {TQEType::YY}},
        {{Pauli::Z, Pauli::Z, Pauli::X, Pauli::Z}, {TQEType::YZ}},
        {{Pauli::Z, Pauli::Z, Pauli::X, Pauli::I}, {TQEType::XZ}},
        {{Pauli::Z, Pauli::Z, Pauli::Y, Pauli::Z}, {TQEType::XZ}},
        {{Pauli::Z, Pauli::Z, Pauli::Y, Pauli::I}, {TQEType::YZ}},
        {{Pauli::Z, Pauli::I, Pauli::X, Pauli::X}, {TQEType::ZX}},
        {{Pauli::Z, Pauli::I, Pauli::X, Pauli::Y}, {TQEType::ZY}},
        {{Pauli::Z, Pauli::I, Pauli::X, Pauli::Z}, {TQEType::ZZ}},
        {{Pauli::Z, Pauli::I, Pauli::Y, Pauli::X}, {TQEType::ZX}},
        {{Pauli::Z, Pauli::I, Pauli::Y, Pauli::Y}, {TQEType::ZY}},
        {{Pauli::Z, Pauli::I, Pauli::Y, Pauli::Z}, {TQEType::ZZ}}};

}  // namespace GreedyPauliSimp

}  // namespace Transforms

}  // namespace tket
