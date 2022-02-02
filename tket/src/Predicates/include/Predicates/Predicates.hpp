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
#include <typeindex>

#include "Architecture/Architecture.hpp"
#include "Transformations/Transform.hpp"

namespace tket {
class Predicate;

typedef std::shared_ptr<Predicate> PredicatePtr;

JSON_DECL(PredicatePtr)

class IncorrectPredicate : public std::logic_error {
 public:
  explicit IncorrectPredicate(const std::string& exception_string)
      : std::logic_error(exception_string) {}
};

class UnsatisfiedPredicate : public std::logic_error {
 public:
  explicit UnsatisfiedPredicate(const std::string& pred_name)
      : std::logic_error(
            "Predicate requirements are not satisfied: " + pred_name) {}
};

const std::string& predicate_name(std::type_index idx);

/////////////////////
// PREDICATE CLASSES//
/////////////////////

// abstract interface class
class Predicate {
 public:
  virtual bool verify(const Circuit& circ) const = 0;

  // implication currently only works between predicates of the same subclass
  virtual bool implies(const Predicate& other) const = 0;
  virtual PredicatePtr meet(const Predicate& other) const = 0;
  virtual std::string to_string() const = 0;
  virtual ~Predicate(){};  // satisfy compiler
};

// all Predicate subclasses must inherit from `Predicate`
class GateSetPredicate : public Predicate {
 public:
  explicit GateSetPredicate(const OpTypeSet& allowed_types)
      : allowed_types_(allowed_types) {}
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;

  std::string to_string() const override;
  const OpTypeSet& get_allowed_types() const { return allowed_types_; }

 private:
  const OpTypeSet allowed_types_;
};

/**
 * Asserts that there are no conditional gates in the circuit.
 */
class NoClassicalControlPredicate : public Predicate {
 public:
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
};

class NoFastFeedforwardPredicate : public Predicate {
  // verifies Circuit has no classical bits which are written from quantum gates
  // and then read in later in the Circuit
 public:
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
};

class NoClassicalBitsPredicate : public Predicate {
  // this verifies that the Circuit uses no classical bits
  // (read or write -- so no measures and no classical controls)
 public:
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
};

class NoWireSwapsPredicate : public Predicate {
  /**
   * Verifies that you can follow the paths of each qubit/bit in the circuit
   * and finish on the same qubit/bit you started with
   */
 public:
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
};

class MaxTwoQubitGatesPredicate : public Predicate {
  // this verifies that the Circuit uses no gates with greater than 2 qubits
  // Barriers are ignored
 public:
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
};

class PlacementPredicate : public Predicate {
 public:
  explicit PlacementPredicate(const Architecture& arch)
      : nodes_(arch.nodes()) {}
  explicit PlacementPredicate(const node_set_t& nodes) : nodes_(nodes) {}
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
  const node_set_t& get_nodes() const { return nodes_; }

 private:
  const node_set_t nodes_;
};

class ConnectivityPredicate : public Predicate {
 public:
  explicit ConnectivityPredicate(const Architecture& arch) : arch_(arch) {}
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
  const Architecture& get_arch() const { return arch_; }

 private:
  const Architecture arch_;
};

class DirectednessPredicate : public Predicate {
 public:
  explicit DirectednessPredicate(const Architecture& arch) : arch_(arch) {}
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
  const Architecture& get_arch() const { return arch_; }

 private:
  const Architecture arch_;
};

class CliffordCircuitPredicate : public Predicate {
 public:
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
};

class UserDefinedPredicate : public Predicate {
 public:
  explicit UserDefinedPredicate(const std::function<bool(const Circuit&)>& func)
      : func_(func) {}
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate&) const override;
  PredicatePtr meet(const Predicate&) const override;
  std::string to_string() const override;

 private:
  const std::function<bool(const Circuit&)> func_;
};

class DefaultRegisterPredicate : public Predicate {
 public:
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
};

class MaxNQubitsPredicate : public Predicate {
 public:
  explicit MaxNQubitsPredicate(unsigned n_qubits) : n_qubits_(n_qubits) {}
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
  unsigned get_n_qubits() const { return n_qubits_; }

 private:
  const unsigned n_qubits_;
};

/**
 * Asserts that the circuit contains no \ref OpType::Barrier
 */
class NoBarriersPredicate : public Predicate {
 public:
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
};

/**
 * Asserts that any measurements occur at the end of the circuit
 */
class NoMidMeasurePredicate : public Predicate {
 public:
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
};

/**
 * Asserts that no gates in the circuit have symbolic parameters
 */
class NoSymbolsPredicate : public Predicate {
 public:
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
};

/**
 * Asserts that all NPhasedX gates act on all qubits
 * In the future, it might be useful to have a generic GlobalGatePredicate
 * for other global gates, or flag some gates as global
 */
class GlobalPhasedXPredicate : public Predicate {
 public:
  bool verify(const Circuit& circ) const override;
  bool implies(const Predicate& other) const override;
  PredicatePtr meet(const Predicate& other) const override;
  std::string to_string() const override;
};

}  // namespace tket
