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

/**
 * @file
 * @brief Classical operations
 */

#include "Op.hpp"
#include "Utils/Json.hpp"

namespace tket {

/**
 * A purely classical operation.
 */
class ClassicalOp : public Op {
 public:
  /**
   * Construct a ClassicalOp of specified shape
   *
   * this is a classical operation, only acting on the classical parts of the
   * circuit
   *
   * @param type operation type
   * @param n_i number of input-only bits
   * @param n_io number of input/output bits
   * @param n_o number of output-only bits
   * @param name name of operation
   */
  ClassicalOp(
      OpType type, unsigned n_i, unsigned n_io, unsigned n_o,
      const std::string &name = "");

  // Trivial overrides
  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &) const override {
    return Op_ptr();
  }
  SymSet free_symbols() const override { return {}; }
  unsigned n_qubits() const override { return 0; }

  op_signature_t get_signature() const override { return sig_; }

  nlohmann::json serialize() const override;

  static Op_ptr deserialize(const nlohmann::json &j);

  std::string get_name(bool latex = false) const override;

  /** Number of input-only bits. */
  unsigned get_n_i() const { return n_i_; }

  /** Number of input-output bits. */
  unsigned get_n_io() const { return n_io_; }

  /** Number of output-only bits. */
  unsigned get_n_o() const { return n_o_; }

  /**
   * Equality check between two ClassicalEvalOp instances
   */
  bool is_equal(const Op &other) const override;

 protected:
  const unsigned n_i_;
  const unsigned n_io_;
  const unsigned n_o_;
  const std::string name_;
  std::vector<EdgeType> sig_;
};

class ClassicalEvalOp : public ClassicalOp {
 public:
  /**
   * Construct a ClassicalEvalOp of specified shape
   *
   * this is a classical operation, only acting on the classical parts of the
   * Circuit In addition to the ClassicalOp  has the class the eval function in
   * the signature, which allows a evaluation of this op
   *
   * @param type operation type
   * @param n_i number of input-only bits
   * @param n_io number of input/output bits
   * @param n_o number of output-only bits
   * @param name name of operation
   */
  ClassicalEvalOp(
      OpType type, unsigned n_i, unsigned n_io, unsigned n_o,
      const std::string &name = "");

  /**
   * Evaluation
   *
   * @param x vector of input bits
   *
   * @return vector of outbut bits
   */
  virtual std::vector<bool> eval(const std::vector<bool> &x) const = 0;

  /**
   * Equality check between two ClassicalEvalOp instances
   */
  bool is_equal(const Op &other) const override;
};

/**
 * A general classical operation where all inputs are also outputs
 */
class ClassicalTransformOp : public ClassicalEvalOp {
 public:
  /**
   * Construct from a truth table.
   *
   * The truth table is represented by a vector of integers such that the j^th
   * bit (in little-endian order) of the (sum_i a_i 2^i)^th term is the j^th
   * output of the function applied to (a_i).
   *
   * @param n number of input/output bits
   * @param values table of binary-encoded values
   * @param name name of operation
   *
   * @pre n <= 32
   */
  ClassicalTransformOp(
      unsigned n, const std::vector<uint32_t> &values,
      const std::string &name = "ClassicalTransform");

  std::vector<bool> eval(const std::vector<bool> &x) const override;

  std::vector<uint32_t> get_values() const { return values_; }

 private:
  const std::vector<uint32_t> values_;
};

/**
 * Op containing a classical wasm function call
 */
class WASMOp : public ClassicalOp {
 public:
  /**
   * contains a wasm op that could be added to a circuit.
   * This op stores in its signatures which bits are interacting as input and
   * output with the call to which function of the wasm file
   *
   * @param _n total number bits it is interacting with
   * @param _ni_vec vector of bits for each input i32
   * @param _no_vec vector of bits for each output i32
   * @param _func_name name of the function
   * @param _wasm_uid uid of the wasm file to be called
   */
  WASMOp(
      unsigned _n, std::vector<unsigned> _ni_vec, std::vector<unsigned> _no_vec,
      const std::string &_func_name, const std::string &_wasm_uid);

  /**
   * return if the op is external
   */
  bool is_extern() const override { return true; }

  /**
   * serialize wasmop to json
   */
  nlohmann::json serialize() const override;

  /**
   * deserialize json to wasmop
   */
  static Op_ptr deserialize(const nlohmann::json &j);

  /**
   * Equality check between two WASMOp instances
   */
  bool is_equal(const Op &other) const override;

  /**
   * returns the number of classical bits the wasm op is acting on
   */
  unsigned get_n() const { return n_; }

  /**
   * returns the number of i32 the function is acting on
   */
  unsigned get_n_i32() const { return n_i32_; }

  /**
   * returns the vector of number of bit used for each of the input i32
   * variables
   */
  std::vector<unsigned> get_ni_vec() const { return ni_vec_; }

  /**
   * returns the vector of number of bit used for each of the output i32
   * variables
   */
  std::vector<unsigned> get_no_vec() const { return no_vec_; }

  /**
   * returns the name of the function the wasm op is using
   */
  std::string get_func_name() const { return func_name_; }

  /**
   * returns the uid of the wasm file the op is using, the file is stored on the
   * python layer
   */
  std::string get_wasm_uid() const { return wasm_uid_; }

 private:
  /**
   * total number of classical bits the op is interacting with
   */
  const unsigned n_;

  /**
   * total number of i32 inut and output variables
   */
  const unsigned n_i32_;

  /**
   * vector of bits for each input i32
   */
  const std::vector<unsigned> ni_vec_;

  /**
   * vector of bits for each output i32
   */
  const std::vector<unsigned> no_vec_;

  /**
   * name of the called function
   */
  const std::string func_name_;

  /**
   * uid of the wasm file the op is using
   */
  const std::string wasm_uid_;
};

/**
 * An operation to set some bits to specified values
 */
class SetBitsOp : public ClassicalEvalOp {
 public:
  /**
   * Construct from values.
   *
   * @param values values to set
   */
  explicit SetBitsOp(const std::vector<bool> &values)
      : ClassicalEvalOp(OpType::SetBits, 0, 0, values.size(), "SetBits"),
        values_(values) {}

  std::string get_name(bool latex) const override;

  std::vector<bool> get_values() const { return values_; }

  std::vector<bool> eval(const std::vector<bool> &x) const override;

 private:
  std::vector<bool> values_;
};

/**
 * An operation to copy some bit values
 *
 * @param n number of bits copied
 */
class CopyBitsOp : public ClassicalEvalOp {
 public:
  explicit CopyBitsOp(unsigned n)
      : ClassicalEvalOp(OpType::CopyBits, n, 0, n, "CopyBits") {}

  std::vector<bool> eval(const std::vector<bool> &x) const override;
};

/**
 * A classical operation with single output bit.
 *
 * There may be any number of input bits. The output bit is distinct from these.
 */
class PredicateOp : public ClassicalEvalOp {
 public:
  /**
   * Construct a PredicateOp of specified arity
   *
   * @param type type of classical operation
   * @param n number of input bits
   * @param name name of operation
   */
  PredicateOp(OpType type, unsigned n, const std::string &name = "")
      : ClassicalEvalOp(type, n, 0, 1, name) {}
};

/**
 * A predicate defined by a range of values in binary encoding
 */
class RangePredicateOp : public PredicateOp {
 public:
  /**
   * Construct from a lower and upper bound
   *
   * The lower and upper bounds are both inclusive. The output is set to 1 if
   * and only the encoded number is in the specified range.
   *
   * @param n number of inputs to predicate
   * @param a lower bound in little-endian encoding
   * @param b upper bound in little-endian encoding
   */
  RangePredicateOp(
      unsigned n, uint32_t a = 0,
      uint32_t b = std::numeric_limits<uint32_t>::max())
      : PredicateOp(OpType::RangePredicate, n, "RangePredicate"), a(a), b(b) {}

  std::string get_name(bool latex) const override;

  uint32_t upper() const { return b; }

  uint32_t lower() const { return a; }

  std::vector<bool> eval(const std::vector<bool> &x) const override;

  /**
   * Equality check between two RangePredicateOp instances
   */
  bool is_equal(const Op &other) const override;

 private:
  uint32_t a;
  uint32_t b;
};

/**
 * A predicate defined explicitly by a truth table
 */
class ExplicitPredicateOp : public PredicateOp {
 public:
  /**
   * @brief Construct from a table of values
   *
   * The truth table is represented by a vector of bool whose
   * (sum_i a_i 2^i)^th term is the predicate applied to (a_i).
   *
   * @param n number of inputs to predicate
   * @param values table of values
   * @param name name of operation
   *
   * @pre n <= 32
   */
  ExplicitPredicateOp(
      unsigned n, const std::vector<bool> &values,
      const std::string &name = "ExplicitPredicate");

  std::vector<bool> eval(const std::vector<bool> &x) const override;

  std::vector<bool> get_values() const { return values_; }

 private:
  const std::vector<bool> values_;
};

/**
 * A classical operation with one output bit which is also an input bit
 */
class ModifyingOp : public ClassicalEvalOp {
 public:
  /**
   * Construct a ModifyingOp of specified arity
   *
   * @param type type of classical operation
   * @param n number of input bits in addition to the modified bit
   * @param name name of operation
   */
  ModifyingOp(OpType type, unsigned n, const std::string &name)
      : ClassicalEvalOp(type, n, 1, 0, name) {}
};

/**
 * A modifying operation defined explicitly by a truth table
 */
class ExplicitModifierOp : public ModifyingOp {
 public:
  /**
   * @brief Construct from a table of values
   *
   * The truth table is represented by a vector of bool whose
   * (sum_i a_i 2^i)^th term is the predicate applied to (a_i), where a_m is
   * the initial value of the modified bit.
   *
   * @param n number of inputs to predicate in addition to the modified bit
   * @param values table of values
   * @param name name of operation
   *
   * @pre n <= 31
   */
  ExplicitModifierOp(
      unsigned n, const std::vector<bool> &values,
      const std::string &name = "ExplicitModifier");

  std::vector<bool> eval(const std::vector<bool> &x) const override;

  std::vector<bool> get_values() const { return values_; }

 private:
  const std::vector<bool> values_;
};

/**
 * A classical operation applied simultaneously to multiple bits.
 *
 * The order of arguments is: all arguments to first operation, then all
 * arguments to second operation, and so on.
 */
class MultiBitOp : public ClassicalEvalOp {
 public:
  MultiBitOp(std::shared_ptr<const ClassicalEvalOp> op, unsigned n);

  std::string get_name(bool latex) const override;

  std::shared_ptr<const ClassicalEvalOp> get_op() const { return op_; }

  unsigned get_n() const { return n_; }

  std::vector<bool> eval(const std::vector<bool> &x) const override;

  /**
   * Equality check between two MultiBitOp instances
   */
  bool is_equal(const Op &other) const override;

 private:
  std::shared_ptr<const ClassicalEvalOp> op_;
  unsigned n_;
};

/**
 * Classical NOT transform
 */
std::shared_ptr<ClassicalTransformOp> ClassicalX();

/**
 * Classical CNOT transform
 */
std::shared_ptr<ClassicalTransformOp> ClassicalCX();

/**
 * Unary NOT operator
 */
std::shared_ptr<ExplicitPredicateOp> NotOp();

/**
 * Binary AND operator
 */
std::shared_ptr<ExplicitPredicateOp> AndOp();

/**
 * Binary OR operator
 */
std::shared_ptr<ExplicitPredicateOp> OrOp();

/**
 * Binary XOR operator
 */
std::shared_ptr<ExplicitPredicateOp> XorOp();

/**
 * In-place AND with another input
 */
std::shared_ptr<ExplicitModifierOp> AndWithOp();

/**
 * In-place OR with another input
 */
std::shared_ptr<ExplicitModifierOp> OrWithOp();

/**
 * In-place XOR with another input
 */
std::shared_ptr<ExplicitModifierOp> XorWithOp();

}  // namespace tket
