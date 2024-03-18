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

/**
 * @file
 * @brief Holding box for abstract expressions on Bits
 */

#include "tket/Circuit/Boxes.hpp"

namespace tket {

// Used as an intermediate class in between ClassicalExpBox and Box
class ClassicalExpBoxBase : public Box {
 public:
  ClassicalExpBoxBase() : Box(OpType::ClassicalExpBox) {}
  ClassicalExpBoxBase(const ClassicalExpBoxBase &other) : Box(other) {}
  ~ClassicalExpBoxBase() override {}

  /**
   * @brief Rename the units in the logic expression according to the
   * given Bit map. We can still mutate python objects even though the
   * function is const.
   *
   * @param bm
   * @return bool
   */
  virtual bool rename_units(const std::map<Bit, Bit> &bm) const {
    (void)bm;
    return false;
  }
};

/**
 * Holding box for abstract expressions on Bits
 * Templated by type T which holds expression.
 * T must implement the functions listed in
 * https://pybind11.readthedocs.io/en/stable/reference.html#_CPPv4I0E10object_api
 */

template <typename T>
class ClassicalExpBox : public ClassicalExpBoxBase {
 public:
  /**
   * Construct a ClassicalExpBox of specified shape with expression
   *
   * @param n_i number of input-only bits
   * @param n_io number of input/output bits
   * @param n_o number of output-only bits
   * @param exp expression stored
   */
  explicit ClassicalExpBox(unsigned n_i, unsigned n_io, unsigned n_o, T exp)
      : n_i_(n_i), n_io_(n_io), n_o_(n_o), exp_(exp) {
    for (unsigned i = 0; i < n_i; i++) {
      sig_.push_back(EdgeType::Boolean);
    }
    for (unsigned j = 0; j < n_io + n_o; j++) {
      sig_.push_back(EdgeType::Classical);
    }
  };

  /**
   * Copy constructor
   */
  ClassicalExpBox(const ClassicalExpBox &other)
      : ClassicalExpBoxBase(other),
        n_i_(other.n_i_),
        n_io_(other.n_io_),
        n_o_(other.n_o_),
        exp_(other.exp_),
        sig_(other.sig_){};

  ~ClassicalExpBox() override {}
  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &) const override {
    return Op_ptr(this);
  };

  SymSet free_symbols() const override { return SymSet(); }

  /**
   * Equality check between two ClassicalExpBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const ClassicalExpBox &other =
        dynamic_cast<const ClassicalExpBox &>(op_other);
    if (id_ == other.get_id()) return true;
    return content_equality(other);
  }

  op_signature_t get_signature() const override { return sig_; }

  /** Number of input-only bits. */
  unsigned get_n_i() const { return n_i_; }

  /** Number of input-output bits. */
  unsigned get_n_io() const { return n_io_; }

  /** Number of output-only bits. */
  unsigned get_n_o() const { return n_o_; }
  T get_exp() const { return exp_; }

  bool content_equality(const ClassicalExpBox &other) const {
    if (this->get_type() != other.get_type()) return false;
    const ClassicalExpBox &other_box =
        static_cast<const ClassicalExpBox &>(other);
    return (n_i_ == other_box.n_i_) && (n_io_ == other_box.n_io_) &&
           (n_o_ == other_box.n_o_) && (sig_ == other_box.sig_) &&
           (exp_.equal(other_box.exp_));
  }
  /**
   * @brief Rename the units in the logic expression according to the
   * given Bit map. Invokes the python function rename_args.
   *
   * @param bm
   * @return bool
   */
  bool rename_units(const std::map<Bit, Bit> &bm) const override {
    return (bool)exp_.attr("rename_args")(bm);
  }

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  // TODO decomposition TKET-1070
  void generate_circuit() const override {
    throw std::logic_error(
        "ClassicalExpBox cannot be decomposed to Circuit. Try the "
        "DecomposeClassicalExp compiler pass.");
  };

  ClassicalExpBox() : n_i_(0), n_io_(0), n_o_(0), exp_(), sig_() {}

 private:
  const unsigned n_i_;
  const unsigned n_io_;
  const unsigned n_o_;
  const T exp_;
  op_signature_t sig_;
};

}  // namespace tket
