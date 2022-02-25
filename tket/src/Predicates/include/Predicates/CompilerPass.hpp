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

#include "CompilationUnit.hpp"
#include "Predicates.hpp"
#include "Utils/Json.hpp"

namespace tket {

enum class Guarantee;
struct PostConditions;
class BasePass;
class StandardPass;
class SequencePass;
class RepeatPass;
typedef std::shared_ptr<BasePass> PassPtr;
typedef std::map<std::type_index, Guarantee> PredicateClassGuarantees;
typedef std::pair<PredicatePtrMap, PostConditions> PassConditions;
typedef std::function<void(const CompilationUnit&, const nlohmann::json&)>
    PassCallback;

JSON_DECL(PassPtr)

class IncompatibleCompilerPasses : public std::logic_error {
 public:
  explicit IncompatibleCompilerPasses(const std::type_index& typeid1)
      : std::logic_error(
            "Cannot compose these Compiler Passes due to mismatching "
            "Predicates of type: " +
            predicate_name(typeid1)) {}
};

// dictate how a `CompilerPass` should affect a class of Predicates
// PredicateClassGuarantees can't `Set` Predicates, they can only Clear or
// Preserve
enum class Guarantee { Clear, Preserve };

enum class SafetyMode {
  Audit, /* Check after every pass that the cache is being updated correctly */
  Default, /* Check composition and check precons and postcons */
  Off      /* Check only composition (this should perhaps be removed as well) */
};

/* Guarantee composition for CompilerPass composition:
    Given A, B, ..., Z CompilerPasss, the Guarantee of Predicate 'p' for A >> B
   >> ... >> Z is the last non-"Preserve" Guarantee for 'p' in list. If all
   "Preserve", Guarantee 'p' := "Preserve"
 */

// Priority hierarchy: 1) Specific, 2) Generic, 3) Default
struct PostConditions {
  PredicatePtrMap specific_postcons_;
  PredicateClassGuarantees generic_postcons_;
  Guarantee default_postcon_;
  PostConditions(
      const PredicatePtrMap& specific_postcons = {},
      const PredicateClassGuarantees& generic_postcons = {},
      Guarantee default_postcon = Guarantee::Clear)
      : specific_postcons_(specific_postcons),
        generic_postcons_(generic_postcons),
        default_postcon_(default_postcon) {}
};

/**
 * @brief Default callback when applying a pass (does nothing)
 */
void trivial_callback(const CompilationUnit&, const nlohmann::json&);

/* Passes are used to generate full sequences of rewrite rules for Circuits. It
   internally stores pre and postcons which are composed together. Whenever a
   CompilationUnit is passed through a Pass it has its cache of Predicates
    updated accordingly. */

// Passes can be run in safe or unsafe mode (bool flag to dictate).
// Safe := runs Predicates at every StandardPass to make sure cache is being
// updated correctly Unsafe := When possible, updates the cache without running
// Predicates
// TODO: Super Unsafe AKA Cowabunga Mode := check nothing, allow everything

class BasePass {
 public:
  BasePass() {}

  /**
   * @brief Apply the pass and invoke callbacks
   * @param c_unit
   * @param before_apply Called at the start of the apply procedure.
   * The parameters are the CompilationUnit and a summary of the pass
   * configuration.
   * @param after_apply Called at the end of the apply procedure.
   * The parameters are the CompilationUnit and a summary of the pass
   * configuration.
   * @param safe_mode
   * @return True if pass modified the circuit, else False
   */
  virtual bool apply(
      CompilationUnit& c_unit, SafetyMode safe_mode = SafetyMode::Default,
      const PassCallback& before_apply = trivial_callback,
      const PassCallback& after_apply = trivial_callback) const = 0;

  friend PassPtr operator>>(const PassPtr& lhs, const PassPtr& rhs);

  virtual std::string to_string() const = 0;

  /**
   * @brief Get the config object
   *
   * @return json containing the name, and params for the pass.
   */
  virtual nlohmann::json get_config() const = 0;
  PassConditions get_conditions() const;
  Guarantee get_guarantee(const std::type_index& ti) const;
  static Guarantee get_guarantee(
      const std::type_index& ti, const PassConditions& conditions);

  virtual ~BasePass(){};

 protected:
  BasePass(const PredicatePtrMap& precons, const PostConditions& postcons)
      : precons_(precons), postcons_(postcons) {}
  PredicatePtrMap precons_;
  PostConditions postcons_;

  /**
   * Check whether any preconditions of the compilation unit are unsatisfied.
   *
   * Returns the first unsatisfied precondition found.
   *
   * @param c_unit compilation unit
   * @param[in] safe_mode safety mode
   *
   * @return unsatisfied precondition, if any
   */
  std::optional<PredicatePtr> unsatisfied_precondition(
      const CompilationUnit& c_unit, SafetyMode safe_mode) const;

  void update_cache(const CompilationUnit& c_unit, SafetyMode safe_mode) const;
  static PassConditions match_passes(const PassPtr& lhs, const PassPtr& rhs);
  static PassConditions match_passes(
      const PassConditions& lhs, const PassConditions& rhs);
};

/* Basic Pass that all combinators can be used on */
class StandardPass : public BasePass {
 public:
  /**
   * @brief Construct a new StandardPass object with info about the pass.
   *
   * @param precons
   * @param trans
   * @param postcons
   * @param pass_config A nlohmann::json object containing the name, and params
   * for the pass.
   */
  StandardPass(
      const PredicatePtrMap& precons, const Transform& trans,
      const PostConditions& postcons, const nlohmann::json& pass_config)
      : BasePass(precons, postcons), trans_(trans), pass_config_(pass_config) {}

  bool apply(
      CompilationUnit& c_unit, SafetyMode = SafetyMode::Default,
      const PassCallback& before_apply = trivial_callback,
      const PassCallback& after_apply = trivial_callback) const override;
  std::string to_string() const override;
  nlohmann::json get_config() const override;

 private:
  Transform trans_;
  nlohmann::json pass_config_ = "{\"name\": \"StandardPass\"}"_json;
};

/* Runs a sequence of Passes */
class SequencePass : public BasePass {
 public:
  SequencePass() {}
  explicit SequencePass(const std::vector<PassPtr>& ptvec);
  bool apply(
      CompilationUnit& c_unit, SafetyMode safe_mode = SafetyMode::Default,
      const PassCallback& before_apply = trivial_callback,
      const PassCallback& after_apply = trivial_callback) const override {
    before_apply(c_unit, this->get_config());
    bool success = false;
    for (const PassPtr& b : seq_)
      success |= b->apply(c_unit, safe_mode, before_apply, after_apply);
    after_apply(c_unit, this->get_config());
    return success;
  }
  std::string to_string() const override;
  nlohmann::json get_config() const override;
  std::vector<PassPtr> get_sequence() const { return seq_; }

  friend PassPtr operator>>(const PassPtr& lhs, const PassPtr& rhs);

 private:
  std::vector<PassPtr> seq_;
};

/* Repeats a Pass until it returns `false` */
class RepeatPass : public BasePass {
 public:
  explicit RepeatPass(const PassPtr& pass);
  bool apply(
      CompilationUnit& c_unit, SafetyMode safe_mode = SafetyMode::Default,
      const PassCallback& before_apply = trivial_callback,
      const PassCallback& after_apply = trivial_callback) const override {
    before_apply(c_unit, this->get_config());
    bool success = false;
    while (pass_->apply(c_unit, safe_mode, before_apply, after_apply))
      success = true;
    after_apply(c_unit, this->get_config());
    return success;
  }
  std::string to_string() const override;
  nlohmann::json get_config() const override;
  PassPtr get_pass() const { return pass_; }

 private:
  PassPtr pass_;
};

class RepeatWithMetricPass : public BasePass {
 public:
  RepeatWithMetricPass(const PassPtr& pass, const Transform::Metric& metric);
  bool apply(
      CompilationUnit& c_unit, SafetyMode safe_mode = SafetyMode::Default,
      const PassCallback& before_apply = trivial_callback,
      const PassCallback& after_apply = trivial_callback) const override;
  std::string to_string() const override;
  nlohmann::json get_config() const override;
  PassPtr get_pass() const { return pass_; }
  Transform::Metric get_metric() const { return metric_; }

 private:
  PassPtr pass_;
  Transform::Metric metric_;
};

class RepeatUntilSatisfiedPass : public BasePass {
 public:
  RepeatUntilSatisfiedPass(const PassPtr& pass, const PredicatePtr& to_satisfy);
  RepeatUntilSatisfiedPass(
      const PassPtr& pass, const std::function<bool(const Circuit&)>& func) {
    PredicatePtr custom_pred = std::make_shared<UserDefinedPredicate>(func);
    *this = RepeatUntilSatisfiedPass(pass, custom_pred);
  }
  /* Careful: If the predicate is never satisfied this will not terminate */
  bool apply(
      CompilationUnit& c_unit, SafetyMode safe_mode = SafetyMode::Default,
      const PassCallback& before_apply = trivial_callback,
      const PassCallback& after_apply = trivial_callback) const override;
  std::string to_string() const override;
  nlohmann::json get_config() const override;
  PassPtr get_pass() const { return pass_; }
  PredicatePtr get_predicate() const { return pred_; }

 private:
  PassPtr pass_;
  PredicatePtr pred_;
};

// TODO: Repeat with a metric, repeat until a Predicate is satisfied...

}  // namespace tket
