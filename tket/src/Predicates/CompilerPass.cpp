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

#include "tket/Predicates/CompilerPass.hpp"

#include <memory>
#include <optional>
#include <tklog/TketLog.hpp>

#include "tket/Mapping/RoutingMethodJson.hpp"
#include "tket/Predicates/PassGenerators.hpp"
#include "tket/Predicates/PassLibrary.hpp"
#include "tket/Transformations/ContextualReduction.hpp"
#include "tket/Transformations/PauliOptimisation.hpp"
#include "tket/Utils/Json.hpp"
#include "tket/Utils/UnitID.hpp"

namespace tket {

void trivial_callback(const CompilationUnit&, const nlohmann::json&) {}

PassConditions BasePass::get_conditions() const {
  return {precons_, postcons_};
}

std::string BasePass::to_string() const {
  std::string str = "Preconditions:\n";
  for (const TypePredicatePair& pp : precons_) {
    str += ("  " + pp.second->to_string() + "\n");
  }
  str += "Specific Postconditions:\n";
  for (const TypePredicatePair& pp : postcons_.specific_postcons_) {
    str += ("  " + pp.second->to_string() + "\n");
  }
  str += "Generic Postconditions:\n";
  for (const std::pair<const std::type_index, Guarantee>& pg :
       postcons_.generic_postcons_) {
    str += ("  " + predicate_name(pg.first) + " ");
    str += (pg.second == Guarantee::Clear) ? "Clear\n" : "Preserve\n";
  }
  str += "Default Postcondition: ";
  str += (postcons_.default_postcon_ == Guarantee::Clear) ? "Clear\n"
                                                          : "Preserve\n";
  return str;
};

std::optional<PredicatePtr> BasePass::unsatisfied_precondition(
    const CompilationUnit& c_unit, SafetyMode safe_mode) const {
  for (const TypePredicatePair& pp : precons_) {
    PredicateCache::const_iterator cache_iter = c_unit.cache_.find(pp.first);
    if (cache_iter == c_unit.cache_.end()) {  // cache does not contain
                                              // predicate
      if (!c_unit.calc_predicate(*pp.second)) return pp.second;
      c_unit.cache_.insert({pp.first, {pp.second, true}});
    } else {
      /* if a Predicate is not `true` in the cache or implied by a set Predicate
         in the cache then it is assumed to be `false` */
      if (cache_iter->second.second) {
        if (!cache_iter->second.first->implies(*pp.second)) {
          if (!c_unit.calc_predicate(*pp.second)) return pp.second;
        }
      } else {
        if (!c_unit.calc_predicate(*pp.second)) return pp.second;
      }
    }
  }
  if (safe_mode == SafetyMode::Audit) {
    for (const TypePredicatePair& pp : precons_) {
      if (!c_unit.calc_predicate(*pp.second)) return pp.second;
    }
  }
  return {};
}

void BasePass::update_cache(
    const CompilationUnit& c_unit, SafetyMode safe_mode) const {
  if (postcons_.default_postcon_ == Guarantee::Clear) {
    for (PredicateCache::iterator it = c_unit.cache_.begin();
         it != c_unit.cache_.end(); ++it) {
      it->second.second = false;
    }
  }
  for (const std::pair<const std::type_index, Guarantee>& pg :
       postcons_.generic_postcons_) {
    if (pg.second == Guarantee::Clear) {
      PredicateCache::iterator cache_iter = c_unit.cache_.find(pg.first);
      if (cache_iter != c_unit.cache_.end()) cache_iter->second.second = false;
    }
  }
  for (const TypePredicatePair& pp : postcons_.specific_postcons_) {
    if (safe_mode == SafetyMode::Audit && !pp.second->verify(c_unit.circ_))
      throw UnsatisfiedPredicate(pp.second->to_string());
    std::pair<PredicatePtr, bool> cache_pair{pp.second, true};
    c_unit.cache_[pp.first] = cache_pair;
  }
}

Guarantee BasePass::get_guarantee(const std::type_index& ti) const {
  return get_guarantee(ti, this->get_conditions());
}

Guarantee BasePass::get_guarantee(
    const std::type_index& ti, const PassConditions& conditions) {
  PredicateClassGuarantees::const_iterator class_guar_iter =
      conditions.second.generic_postcons_.find(ti);
  if (class_guar_iter == conditions.second.generic_postcons_.end()) {
    return conditions.second.default_postcon_;
  } else {
    return class_guar_iter->second;
  }
}

static PredicateClassGuarantees match_class_guarantees(
    const PassConditions& lhs, const PassConditions& rhs) {
  PredicateClassGuarantees cgs;
  for (const std::pair<const std::type_index, Guarantee>& postcon :
       rhs.second.generic_postcons_) {
    if (postcon.second == Guarantee::Preserve) {
      if (BasePass::get_guarantee(postcon.first, lhs) == Guarantee::Preserve)
        cgs.insert(postcon);
      else
        cgs.insert({postcon.first, Guarantee::Clear});
    } else
      cgs.insert(postcon);
  }
  return cgs;
}

PassConditions BasePass::match_passes(
    const PassConditions& lhs, const PassConditions& rhs, bool strict) {
  PredicatePtrMap new_precons = lhs.first;
  for (const TypePredicatePair& precon : rhs.first) {
    PredicatePtrMap::const_iterator data_guar_iter =
        lhs.second.specific_postcons_.find(precon.first);
    if (data_guar_iter == lhs.second.specific_postcons_.end()) {
      if (strict && get_guarantee(precon.first, lhs) == Guarantee::Clear) {
        throw IncompatibleCompilerPasses(precon.first);
      } else {
        PredicatePtrMap::iterator new_pre_it = new_precons.find(precon.first);
        if (new_pre_it == new_precons.end())
          new_precons.insert(precon);
        else {
          PredicatePtr to_put_in = new_pre_it->second->meet(*precon.second);
          TypePredicatePair tpp = CompilationUnit::make_type_pair(to_put_in);
          new_precons[tpp.first] = tpp.second;
        }
      }
    } else {
      if (strict && !data_guar_iter->second->implies(*precon.second)) {
        throw IncompatibleCompilerPasses(precon.first);
      }
    }
  }
  PostConditions new_postcons;
  new_postcons.specific_postcons_ = rhs.second.specific_postcons_;
  for (const TypePredicatePair& postcon : lhs.second.specific_postcons_) {
    PredicatePtrMap::iterator specific_post_it =
        new_postcons.specific_postcons_.find(postcon.first);
    if (specific_post_it == new_postcons.specific_postcons_.end()) {
      if (get_guarantee(postcon.first, rhs) == Guarantee::Preserve)
        new_postcons.specific_postcons_.insert(postcon);
    }
  }
  new_postcons.generic_postcons_ = match_class_guarantees(lhs, rhs);
  PredicateClassGuarantees others = match_class_guarantees(rhs, lhs);
  new_postcons.generic_postcons_.insert(others.begin(), others.end());

  if (rhs.second.default_postcon_ == Guarantee::Clear)
    new_postcons.default_postcon_ = Guarantee::Clear;
  else
    new_postcons.default_postcon_ = lhs.second.default_postcon_;

  return {new_precons, new_postcons};
}

PassConditions BasePass::match_passes(const PassPtr& lhs, const PassPtr& rhs) {
  return match_passes(lhs->get_conditions(), rhs->get_conditions());
}

bool StandardPass::apply(
    CompilationUnit& c_unit, SafetyMode safe_mode,
    const PassCallback& before_apply, const PassCallback& after_apply) const {
  before_apply(c_unit, this->get_config());
  std::optional<PredicatePtr> unsatisfied_precon =
      unsatisfied_precondition(c_unit, safe_mode);
  if (unsatisfied_precon)
    throw UnsatisfiedPredicate(
        unsatisfied_precon.value()
            ->to_string());  // just raise warning in super-unsafe mode
  // Allow trans_ to update the initial and final map
  bool changed = trans_.apply_fn(c_unit.circ_, c_unit.maps);
  update_cache(c_unit, safe_mode);
  after_apply(c_unit, this->get_config());
  return changed;
}

std::string StandardPass::to_string() const {
  std::string str = "***PassType: StandardPass***\n";
  str += BasePass::to_string();
  return str;
}

nlohmann::json StandardPass::get_config() const {
  nlohmann::json j;
  j["pass_class"] = "StandardPass";
  j["StandardPass"] = pass_config_;
  return j;
}

PassPtr operator>>(const PassPtr& lhs, const PassPtr& rhs) {
  PassConditions pre_post_cons = BasePass::match_passes(lhs, rhs);
  SequencePass new_pass;
  new_pass.precons_ = pre_post_cons.first;
  new_pass.postcons_ = pre_post_cons.second;
  new_pass.seq_ = {lhs, rhs};
  PassPtr sequence = std::make_shared<SequencePass>(new_pass);
  return sequence;
}

SequencePass::SequencePass(const std::vector<PassPtr>& ptvec, bool strict) {
  if (ptvec.size() == 0)
    throw std::logic_error("Cannot generate CompilerPass from empty list");
  std::vector<PassPtr>::const_iterator iter = ptvec.begin();
  PassConditions conditions = (*iter)->get_conditions();
  for (++iter; iter != ptvec.end(); ++iter) {
    const PassConditions next_cons = (*iter)->get_conditions();
    conditions = match_passes(conditions, next_cons, strict);
  }
  this->precons_ = conditions.first;
  this->postcons_ = conditions.second;
  this->seq_ = ptvec;
}

std::string SequencePass::to_string() const {
  std::string str = "***PassType: SequencePass***\n";
  str += BasePass::to_string();
  return str;
}

nlohmann::json SequencePass::get_config() const {
  nlohmann::json j;
  j["pass_class"] = "SequencePass";
  j["SequencePass"]["sequence"] = seq_;
  return j;
}

RepeatPass::RepeatPass(const PassPtr& pass, bool strict_check)
    : pass_(pass), strict_check_(strict_check) {
  /* check precons and postcons are compatible for repetition */
  std::tie(precons_, postcons_) = BasePass::match_passes(pass, pass);
}

std::string RepeatPass::to_string() const {
  std::string str = "***PassType: RepeatPass***\n";
  str += BasePass::to_string();
  return str;
}

nlohmann::json RepeatPass::get_config() const {
  nlohmann::json j;
  j["pass_class"] = "RepeatPass";
  j["RepeatPass"]["body"] = pass_;
  return j;
}

RepeatWithMetricPass::RepeatWithMetricPass(
    const PassPtr& pass, const Transform::Metric& metric)
    : pass_(pass), metric_(metric) {
  std::tie(precons_, postcons_) = BasePass::match_passes(pass, pass);
}

bool RepeatWithMetricPass::apply(
    CompilationUnit& c_unit, SafetyMode safe_mode,
    const PassCallback& before_apply, const PassCallback& after_apply) const {
  before_apply(c_unit, this->get_config());
  bool success = false;
  unsigned currentVal = metric_(c_unit.get_circ_ref());
  CompilationUnit* c_unit_current = &c_unit;
  CompilationUnit c_unit_new = c_unit;
  pass_->apply(
      c_unit_new, safe_mode);  // I can't make it apply the pass to a copy
                               // without copying the whole CompilationUnit
  unsigned newVal = metric_(c_unit_new.get_circ_ref());
  while (newVal < currentVal) {
    c_unit_current = &c_unit_new;
    currentVal = newVal;
    success = true;
    pass_->apply(c_unit_new, safe_mode, before_apply, after_apply);
    newVal = metric_(c_unit_new.get_circ_ref());
  }
  if (&c_unit != c_unit_current) c_unit = *c_unit_current;
  after_apply(c_unit, this->get_config());
  return success;
}

std::string RepeatWithMetricPass::to_string() const {
  std::string str = "***PassType: RepeatWithMetricPass***\n";
  str += BasePass::to_string();
  return str;
}

nlohmann::json RepeatWithMetricPass::get_config() const {
  nlohmann::json j;
  j["pass_class"] = "RepeatWithMetricPass";
  j["RepeatWithMetricPass"]["body"] = pass_;
  j["RepeatWithMetricPass"]["metric"] =
      "SERIALIZATION OF METRICS NOT YET IMPLEMENTED";
  return j;
}

RepeatUntilSatisfiedPass::RepeatUntilSatisfiedPass(
    const PassPtr& pass, const PredicatePtr& to_satisfy)
    : pass_(pass), pred_(to_satisfy) {
  std::tie(precons_, postcons_) = BasePass::match_passes(pass, pass);
}

bool RepeatUntilSatisfiedPass::apply(
    CompilationUnit& c_unit, SafetyMode safe_mode,
    const PassCallback& before_apply, const PassCallback& after_apply) const {
  before_apply(c_unit, this->get_config());
  bool success = false;
  while (!pred_->verify(c_unit.get_circ_ref())) {
    pass_->apply(c_unit, safe_mode, before_apply, after_apply);
    success = true;
  }
  after_apply(c_unit, this->get_config());
  return success;
}

std::string RepeatUntilSatisfiedPass::to_string() const {
  std::string str = "***PassType: RepeatUntilSatisfiedPass***\n";
  str += BasePass::to_string();
  return str;
}

nlohmann::json RepeatUntilSatisfiedPass::get_config() const {
  nlohmann::json j;
  j["pass_class"] = "RepeatUntilSatisfiedPass";
  j["RepeatUntilSatisfiedPass"]["body"] = pass_;
  j["RepeatUntilSatisfiedPass"]["predicate"] = pred_;
  return j;
}

void to_json(nlohmann::json& j, const PassPtr& pp) { j = pp->get_config(); }

}  // namespace tket
