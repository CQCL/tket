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
#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <vector>

#include "tket/Circuit/Circuit.hpp"
#include "typecast.hpp"
namespace nb = nanobind;

namespace tket {
template <typename ID>
static Circuit *add_gate_method_sequence(
    Circuit *circ, const Op_ptr &op,
    const nb::tket_custom::SequenceVec<ID> &args_seq,
    const nb::kwargs &kwargs) {
  std::vector<ID> args = args_seq;
  return add_gate_method(circ, op, args, kwargs);
}

template <typename ID>
static Circuit *add_gate_method(
    Circuit *circ, const Op_ptr &op, const std::vector<ID> &args,
    const nb::kwargs &kwargs) {
  if (op->get_desc().is_meta()) {
    throw CircuitInvalidity("Cannot add metaop to a circuit.");
  }
  if (op->get_desc().is_barrier()) {
    throw CircuitInvalidity(
        "Please use `add_barrier` to add a "
        "barrier to a circuit.");
  }
  static const std::set<std::string> allowed_kwargs = {
      "opgroup", "condition", "condition_bits", "condition_value"};
  for (const auto kwarg : kwargs) {
    const std::string kwargstr = nb::cast<std::string>(kwarg.first);
    if (!allowed_kwargs.contains(kwargstr)) {
      std::stringstream msg;
      msg << "Unsupported keyword argument '" << kwargstr << "'";
      throw CircuitInvalidity(msg.str());
    }
  }
  std::optional<std::string> opgroup;
  if (kwargs.contains("opgroup")) {
    opgroup = nb::cast<std::string>(kwargs["opgroup"]);
  }
  bool condition_given = kwargs.contains("condition");
  bool condition_bits_given = kwargs.contains("condition_bits");
  bool condition_value_given = kwargs.contains("condition_value");
  if (condition_given && condition_bits_given) {
    throw CircuitInvalidity("Both `condition` and `condition_bits` specified");
  }
  if (condition_value_given && !condition_bits_given) {
    throw CircuitInvalidity(
        "`condition_value` specified without `condition_bits`");
  }
  if (condition_given) {
    nb::module_ condition =
        nb::module_::import_("pytket.circuit.add_condition");
    nb::object add_condition = condition.attr("_add_condition");
    auto conditions = nb::cast<std::pair<Bit, bool>>(
        add_condition(circ, kwargs["condition"]));
    unit_vector_t new_args = {conditions.first};
    unsigned n_new_args = new_args.size();
    Op_ptr cond =
        std::make_shared<Conditional>(op, n_new_args, int(conditions.second));
    op_signature_t sig = op->get_signature();

    for (unsigned i = 0; i < args.size(); ++i) {
      switch (sig.at(i)) {
        case EdgeType::Quantum: {
          new_args.push_back(Qubit(args[i]));
          break;
        }
        case EdgeType::WASM: {
          new_args.push_back(WasmState(args[i]));
          break;
        }
        case EdgeType::Classical:
        case EdgeType::Boolean: {
          new_args.push_back(Bit(args[i]));
          break;
        }
        default: {
          TKET_ASSERT(!"add_gate_method found invalid edge type in signature");
        }
      }
    }
    circ->add_op(cond, new_args, opgroup);
  } else if (condition_bits_given) {
    std::vector<ID> new_args =
        nb::cast<std::vector<ID>>(kwargs["condition_bits"]);
    unsigned n_new_args = new_args.size();
    unsigned value = condition_value_given
                         ? nb::cast<unsigned>(kwargs["condition_value"])
                         : (1u << n_new_args) - 1;
    Op_ptr cond = std::make_shared<Conditional>(op, n_new_args, value);
    new_args.insert(new_args.end(), args.begin(), args.end());
    circ->add_op(cond, new_args, opgroup);
  } else
    circ->add_op(op, args, opgroup);
  return circ;
}
}  // namespace tket
