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

#include <Circuit/Circuit.hpp>
#include <Simulation/CircuitSimulator.hpp>
#include <Transformations/BasicOptimisation.hpp>
#include <Transformations/OptimisationPass.hpp>
#include <iostream>
#include <tkassert/Assert.hpp>

using namespace tket;

int main() {
  {
    Circuit circ(2, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::Rz, {0.2}, {1}, {0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 1);

    circ.add_conditional_gate<unsigned>(OpType::H, {}, {0}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::H, {}, {1}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);

    Transforms::two_qubit_squash(OpType::TK2).apply(circ);
    std::cout << circ << std::endl;

    // Expected results
    Circuit circ_0(2);
    circ_0.add_op<unsigned>(OpType::H, {0});
    circ_0.add_op<unsigned>(OpType::H, {1});
    circ_0.add_op<unsigned>(OpType::CX, {0, 1});

    Circuit circ_1(2);
    circ_1.add_op<unsigned>(OpType::ZZPhase, 0.2, {0, 1});
    std::vector<Circuit> exp_circs{circ_0, circ_1};

    for (unsigned i = 0; i < (1 << circ.n_bits()); ++i) {
      Circuit condcirc(circ.all_qubits(), {});
      condcirc.add_phase(circ.get_phase());
      for (const Command& cmd : circ) {
        Op_ptr op = cmd.get_op_ptr();
        if (op->get_type() == OpType::Conditional) {
          const Conditional& cond = static_cast<const Conditional&>(*op);
          if (cond.get_width() != 1) {
            throw std::runtime_error(
                "Only conditionals with one bit are supported");
          }
          // Only support single bit controls
          unsigned b = cmd.get_args()[0].index()[0];
          // Get b-th bit
          unsigned val = (i >> b) % 2;
          if (val == cond.get_value()) {
            condcirc.add_op(cond.get_op(), cmd.get_qubits());
          }
        } else {
          condcirc.add_op(op, cmd.get_qubits());
        }
      }
      auto u = tket_sim::get_unitary(condcirc);
      auto exp_u = tket_sim::get_unitary(exp_circs[i]);
      std::cout << u << std::endl;
      std::cout << exp_u << std::endl;
      std::cout << "=============" << std::endl;
    }
  }
  {
    std::cout << std::endl << std::endl << "with barrier!" << std::endl;
    Circuit circ(2, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::Rz, {0.2}, {1}, {0}, 1);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 1);

    circ.add_barrier({0, 1});
    circ.add_conditional_gate<unsigned>(OpType::H, {}, {0}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::H, {}, {1}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);

    Transforms::two_qubit_squash(OpType::TK2).apply(circ);
    std::cout << circ << std::endl;

    // Expected results
    Circuit circ_0(2);
    circ_0.add_op<unsigned>(OpType::H, {0});
    circ_0.add_op<unsigned>(OpType::H, {1});
    circ_0.add_op<unsigned>(OpType::CX, {0, 1});

    Circuit circ_1(2);
    circ_1.add_op<unsigned>(OpType::ZZPhase, 0.2, {0, 1});
    std::vector<Circuit> exp_circs{circ_0, circ_1};

    for (unsigned i = 0; i < (1 << circ.n_bits()); ++i) {
      Circuit condcirc(circ.all_qubits(), {});
      condcirc.add_phase(circ.get_phase());
      for (const Command& cmd : circ) {
        Op_ptr op = cmd.get_op_ptr();
        if (op->get_type() == OpType::Conditional) {
          const Conditional& cond = static_cast<const Conditional&>(*op);
          if (cond.get_width() != 1) {
            throw std::runtime_error(
                "Only conditionals with one bit are supported");
          }
          // Only support single bit controls
          unsigned b = cmd.get_args()[0].index()[0];
          // Get b-th bit
          unsigned val = (i >> b) % 2;
          if (val == cond.get_value()) {
            condcirc.add_op(cond.get_op(), cmd.get_qubits());
          }
        } else {
          condcirc.add_op(op, cmd.get_qubits());
        }
      }
      auto u = tket_sim::get_unitary(condcirc);
      auto exp_u = tket_sim::get_unitary(exp_circs[i]);
      std::cout << u << std::endl;
      std::cout << exp_u << std::endl;
      std::cout << "=============" << std::endl;
    }
  }
}
