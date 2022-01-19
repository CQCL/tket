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

#include "GateNodesBuffer.hpp"

#include "Utils/Assert.hpp"
#include "Utils/Exceptions.hpp"

namespace tket {
namespace tket_sim {
namespace internal {

struct GateNodesBuffer::Impl {
  Eigen::MatrixXcd& matrix;
  const double abs_epsilon;
  const unsigned number_of_qubits;
  double global_phase;

  Impl(Eigen::MatrixXcd& matr, double abs_eps)
      : matrix(matr),
        abs_epsilon(abs_eps),
        number_of_qubits(get_number_of_qubits(matr.rows())),
        global_phase(0.0) {
    if (matr.cols() == 0) {
      throw NotValid("Matrix has zero cols");
    }
  }

  void push(const GateNode&);

  void add_global_phase(double ph) { global_phase += ph; }

  void flush();
};

void GateNodesBuffer::Impl::push(const GateNode& node) {
  // Later, we might add fancy optimisation here:
  // storing the gate for later use, looking for other compatible gates
  // acting on the same qubits to merge with this, etc. etc.
  node.apply_full_unitary(matrix, number_of_qubits);
}

void GateNodesBuffer::Impl::flush() {
  if (global_phase != 0.0) {
    const auto factor = std::polar(1.0, PI * global_phase);
    matrix *= factor;
    global_phase = 0.0;
  }
}

GateNodesBuffer::GateNodesBuffer(Eigen::MatrixXcd& matrix, double abs_epsilon)
    : pimpl(std::make_unique<Impl>(matrix, abs_epsilon)) {}

GateNodesBuffer::~GateNodesBuffer() {}

void GateNodesBuffer::push(const GateNode& node) { pimpl->push(node); }

void GateNodesBuffer::add_global_phase(double ph) {
  pimpl->add_global_phase(ph);
}

void GateNodesBuffer::flush() { pimpl->flush(); }

}  // namespace internal
}  // namespace tket_sim
}  // namespace tket
