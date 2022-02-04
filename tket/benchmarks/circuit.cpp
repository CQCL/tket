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

#include <benchmark/benchmark.h>

// tket includes
#include "Circuit/Circuit.hpp"

class FX_Circuit : public ::benchmark::Fixture {
 public:
  FX_Circuit() {}
  void SetUp(const ::benchmark::State& state) {
    for (int nb_gates = 1; nb_gates < state.range(0); ++nb_gates) {
      circuit.add_op<unsigned>(tket::OpType::X, {0});
    }
  }
  void TearDown(const ::benchmark::State&) {}
  ~FX_Circuit() {}
  tket::Circuit circuit = tket::Circuit(4);
};

BENCHMARK_DEFINE_F(FX_Circuit, BM_Circuit)(benchmark::State& state) {
  // Benchmark timing the Circuit::get_OpType_slices method
  for (auto _ : state) {
    circuit.get_OpType_slices(tket::OpType::X);
  }
}

BENCHMARK_REGISTER_F(FX_Circuit, BM_Circuit)
    ->DenseRange(0, 1000, 100)
    ->Iterations(1000)
    ->Repetitions(4)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();
