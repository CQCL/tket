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

#include <cstdlib>
#include <fstream>
#include <iostream>

// tket includes
#include "Circuit/Circuit.hpp"

class FX_Circuit_Random : public ::benchmark::Fixture {
 public:
  FX_Circuit_Random() {}
  void SetUp(const ::benchmark::State&) {}
  void TearDown(const ::benchmark::State&) {}
  ~FX_Circuit_Random() {}
  // TODO: Find a way to pass file path from benchmark arguments rather than
  // hardcoded.
  tket::Circuit circuit_random = tket::Circuit(
      "./input_files/circuit_random_nb_qubits=20_nb_layers=200_example.tkc");
};

BENCHMARK_DEFINE_F(FX_Circuit_Random, BM_Circuit_Random)
(benchmark::State& state) {
  // Benchmark timing the Circuit::get_OpType_slices method
  for (auto _ : state) {
    circuit_random.get_OpType_slices(tket::OpType::CX);
  }
}

BENCHMARK_REGISTER_F(FX_Circuit_Random, BM_Circuit_Random)
    // ->DenseRange(0, 100, 100)
    ->Iterations(1000)
    ->Repetitions(4)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();
