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

#include <algorithm>
#include <array>
#include <catch2/catch.hpp>
#include <istream>
#include <sstream>

#include "../Graphs/RNG.hpp"
#include "../WeightSubgrMono/TestUtils/FixedArchitectures.hpp"
#include "PlacementWithWSM/FullPlacementResult.hpp"
#include "PlacementWithWSM/PatternGraphTimeSlices.hpp"
#include "PlacementWithWSM/TargetGraphData.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

using tket::graphs::tests::RNG;

typedef std::vector<std::set<VertexWSM>> Gates;

static Gates get_gates(
    unsigned number_of_vertices, unsigned approx_number_of_2_qubit_gates,
    unsigned approx_number_of_3_qubit_gates, RNG& rng) {
  Gates gates;

  const auto insert_gates = [&gates, &rng, number_of_vertices](
                                unsigned number_of_gates,
                                unsigned number_of_qubits) {
    for (unsigned ii = 0; ii < number_of_gates; ++ii) {
      gates.emplace_back();
      for (unsigned jj = 0; jj < number_of_qubits; ++jj) {
        gates.back().insert(rng.get_size_t(number_of_vertices - 1));
      }
    }
  };
  // 1 qubit gates are, of course, just ignored.
  insert_gates(approx_number_of_2_qubit_gates, 1);
  insert_gates(approx_number_of_2_qubit_gates, 2);
  insert_gates(approx_number_of_3_qubit_gates, 3);
  rng.do_shuffle(gates);
  return gates;
}

static std::string solve_problem(
    unsigned timeout_ms, const Gates& gates,
    const GraphEdgeWeights& target_initial_graph,
    std::optional<FullPlacementResult::Pass> pass_opt, unsigned iterations,
    unsigned expected_time_ms) {
  std::stringstream ss;
  ss << "Timeout " << timeout_ms;
  const PatternGraphTimeSlices slices(gates);
  const PatternGraphTimeSlices::WeightParameters p_parameters;
  const auto pattern_graph = slices.get_weights(p_parameters);
  const TargetGraphData::Parameters t_parameters;
  const TargetGraphData target_full_graph(target_initial_graph, t_parameters);
  unsigned two_qubit_gate_count = 0;
  unsigned multi_qubit_gate_count = 0;

  for (const auto& gate : gates) {
    if (gate.size() == 2) {
      ++two_qubit_gate_count;
    } else {
      if (gate.size() > 2) {
        ++multi_qubit_gate_count;
      }
    }
  }

  const auto p_vertices = get_vertices(pattern_graph);

  ss << "; " << gates.size() << " gates (" << two_qubit_gate_count
     << " two qubit gates; " << multi_qubit_gate_count << " >two qubit gates); "
     << slices.time_sliced_data.size()
     << " slices. P-graph (V=" << p_vertices.size()
     << ", E=" << pattern_graph.size()
     << "), T-graph (V=" << target_full_graph.sorted_vertices.size()
     << ", E=" << target_full_graph.final_data.size() << ": original "
     << target_initial_graph.size() << ").\n";

  FullPlacementResult::Parameters parameters;
  parameters.timeout_ms = timeout_ms;
  if (pass_opt) {
    parameters.pass_data_opt.emplace(pass_opt.value(), iterations);
  }
  parameters.max_iterations_opt = 5 * iterations;

  const FullPlacementResult full_result(
      pattern_graph, target_initial_graph, target_full_graph.final_data, gates,
      parameters);

  // Check the returned solution for validity.
  {
    std::set<VertexWSM> t_vertices_used;
    // It's already a map.
    for (const auto& entry : full_result.result.valid_assignments) {
      const auto& pv = entry.first;
      const auto& tv = entry.second;
      t_vertices_used.insert(tv);
      CHECK(std::binary_search(p_vertices.cbegin(), p_vertices.cend(), pv));

      CHECK(std::binary_search(
          target_full_graph.sorted_vertices.cbegin(),
          target_full_graph.sorted_vertices.cend(), tv));
    }
    CHECK(
        t_vertices_used.size() == full_result.result.valid_assignments.size());
  }
  if (full_result.result.valid_assignments.size() == p_vertices.size()) {
    ss << "Complete assignment; ";
  } else {
    ss << "Assigned " << full_result.result.valid_assignments.size() << "/"
       << p_vertices.size() << " vertices; ";
  }
  ss << full_result.str();
  const auto total_time =
      full_result.total_init_time_ms + full_result.total_search_time_ms;
  CHECK(total_time < timeout_ms + 10);
  if (expected_time_ms < 20) {
    CHECK(total_time < 10 + 2 * expected_time_ms);
  } else {
    CHECK(total_time <= 2 * expected_time_ms);
    CHECK(total_time >= expected_time_ms / 2);
  }
  if (pass_opt) {
    CHECK(full_result.pass == pass_opt.value());
  }
  CHECK(full_result.iterations_for_pass == iterations);
  return ss.str();
}

std::string get_stripped_string(const std::string& str) {
  std::istringstream iss(str);
  std::string fragment;
  std::stringstream ss;
  while (iss) {
    iss >> fragment;
    ss << fragment;
  }
  return ss.str();
}

SCENARIO("random gates; smaller") {
  const std::map<unsigned, GraphEdgeWeights> architectures{
      {7, FixedArchitectures::get_ibm_perth_7_qubits()},
      {16, FixedArchitectures::get_ibm_guadalupe_16_qubits()},
      {27, FixedArchitectures::get_ibm_montreal_27_qubits()},
      {65, FixedArchitectures::get_ibm_brooklyn_65_qubits()}};

  // In each vector:
  //  list the 3 arguments to get_gates;
  //  the number of target qubits;
  //  the best pass to use (0 for none: if it didn't time out,
  //            1 for INITIAL, 2 for COMPLETE_TARGET_GRAPH);
  //  the expected number of iterations (will be used as an upper bound
  //       only if "best pass to use" is not none);
  //  an estimate of the total time in ms.
  const std::vector<std::array<unsigned, 7>> encoded_problems{
      {5, 10, 3, 7, 0, 818, 2},
      {5, 20, 5, 7, 0, 1614, 4},
      // These are the best results with 1 second of computation.
      {10, 30, 5, 16, 1, 136233, 1000},
      {10, 50, 5, 27, 1, 181123, 1000},
      {10, 50, 5, 27, 1, 141871, 1000},
      {20, 200, 5, 27, 2, 286410, 1000},
      {20, 200, 5, 65, 2, 127565, 1000}};

  std::vector<std::string> calc_messages;
  RNG rng;
  // Note: the original results were for a timeout of 1 second,
  // so we set the timeout to be much higher, but ALSO set
  // the exact number of iterations.
  const unsigned timeout_ms = 10000;

  for (const auto& code : encoded_problems) {
    const auto gates = get_gates(code[0], code[1], code[2], rng);
    const auto& target_graph = architectures.at(code[3]);
    std::optional<FullPlacementResult::Pass> pass_opt;
    switch (code[4]) {
      case 0:
        break;
      case 1:
        pass_opt = FullPlacementResult::Pass::INITIAL;
        break;
      case 2:
        pass_opt = FullPlacementResult::Pass::COMPLETE_TARGET_GRAPH;
        break;
      default:
        REQUIRE(false);
    }
    const auto& iterations = code[5];
    const auto& expected_time_ms = code[6];

    calc_messages.emplace_back(solve_problem(
        timeout_ms, gates, target_graph, pass_opt, iterations,
        expected_time_ms));
  }

  // Note: whitespace is ignored in the comparison test,
  // for easier copy/paste.
  std::vector<std::string> expected_messages{
      "Timeout 10000; 23 gates (10 two qubit gates; 1 >two qubit gates); 10 "
      "slices. "
      "P-graph (V=5, E=9), T-graph (V=7, E=21: original 6). "
      "Complete assignment; assigned 5 qubits; 4 twoQ gates in place; 6 twoQ "
      "gates "
      "nearby; 27 total swap weights; 0 twoQ bad gates; 0 twoQ gates "
      "unassigned; 12 "
      "oneQ gates; 1 nQ gates; 0 nQ gates unassigned. "
      "Passes: 1; best: INITIAL; iterations: 818",

      "Timeout 10000; 45 gates (18 two qubit gates; 0 >two qubit gates); 14 "
      "slices. "
      "P-graph (V=5, E=8), T-graph (V=7, E=21: original 6). "
      "Complete assignment; assigned 5 qubits; 10 twoQ gates in place; 8 twoQ "
      "gates "
      "nearby; 32 total swap weights; 0 twoQ bad gates; 0 twoQ gates "
      "unassigned; 27 "
      "oneQ gates; 0 nQ gates; 0 nQ gates unassigned. "
      "Passes: 1; best: INITIAL; iterations: 1614",

      "Timeout 10000; 65 gates (31 two qubit gates; 3 >two qubit gates); 26 "
      "slices. "
      "P-graph (V=10, E=22), T-graph (V=16, E=90: original 16). "
      "Complete assignment; assigned 10 qubits; 15 twoQ gates in place; 16 "
      "twoQ "
      "gates nearby; 88 total swap weights; 0 twoQ bad gates; 0 twoQ gates "
      "unassigned; 31 oneQ gates; 3 nQ gates; 0 nQ gates unassigned. "
      "Passes: 1; best: INITIAL; iterations: 136233",

      "Timeout 10000; 105 gates (46 two qubit gates; 4 >two qubit gates); 22 "
      "slices. "
      "P-graph (V=10, E=29), T-graph (V=27, E=193: original 28). "
      "Complete assignment; assigned 10 qubits; 21 twoQ gates in place; 25 "
      "twoQ "
      "gates nearby; 199 total swap weights; 0 twoQ bad gates; 0 twoQ gates "
      "unassigned; 55 oneQ gates; 4 nQ gates; 0 nQ gates unassigned. "
      "Passes: 1; best: INITIAL; iterations: 181123",

      "Timeout 10000; 105 gates (49 two qubit gates; 3 >two qubit gates); 27 "
      "slices. "
      "P-graph (V=10, E=28), T-graph (V=27, E=193: original 28). "
      "Complete assignment; assigned 10 qubits; 20 twoQ gates in place; 29 "
      "twoQ "
      "gates nearby; 164 total swap weights; 0 twoQ bad gates; 0 twoQ gates "
      "unassigned; 53 oneQ gates; 3 nQ gates; 0 nQ gates unassigned. "
      "Passes: 1; best: INITIAL; iterations: 141871",

      "Timeout 10000; 405 gates (189 two qubit gates; 5 >two qubit gates); 49 "
      "slices. P-graph (V=20, E=123), T-graph (V=27, E=193: original 28). "
      "Complete assignment; assigned 20 qubits; 25 twoQ gates in place; 113 "
      "twoQ "
      "gates nearby; 1043 total swap weights; 51 twoQ bad gates; 0 twoQ gates "
      "unassigned; 211 oneQ gates; 5 nQ gates; 0 nQ gates unassigned. "
      "Passes: 1; best: COMPLETE_TARGET_GRAPH; iterations: 286410",

      "Timeout 10000; 405 gates (188 two qubit gates; 5 >two qubit gates); 51 "
      "slices. P-graph (V=20, E=127), T-graph (V=65, E=569: original 72). "
      "Complete assignment; assigned 20 qubits; 17 twoQ gates in place; 109 "
      "twoQ "
      "gates nearby; 973 total swap weights; 62 twoQ bad gates; 0 twoQ gates "
      "unassigned; 212 oneQ gates; 5 nQ gates; 0 nQ gates unassigned. "
      "Passes: 1; best: COMPLETE_TARGET_GRAPH; iterations: 127565"};

  CHECK(expected_messages.size() == calc_messages.size());
  if (expected_messages.size() < calc_messages.size()) {
    expected_messages.resize(calc_messages.size());
  }
  for (unsigned ii = 0; ii < calc_messages.size(); ++ii) {
    const auto expect_stripped = get_stripped_string(expected_messages[ii]);
    const auto calc_stripped = get_stripped_string(calc_messages[ii]);
    CHECK(expect_stripped.size() > expected_messages[ii].size() / 2);
    CHECK(calc_stripped.size() > calc_messages[ii].size() / 2);
    if (expect_stripped != calc_stripped) {
      // Let's strip all whitespace and see if they're still the same.
      CHECK(expected_messages[ii] == calc_messages[ii]);
    }
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
