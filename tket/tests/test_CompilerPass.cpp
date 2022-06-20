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
#include <catch2/catch_test_macros.hpp>

#include "Circuit/CircPool.hpp"
#include "Circuit/Circuit.hpp"
#include "Mapping/LexiLabelling.hpp"
#include "Mapping/LexiRoute.hpp"
#include "OpType/OpType.hpp"
#include "OpType/OpTypeFunctions.hpp"
#include "Placement/Placement.hpp"
#include "Predicates/CompilationUnit.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Predicates/PassLibrary.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Transformations/ContextualReduction.hpp"
#include "Transformations/PauliOptimisation.hpp"
#include "Transformations/Rebase.hpp"
#include "testutil.hpp"
namespace tket {
namespace test_CompilerPass {

SCENARIO("Run some basic Compiler Passes") {
  GIVEN(
      "A Compilation Unit that has unsatisfied Predicates at the end of "
      "the Compiler Pass") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    OpTypeSet ots = {OpType::CX};
    PredicatePtr gsp = std::make_shared<GateSetPredicate>(ots);
    PredicatePtrMap ppm{CompilationUnit::make_type_pair(gsp)};

    OpTypeSet ots2 = {OpType::CY};
    PredicatePtr gsp2 = std::make_shared<GateSetPredicate>(ots2);
    PredicatePtrMap ppm2{CompilationUnit::make_type_pair(gsp2)};

    CompilationUnit cu(circ, ppm);
    // safety mode off
    PostConditions pc{ppm2, {}, Guarantee::Preserve};
    PassPtr compass = std::make_shared<StandardPass>(
        ppm, Transforms::id, pc, nlohmann::json{});
    WHEN("Run a basic pass") { REQUIRE(!compass->apply(cu)); }
    // switch safety mode on
    PassPtr compass2 = std::make_shared<StandardPass>(
        ppm, Transforms::id, pc, nlohmann::json{});
    WHEN("Run something with an unsatisfied predicate") {
      REQUIRE_THROWS_AS(
          compass2->apply(cu, SafetyMode::Audit), UnsatisfiedPredicate);
    }
    WHEN("Compose 2 compatible Compiler Passes") {
      PostConditions pc3{ppm2, {}, Guarantee::Preserve};
      PassPtr compass3 = std::make_shared<StandardPass>(
          ppm2, Transforms::id, pc3, nlohmann::json{});
      PassPtr combination = compass >> compass3;

      // safety mode off
      REQUIRE(!combination->apply(cu, SafetyMode::Default));

      PassPtr combination2 = compass2 >> compass3;
      // safety mode on
      REQUIRE_THROWS_AS(
          combination2->apply(cu, SafetyMode::Audit), UnsatisfiedPredicate);
    }
    WHEN("Compose 2 incompatible Compiler Passes") {
      REQUIRE_THROWS_AS(compass2 >> compass, IncompatibleCompilerPasses);
    }
    WHEN("Add a class guarantee that invalidates the GateSetPredicate") {
      PredicateClassGuarantees pcg{
          {CompilationUnit::make_type_pair(gsp2).first, Guarantee::Clear}};
      PostConditions pc_clear{{}, pcg, Guarantee::Preserve};
      PassPtr compass_clear = std::make_shared<StandardPass>(
          ppm2, Transforms::id, pc_clear, nlohmann::json{});
      Circuit circ2(2);
      circ2.add_op<unsigned>(OpType::CY, {0, 1});
      CompilationUnit cu2(circ2, ppm);
      REQUIRE(!compass_clear->apply(cu2));
      REQUIRE(!cu2.check_all_predicates());
    }
  }
}

SCENARIO("Test that qubits added via add_qubit are tracked.") {
  GIVEN("Adding qubit via custom Pass.") {
    Circuit circ(2, 1);
    Qubit weird_qb("weird_q", 3);
    Qubit weird_qb2("weird_q", 5);
    Qubit weird_qb3("weird_qb", 7);
    Bit weird_cb("weird_c", 3, 1);
    circ.add_qubit(weird_qb);
    circ.add_qubit(weird_qb2);
    circ.add_bit(weird_cb);

    CompilationUnit cu(circ);

    // circuit bimaps property wont be changed, nor will compilation unit
    circ.add_qubit(weird_qb3);
    unit_bimap_t cu_initial = cu.get_initial_map_ref();

    auto it = cu_initial.left.find(weird_qb3);
    REQUIRE(it == cu_initial.left.end());

    // Instead add transform for running it
    Transform t =
        Transform([](Circuit& circ, std::shared_ptr<unit_bimaps_t> maps) {
          Qubit weird_qb4("weird_qb", 9);
          circ.add_qubit(weird_qb4);
          if (maps) {
            maps->initial.left.insert({weird_qb4, weird_qb4});
            maps->final.left.insert({weird_qb4, weird_qb4});
          }
          return true;
        });

    // convert to pass
    PredicatePtrMap s_ps;
    PostConditions postcon;
    auto pass =
        std::make_shared<StandardPass>(s_ps, t, postcon, nlohmann::json{});

    // Comparison qubit
    Qubit weird_qb4("weird_qb", 9);
    pass->apply(cu);
    cu_initial = cu.get_initial_map_ref();
    // check all maps to show weird_qb4 is mapped to self in both initial and
    // final
    it = cu_initial.left.find(weird_qb4);
    REQUIRE(it != cu_initial.left.end());
    REQUIRE(it->second == weird_qb4);
    auto cu_final = cu.get_final_map_ref();
    it = cu_final.left.find(weird_qb4);
    REQUIRE(it != cu_final.left.end());
    REQUIRE(it->second == weird_qb4);
  }
}
SCENARIO("Test making (mostly routing) passes using PassGenerators") {
  GIVEN("Correct pass for Predicate") {
    SquareGrid grid(1, 5);

    PassPtr cp_route = gen_default_mapping_pass(grid, false);
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {3, 4}});

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(grid);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);

    WHEN("Ran in unsafe mode") {
      REQUIRE(cp_route->apply(cu, SafetyMode::Default));
      REQUIRE(cu.check_all_predicates());
    }
    WHEN("Ran in safe mode") {
      REQUIRE(cp_route->apply(cu, SafetyMode::Audit));
      REQUIRE(cu.check_all_predicates());
    }
  }
  GIVEN("Incorrect pass for Predicate logs a warning") {
    SquareGrid grid(2, 3);

    PassPtr cp_route = gen_default_mapping_pass(grid, false);
    Circuit circ(6);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {0, 5}, {0, 3}, {1, 2}, {3, 4}});

    SquareGrid grid2(1, 6);
    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(grid2);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);

    WHEN("Ran in unsafe mode") {
      REQUIRE(cp_route->apply(cu));  // warning should be logged
      REQUIRE(!cu.check_all_predicates());
    }
    WHEN("Ran in safe mode") {
      REQUIRE(
          cp_route->apply(cu, SafetyMode::Audit));  // warning should be logged
      REQUIRE(!cu.check_all_predicates());
    }
  }
  GIVEN("Correct sequence of Synthesis, Routing and Rebasing") {
    Circuit circ(6);
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 5});
    circ.add_op<unsigned>(OpType::CZ, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Z, {3});
    circ.add_op<unsigned>(OpType::CY, {3, 4});

    SquareGrid grid(2, 3);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(grid);
    OpTypeSet ots = {OpType::CX, OpType::PhasedX, OpType::Rz};
    PredicatePtr gsp = std::make_shared<GateSetPredicate>(ots);

    PredicatePtrMap preds = {
        CompilationUnit::make_type_pair(routed_correctly),
        CompilationUnit::make_type_pair(gsp)};
    CompilationUnit cu(circ, preds);

    PassPtr cp_route = gen_default_mapping_pass(grid, false);

    Circuit cx(2);
    cx.add_op<unsigned>(OpType::CX, {0, 1});
    PassPtr pz_rebase = gen_rebase_pass(
        {OpType::CX, OpType::PhasedX, OpType::Rz}, cx,
        CircPool::tk1_to_PhasedXRz);
    PassPtr all_passes = SynthesiseTK() >> cp_route >> pz_rebase;

    REQUIRE(all_passes->apply(cu));
    REQUIRE(cu.check_all_predicates());
    WHEN("Ran in safe mode") {
      REQUIRE(all_passes->apply(cu, SafetyMode::Audit));
      REQUIRE(cu.check_all_predicates());
    }
    WHEN("Make incorrect sequence") {
      PassPtr bad_pass = cp_route >> SynthesiseTK();
      bad_pass->apply(cu);
      REQUIRE(!cu.check_all_predicates());
    }
  }
  GIVEN("Synthesise Passes in a row then routing") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    circ.add_op<unsigned>(OpType::CH, {0, 2});
    circ.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3});
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    circ.add_op<unsigned>(OpType::X, {4});
    OpTypeSet ots = {OpType::TK2, OpType::TK1, OpType::SWAP};
    PredicatePtr gsp = std::make_shared<GateSetPredicate>(ots);
    SquareGrid grid(2, 3);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(grid);
    PredicatePtrMap preds = {
        CompilationUnit::make_type_pair(routed_correctly),
        CompilationUnit::make_type_pair(gsp)};

    CompilationUnit cu(circ, preds);

    PlacementPtr pp = std::make_shared<GraphPlacement>(grid);
    PassPtr cp_route = gen_full_mapping_pass(
        grid, pp,
        {std::make_shared<LexiLabellingMethod>(),
         std::make_shared<LexiRouteRoutingMethod>()});

    PassPtr all_passes = SynthesiseHQS() >> SynthesiseOQC() >>
                         SynthesiseUMD() >> SynthesiseTK() >> cp_route;
    REQUIRE(all_passes->apply(cu));
    REQUIRE(cu.check_all_predicates());
  }
  GIVEN("A gen_euler_pass test with strict decomposition") {
    PassPtr squash = gen_euler_pass(OpType::Rz, OpType::Rx, true);
    Circuit circ(1);
    for (unsigned i = 0; i < 9; ++i) {
      circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
      circ.add_op<unsigned>(OpType::Rx, 0.5, {0});
    }
    CompilationUnit cu(circ);
    squash->apply(cu);
    const Circuit& c = cu.get_circ_ref();
    c.assert_valid();
    REQUIRE(c.n_gates() == 3);
  }
  GIVEN("A gen_euler_pass test") {
    PassPtr squash = gen_euler_pass(OpType::Rz, OpType::Rx, false);
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    for (unsigned i = 0; i < 9; ++i) {
      circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
      circ.add_op<unsigned>(OpType::Rx, 0.5, {0});
      circ.add_op<unsigned>(OpType::Rx, 0.3, {1});
      circ.add_op<unsigned>(OpType::Rz, 0.5, {1});
    }
    CompilationUnit cu(circ);
    squash->apply(cu);
    const Circuit& c = cu.get_circ_ref();
    c.assert_valid();
    REQUIRE(c.n_gates() == 3 + 3 + 1);
    auto cmds = c.get_commands();
    std::vector<OpType> expected_optypes{
        OpType::Rz, OpType::Rx,              // before CX
        OpType::CX, OpType::Rx, OpType::Rz,  // qubit 0 after CX
        OpType::Rz, OpType::Rx               // qubit 1 after CX
    };
    for (unsigned i = 0; i < expected_optypes.size(); ++i) {
      REQUIRE(cmds[i].get_op_ptr()->get_type() == expected_optypes[i]);
    }
  }
  GIVEN("A gen_euler_pass test with 2CX") {
    PassPtr squash = gen_euler_pass(OpType::Rz, OpType::Rx, false);
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    for (unsigned i = 0; i < 9; ++i) {
      circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
      circ.add_op<unsigned>(OpType::Rx, 0.5, {0});
      circ.add_op<unsigned>(OpType::Rx, 0.3, {1});
      circ.add_op<unsigned>(OpType::Rz, 0.5, {1});
    }
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    for (unsigned i = 0; i < 9; ++i) {
      circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
      circ.add_op<unsigned>(OpType::Rx, 0.5, {0});
      circ.add_op<unsigned>(OpType::Rx, 0.3, {1});
      circ.add_op<unsigned>(OpType::Rz, 0.5, {1});
    }
    CompilationUnit cu(circ);
    squash->apply(cu);
    const Circuit& c = cu.get_circ_ref();
    c.assert_valid();
    REQUIRE(c.n_gates() == 4 * 2 + 1 + 1 + 2);
    auto cmds = c.get_commands();
    std::vector<OpType> expected_optypes{
        OpType::Rx,                          // qubit 0 before CXs
        OpType::Rz,                          // qubit 1 before CXs
        OpType::CX, OpType::Rz, OpType::Rx,  // qubit 0 between CXs
        OpType::Rx, OpType::Rz,              // qubit 1 between CXs
        OpType::CX, OpType::Rx, OpType::Rz,  // qubit 0 after CXs
        OpType::Rz, OpType::Rx,              // qubit 1 after CXs
    };
    for (unsigned i = 0; i < expected_optypes.size(); ++i) {
      REQUIRE(cmds[i].get_op_ptr()->get_type() == expected_optypes[i]);
    }
  }

  GIVEN("gen_euler_pass commuting conditionals through CX") {
    PassPtr squash = gen_euler_pass(OpType::Rz, OpType::Rx);
    Circuit circ(2, 1);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_conditional_gate<unsigned>(OpType::Rz, {0.142}, {0}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::Rz, {0.143}, {0}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::Rx, {0.528}, {1}, {0}, 0);
    CompilationUnit cu(circ);
    squash->apply(cu);
    const Circuit& c = cu.get_circ_ref();
    c.assert_valid();
    REQUIRE(c.n_gates() == 3);
    std::vector<OpType> expected_optypes{
        OpType::Conditional,  // qubit 0 before CX
        OpType::Conditional,  // qubit 1 before CX
        OpType::CX};
    check_command_types(c, expected_optypes);

    auto cmds = c.get_commands();
    Op_ptr op0 =
        static_cast<const Conditional&>(*cmds[0].get_op_ptr()).get_op();
    Op_ptr op1 =
        static_cast<const Conditional&>(*cmds[1].get_op_ptr()).get_op();

    REQUIRE(op0->get_type() == OpType::Rz);
    REQUIRE(op0->get_params() == std::vector<Expr>{0.285});
    REQUIRE(op1->get_type() == OpType::Rx);
    REQUIRE(op1->get_params() == std::vector<Expr>{0.528});
  }

  GIVEN("Repeat synthesis passes") {
    OpTypeSet ots = {OpType::H};
    PredicatePtr gsp = std::make_shared<GateSetPredicate>(ots);
    PredicatePtrMap ppm{CompilationUnit::make_type_pair(gsp)};
    PostConditions pc{ppm, {}, Guarantee::Preserve};
    PassPtr compass = std::make_shared<StandardPass>(
        ppm, Transforms::id, pc, nlohmann::json{});
    PassPtr rep = std::make_shared<RepeatPass>(compass);
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::H, {0});
    CompilationUnit cu(circ);
    REQUIRE_NOTHROW(rep->apply(cu));
    cu.get_circ_ref().assert_valid();
  }
  GIVEN("Full compilation sequence") {
    SquareGrid grid(1, 5);
    std::vector<PassPtr> passes = {
        DecomposeBoxes(), RebaseTket(), gen_default_mapping_pass(grid, true)};
    REQUIRE_NOTHROW(SequencePass(passes));
  }
}

SCENARIO("Construct sequence pass") {
  std::vector<PassPtr> passes = {CommuteThroughMultis(), KAKDecomposition()};
  PassPtr sequence = std::make_shared<SequencePass>(passes);

  WHEN("Apply to valid CompilationUnit") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    CompilationUnit cu(circ);
    REQUIRE_NOTHROW(sequence->apply(cu));
  }
  WHEN("Apply to invalid CompilationUnit") {
    Circuit circ(2, 1);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
    CompilationUnit cu(circ);
    REQUIRE_THROWS_AS(sequence->apply(cu), UnsatisfiedPredicate);
  }
}

SCENARIO("Construct invalid sequence passes from vector") {
  std::vector<PassPtr> invalid_pass_to_combo{
      SynthesiseHQS(), SynthesiseOQC(), SynthesiseUMD(), SynthesiseTK()};
  for (const PassPtr& pass : invalid_pass_to_combo) {
    std::vector<PassPtr> passes = {pass};
    OpTypeSet ots = {OpType::CX};
    PredicatePtr gsp = std::make_shared<GateSetPredicate>(ots);
    PredicatePtrMap ppm{CompilationUnit::make_type_pair(gsp)};
    PostConditions pc{{}, {}, Guarantee::Preserve};
    PassPtr compass = std::make_shared<StandardPass>(
        ppm, Transforms::id, pc, nlohmann::json{});
    passes.push_back(compass);
    REQUIRE_THROWS_AS((void)SequencePass(passes), IncompatibleCompilerPasses);
  }
}

SCENARIO("Construct invalid sequence of loops") {
  PredicatePtr pp1 = std::make_shared<NoClassicalControlPredicate>();
  PredicatePtrMap ppm{CompilationUnit::make_type_pair(pp1)};
  PostConditions pc{{}, {}, Guarantee::Preserve};
  PassPtr pass1 =
      std::make_shared<StandardPass>(ppm, Transforms::id, pc, nlohmann::json{});
  PassPtr loop1 = std::make_shared<RepeatPass>(pass1);
  PostConditions pc2{{}, {}, Guarantee::Clear};
  PredicatePtrMap empty_ppm{};
  PassPtr pass2 = std::make_shared<StandardPass>(
      empty_ppm, Transforms::id, pc2, nlohmann::json{});
  PassPtr loop2 = std::make_shared<RepeatPass>(pass2);
  std::vector<PassPtr> good_passes{loop1, loop2};
  std::vector<PassPtr> bad_passes{loop2, loop1};
  REQUIRE_NOTHROW((void)SequencePass(good_passes));
  REQUIRE_THROWS_AS((void)SequencePass(bad_passes), IncompatibleCompilerPasses);
}

SCENARIO("Test RepeatWithMetricPass") {
  GIVEN("Monotonically decreasing pass") {
    PassPtr seq_p = RemoveRedundancies() >> CommuteThroughMultis();
    Transform::Metric met = [](const Circuit& circ) {
      return circ.n_vertices();
    };
    PassPtr rwm_p = std::make_shared<RepeatWithMetricPass>(seq_p, met);
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    CompilationUnit cu(circ);
    rwm_p->apply(cu);
    REQUIRE(cu.get_circ_ref().n_gates() == 1);
  }
}

SCENARIO("Track initial and final maps throughout compilation") {
  GIVEN("SynthesiseTK should not affect them") {
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CY, {{0, 3}, {1, 4}, {1, 0}, {2, 1}});
    circ.add_op<unsigned>(OpType::SWAP, {3, 4});
    circ.add_op<unsigned>(OpType::Z, {4});
    circ.replace_SWAPs();
    CompilationUnit cu(circ);
    SynthesiseTK()->apply(cu);
    for (auto pair : cu.get_initial_map_ref().left) {
      REQUIRE(pair.first == pair.second);
    }
    for (auto pair : cu.get_final_map_ref().left) {
      REQUIRE(pair.first == pair.second);
    }
  }
  GIVEN("Routing should modify them") {
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CY, {{0, 3}, {1, 4}, {1, 0}, {2, 1}});
    circ.add_op<unsigned>(OpType::SWAP, {3, 4});
    circ.add_op<unsigned>(OpType::Z, {4});
    circ.replace_SWAPs();
    unit_map_t rename_map = {
        {Qubit(0), Qubit("qa")},
        {Qubit(1), Qubit("qb")},
        {Qubit(2), Qubit("qc")},
        {Qubit(3), Qubit("qd")},
        {Qubit(4), Qubit("qe")}};
    circ.rename_units(rename_map);
    CompilationUnit cu(circ);

    SquareGrid grid(2, 3);
    PassPtr cp_route = gen_default_mapping_pass(grid, false);
    cp_route->apply(cu);
    bool ids_updated = true;
    for (auto pair : cu.get_initial_map_ref().left) {
      ids_updated &= grid.node_exists(Node(pair.second));
    }
    for (auto pair : cu.get_final_map_ref().left) {
      ids_updated &= grid.node_exists(Node(pair.second));
    }
    REQUIRE(ids_updated);
    Circuit res(cu.get_circ_ref());
    Vertex x = res.add_op<Qubit>(
        OpType::X, {Qubit(cu.get_final_map_ref().left.at(Qubit("qe")))});
    Vertex pred = res.get_predecessors(x).front();
    REQUIRE(res.get_OpType_from_Vertex(pred) == OpType::Z);
  }
}

SCENARIO("FlattenRegisters pass") {
  GIVEN("A simple circuit") {
    Circuit circ(3, 2);
    CompilationUnit cu(circ);
    REQUIRE_FALSE(FlattenRegisters()->apply(cu));
  }
  GIVEN("A non-simple circuit") {
    Circuit circ(2, 1);
    Qubit weird_qb("weird_q", 3);
    Qubit weird_qb2("weird_q", 5);
    Bit weird_cb("weird_c", 3, 1);
    circ.add_qubit(weird_qb);
    circ.add_qubit(weird_qb2);
    circ.add_bit(weird_cb);
    CompilationUnit cu(circ);
    REQUIRE(FlattenRegisters()->apply(cu));
    REQUIRE(cu.get_circ_ref().is_simple());
    unit_bimap_t map = cu.get_initial_map_ref();
    REQUIRE(map.left.at(weird_qb) == Qubit(2));
    REQUIRE(map.left.at(weird_qb2) == Qubit(3));
    REQUIRE(map.left.at(weird_cb) == Bit(1));
  }
}

SCENARIO("RemoveBarriers pass") {
  GIVEN("A circuit without a Barrier") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    CompilationUnit cu(circ);
    REQUIRE_FALSE(RemoveBarriers()->apply(cu));
    const Circuit& circ1 = cu.get_circ_ref();
    REQUIRE(circ1 == circ);
  }
  GIVEN("A circuit with a Barrier") {
    Circuit circ(3);
    add_1qb_gates(circ, OpType::H, {0, 1});
    circ.add_barrier({1, 2});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    CompilationUnit cu(circ);
    REQUIRE(RemoveBarriers()->apply(cu));
    const Circuit& circ1 = cu.get_circ_ref();
    REQUIRE(circ1.n_vertices() < circ.n_vertices());
  }
}

SCENARIO("gen_placement_pass test") {
  GIVEN("A simple circuit and device and base Placement.") {
    Circuit circ(4);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {2, 1}, {2, 3}});
    Architecture arc({{0, 1}, {1, 2}, {3, 2}});
    PlacementPtr plptr = std::make_shared<Placement>(arc);
    PassPtr pp_place = gen_placement_pass(plptr);
    CompilationUnit cu(circ);
    pp_place->apply(cu);
    Circuit res(cu.get_circ_ref());
    qubit_vector_t all_res_qbs = res.all_qubits();
    for (unsigned nn = 0; nn <= 3; ++nn) {
      REQUIRE(all_res_qbs[nn] == Qubit(Placement::unplaced_reg(), nn));
    }
  }
  GIVEN("A simple circuit and device and GraphPlacement.") {
    Circuit circ(4);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {2, 1}, {2, 3}});
    Architecture arc({{0, 1}, {1, 2}, {3, 2}});
    PlacementPtr plptr = std::make_shared<GraphPlacement>(arc);
    PassPtr pp_place = gen_placement_pass(plptr);
    CompilationUnit cu(circ);
    pp_place->apply(cu);
    Circuit res(cu.get_circ_ref());
    qubit_vector_t all_res_qbs = res.all_qubits();
    for (unsigned nn = 0; nn <= 3; ++nn) {
      REQUIRE(all_res_qbs[nn] == Node(nn));
    }
  }

  GIVEN("A large circuit and a large architecture.") {
    unsigned N = 150;
    Circuit circ(N);
    for (unsigned i = 0; i < N - 3; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, i + 1});
      circ.add_op<unsigned>(OpType::CX, {i, i + 2});
      circ.add_op<unsigned>(OpType::CX, {i, i + 3});
    }
    // Generate a line architecture
    std::vector<std::pair<unsigned, unsigned>> edges;
    for (unsigned i = 0; i < N - 1; i++) {
      edges.push_back({i, i + 1});
    }
    Architecture line_arc(edges);
    // Get a graph placement
    PassPtr graph_place =
        gen_placement_pass(std::make_shared<GraphPlacement>(line_arc));
    CompilationUnit graph_cu((Circuit(circ)));
    graph_place->apply(graph_cu);
    // Get a noise-aware placement
    PassPtr noise_place =
        gen_placement_pass(std::make_shared<NoiseAwarePlacement>(line_arc));
    CompilationUnit noise_cu((Circuit(circ)));
    noise_place->apply(noise_cu);
    // Get a line placement
    PassPtr line_place =
        gen_placement_pass(std::make_shared<LinePlacement>(line_arc));
    CompilationUnit line_cu((Circuit(circ)));
    line_place->apply(line_cu);
    // Get a fall back placement from a graph placement
    PlacementConfig config(5, line_arc.n_connections(), 10000, 10, 0);
    PassPtr graph_fall_back_place =
        gen_placement_pass(std::make_shared<GraphPlacement>(line_arc, config));
    CompilationUnit graph_fall_back_cu((Circuit(circ)));
    graph_fall_back_place->apply(graph_fall_back_cu);
    // Get a fall back placement from a noise-aware placement
    PassPtr noise_fall_back_place = gen_placement_pass(
        std::make_shared<NoiseAwarePlacement>(line_arc, config));
    CompilationUnit noise_fall_back_cu((Circuit(circ)));
    noise_fall_back_place->apply(noise_fall_back_cu);

    REQUIRE(graph_cu.get_final_map_ref() != line_cu.get_final_map_ref());
    REQUIRE(noise_cu.get_final_map_ref() != line_cu.get_final_map_ref());
    REQUIRE(
        graph_fall_back_cu.get_final_map_ref() == line_cu.get_final_map_ref());
    REQUIRE(
        noise_fall_back_cu.get_final_map_ref() == line_cu.get_final_map_ref());
  }
}

SCENARIO("gen_rename_qubits_pass test") {
  GIVEN("A circuit and a mapping") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    std::map<Qubit, Qubit> qm;
    Qubit newq0("newq0", 0);
    Qubit newq1("newq1", 1);
    Qubit newq2("newq2", 0);
    qm[Qubit(0)] = newq0;
    qm[Qubit(1)] = newq1;
    qm[Qubit(2)] = newq2;
    PassPtr pp = gen_rename_qubits_pass(qm);
    CompilationUnit cu(circ);
    REQUIRE(pp->apply(cu));
    const Circuit& newcirc = cu.get_circ_ref();
    Command cmd = newcirc.get_commands()[0];
    REQUIRE(cmd.get_args()[0] == newq0);
    REQUIRE(cmd.get_args()[1] == newq1);
  }
}

SCENARIO("PeepholeOptimise2Q and FullPeepholeOptimise") {
  GIVEN("A circuit with a Reset operation.") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::ZZMax, {1, 0});
    circ.add_op<unsigned>(OpType::CH, {0, 1});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::Reset, {1});
    CompilationUnit cu(circ);
    REQUIRE(PeepholeOptimise2Q()->apply(cu));
    Circuit circ1 = circ;
    CompilationUnit cu1(circ1);
    REQUIRE(FullPeepholeOptimise()->apply(cu1));
  }
  GIVEN("A symbolic circuit") {
    Sym a = SymEngine::symbol("alpha");
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::ZZMax, {1, 0});
    circ.add_op<unsigned>(OpType::CH, {0, 1});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::Ry, 2 * Expr(a), {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    CompilationUnit cu(circ);
    REQUIRE(PeepholeOptimise2Q()->apply(cu));
    Circuit circ1 = circ;
    CompilationUnit cu1(circ1);
    REQUIRE(FullPeepholeOptimise()->apply(cu1));
  }
  GIVEN("YYPhase") {
    // TKET-1302
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::YYPhase, 1.00378, {0, 1});
    circ.add_op<unsigned>(OpType::CV, {0, 1});
    circ.add_op<unsigned>(OpType::CSX, {1, 0});
    CompilationUnit cu(circ);
    REQUIRE(PeepholeOptimise2Q()->apply(cu));
    REQUIRE(test_unitary_comparison(circ, cu.get_circ_ref()));
    Circuit circ1 = circ;
    CompilationUnit cu1(circ1);
    REQUIRE(FullPeepholeOptimise()->apply(cu1));
    REQUIRE(test_unitary_comparison(circ, cu.get_circ_ref()));
  }
  GIVEN("An X + BRIDGE circuit") {
    // https://github.com/CQCL/tket/issues/9
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::BRIDGE, {0, 1, 2});
    CompilationUnit cu(circ);
    REQUIRE(FullPeepholeOptimise()->apply(cu));
    REQUIRE(test_unitary_comparison(circ, cu.get_circ_ref()));
  }
}

SCENARIO("rebase and decompose PhasePolyBox test") {
  GIVEN("rebase and decompose I") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::Y, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    CompilationUnit cu(circ);
    REQUIRE(ComposePhasePolyBoxes()->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("rebase and decompose II") {
    Circuit circ(2, 2);
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::Measure, {0, 1});
    circ.add_op<unsigned>(OpType::Z, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rx, 0.3, {1});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    CompilationUnit cu(circ);
    REQUIRE(ComposePhasePolyBoxes()->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(result.count_gates(OpType::CX) == 0);
    REQUIRE(result.count_gates(OpType::Rz) == 0);
    REQUIRE(result.count_gates(OpType::X) == 0);
    REQUIRE(result.count_gates(OpType::H) == 4);
    REQUIRE(result.count_gates(OpType::Measure) == 2);
    REQUIRE(result.count_gates(OpType::PhasePolyBox) == 4);
  }
  GIVEN("rebase and decompose III") {
    Circuit circ(8);
    circ.add_op<unsigned>(OpType::Rz, 0.5, {1});
    circ.add_op<unsigned>(OpType::Rx, 1.5, {2});
    circ.add_op<unsigned>(OpType::CX, {1, 7});
    circ.add_op<unsigned>(OpType::CX, {2, 4});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {1});
    circ.add_op<unsigned>(OpType::Rx, 0.5, {2});
    circ.add_op<unsigned>(OpType::CX, {1, 3});
    circ.add_op<unsigned>(OpType::X, {3});
    circ.add_op<unsigned>(OpType::X, {4});
    circ.add_op<unsigned>(OpType::CX, {5, 6});
    circ.add_op<unsigned>(OpType::CX, {6, 7});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {2});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::Rz, 1.5, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 3});
    circ.add_op<unsigned>(OpType::CX, {5, 0});
    circ.add_op<unsigned>(OpType::H, {1});
    add_2qb_gates(circ, OpType::CX, {{1, 4}, {2, 4}, {4, 7}, {3, 0}});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {0});
    circ.add_op<unsigned>(OpType::CX, {6, 3});
    circ.add_op<unsigned>(OpType::Rz, 1.5, {0});
    circ.add_op<unsigned>(OpType::CX, {4, 0});

    CompilationUnit cu(circ);
    REQUIRE(ComposePhasePolyBoxes()->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(circ, result));
  }
  GIVEN("NoWireSwapsPredicate for ComposePhasePolyBoxes") {
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{0, 3}, {1, 4}});
    circ.add_op<unsigned>(OpType::SWAP, {3, 4});
    circ.add_op<unsigned>(OpType::Z, {3});

    REQUIRE(NoWireSwapsPredicate().verify(circ));
    circ.replace_SWAPs();
    REQUIRE(!NoWireSwapsPredicate().verify(circ));

    CompilationUnit cu(circ);
    REQUIRE(ComposePhasePolyBoxes()->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(result, circ));
    REQUIRE(NoWireSwapsPredicate().verify(result));
  }
  GIVEN("NoWireSwapsPredicate for ComposePhasePolyBoxes II") {
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{0, 3}, {1, 4}});
    circ.add_op<unsigned>(OpType::SWAP, {3, 4});
    circ.add_op<unsigned>(OpType::Z, {3});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::Z, {1});
    circ.add_op<unsigned>(OpType::Z, {4});

    REQUIRE(NoWireSwapsPredicate().verify(circ));
    circ.replace_SWAPs();
    REQUIRE(!NoWireSwapsPredicate().verify(circ));

    CompilationUnit cu(circ);
    REQUIRE(ComposePhasePolyBoxes()->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(result, circ));
    REQUIRE(NoWireSwapsPredicate().verify(result));
  }
  GIVEN("NoWireSwapsPredicate for ComposePhasePolyBoxes III") {
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{0, 3}, {1, 4}});
    circ.add_op<unsigned>(OpType::SWAP, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Z, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Z, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::SWAP, {1, 2});
    circ.add_op<unsigned>(OpType::SWAP, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});

    REQUIRE(NoWireSwapsPredicate().verify(circ));
    circ.replace_SWAPs();
    REQUIRE(!NoWireSwapsPredicate().verify(circ));

    CompilationUnit cu(circ);
    REQUIRE(ComposePhasePolyBoxes()->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(result, circ));
    REQUIRE(NoWireSwapsPredicate().verify(result));
  }
  GIVEN("NoWireSwapsPredicate for ComposePhasePolyBoxes min size") {
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{0, 3}, {1, 4}});
    circ.add_op<unsigned>(OpType::SWAP, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Z, {3});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Z, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::SWAP, {1, 2});
    circ.add_op<unsigned>(OpType::SWAP, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});

    REQUIRE(NoWireSwapsPredicate().verify(circ));
    circ.replace_SWAPs();
    REQUIRE(!NoWireSwapsPredicate().verify(circ));

    CompilationUnit cu(circ);
    REQUIRE(ComposePhasePolyBoxes(5)->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(result.count_gates(OpType::H) == 3);
    REQUIRE(result.count_gates(OpType::CX) == 2);
    REQUIRE(result.count_gates(OpType::SWAP) == 0);
    REQUIRE(result.count_gates(OpType::Z) == 0);
    REQUIRE(result.count_gates(OpType::PhasePolyBox) == 2);

    REQUIRE(test_unitary_comparison(result, circ));
    REQUIRE(NoWireSwapsPredicate().verify(result));
  }
  GIVEN("NoWireSwapsPredicate for ComposePhasePolyBoxes min size ii") {
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{0, 3}, {1, 4}});
    circ.add_op<unsigned>(OpType::SWAP, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Z, {3});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Z, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::SWAP, {1, 2});
    circ.add_op<unsigned>(OpType::SWAP, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});

    REQUIRE(NoWireSwapsPredicate().verify(circ));
    circ.replace_SWAPs();
    REQUIRE(!NoWireSwapsPredicate().verify(circ));

    CompilationUnit cu(circ);
    REQUIRE(ComposePhasePolyBoxes(6)->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(result.count_gates(OpType::H) == 3);
    REQUIRE(result.count_gates(OpType::CX) == 7);
    REQUIRE(result.count_gates(OpType::SWAP) == 0);
    REQUIRE(result.count_gates(OpType::Z) == 0);
    REQUIRE(result.count_gates(OpType::PhasePolyBox) == 1);

    REQUIRE(test_unitary_comparison(result, circ));
    REQUIRE(NoWireSwapsPredicate().verify(result));
  }
  GIVEN("NoWireSwapsPredicate for aas I") {
    std::vector<Node> nodes = {Node(0), Node(1), Node(2), Node(3), Node(4)};
    Architecture architecture(
        {{nodes[0], nodes[1]},
         {nodes[1], nodes[2]},
         {nodes[2], nodes[3]},
         {nodes[3], nodes[4]}});

    Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{0, 3}, {1, 4}});
    circ.add_op<unsigned>(OpType::SWAP, {3, 4});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Z, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Z, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::SWAP, {1, 2});
    circ.add_op<unsigned>(OpType::SWAP, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});

    REQUIRE(NoWireSwapsPredicate().verify(circ));
    circ.replace_SWAPs();
    REQUIRE(!NoWireSwapsPredicate().verify(circ));

    CompilationUnit cu(circ);

    REQUIRE(gen_full_mapping_pass_phase_poly(
                architecture, 1, aas::CNotSynthType::Rec)
                ->apply(cu));
    Circuit result = cu.get_circ_ref();

    REQUIRE(test_unitary_comparison(result, circ));
    REQUIRE(NoWireSwapsPredicate().verify(result));
  }
}

SCENARIO("DecomposeArbitrarilyControlledGates test") {
  GIVEN("A CCX gate") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CCX, {2, 0, 1});
    CompilationUnit cu(circ);
    REQUIRE(DecomposeArbitrarilyControlledGates()->apply(cu));
  }
}

SCENARIO("Precomposed passes successfully compose") {
  GIVEN("gen_directed_cx_routing_pass") {
    RingArch arc(6);
    REQUIRE_NOTHROW(gen_directed_cx_routing_pass(
        arc, {std::make_shared<LexiLabellingMethod>(),
              std::make_shared<LexiRouteRoutingMethod>()}));
  }
}

SCENARIO("Test Pauli Graph Synthesis Pass") {
  PassPtr graph_synth = gen_synthesise_pauli_graph(
      Transforms::PauliSynthStrat::Sets, CXConfigType::Star);
  Circuit circ(3);
  PauliExpBox peb({Pauli::Z, Pauli::X, Pauli::Z}, 0.333);
  circ.add_box(peb, {0, 1, 2});
  PauliExpBox peb2({Pauli::Y, Pauli::X, Pauli::X}, 0.174);
  circ.add_box(peb2, {0, 1, 2});

  CompilationUnit cu(circ);
  graph_synth->apply(cu);

  REQUIRE(test_unitary_comparison(circ, cu.get_circ_ref()));
}

SCENARIO("Compose Pauli Graph synthesis Passes") {
  RingArch arc(10);
  PassPtr dir_pass = gen_directed_cx_routing_pass(
      arc, {std::make_shared<LexiLabellingMethod>(),
            std::make_shared<LexiRouteRoutingMethod>()});
  GIVEN("Special UCC Synthesis") {
    PassPtr spec_ucc = gen_special_UCC_synthesis();
    REQUIRE_NOTHROW(spec_ucc >> dir_pass);
  }
  GIVEN("Pauli Graph synthesis") {
    PassPtr graph_synth = gen_synthesise_pauli_graph(
        Transforms::PauliSynthStrat::Sets, CXConfigType::Star);
    REQUIRE_NOTHROW(graph_synth >> dir_pass);
  }
  GIVEN("Pairwise Synthesis") {
    PassPtr pairwise = gen_pairwise_pauli_gadgets(CXConfigType::Tree);
    REQUIRE_NOTHROW(pairwise >> dir_pass);
  }
}

SCENARIO("Commute measurements to the end of a circuit") {
  PassPtr delay_pass = DelayMeasures();
  PredicatePtr mid_meas_pred = std::make_shared<NoMidMeasurePredicate>();
  GIVEN("Measurements already at end") {
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::Z, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Rx, 0.3, {1});
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    c.add_op<unsigned>(OpType::Measure, {1, 1});
    CompilationUnit cu(c);
    REQUIRE_FALSE(delay_pass->apply(cu));
    REQUIRE(mid_meas_pred->verify(cu.get_circ_ref()));
  }
  GIVEN("Gates after measure") {
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::Measure, {0, 1});
    c.add_op<unsigned>(OpType::Z, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::Rx, 0.3, {1});
    c.add_op<unsigned>(OpType::SWAP, {0, 1});
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    CompilationUnit cu(c);
    REQUIRE(delay_pass->apply(cu));
    REQUIRE(mid_meas_pred->verify(cu.get_circ_ref()));
    Circuit expected(2, 2);
    expected.add_op<unsigned>(OpType::Z, {0});
    expected.add_op<unsigned>(OpType::CX, {0, 1});
    expected.add_op<unsigned>(OpType::Rx, 0.3, {1});
    expected.add_op<unsigned>(OpType::SWAP, {0, 1});
    expected.add_op<unsigned>(OpType::Measure, {0, 0});
    expected.add_op<unsigned>(OpType::Measure, {1, 1});
    REQUIRE(cu.get_circ_ref() == expected);
  }
  GIVEN("Measure blocked by quantum gate") {
    Circuit c(1, 1);
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    c.add_op<unsigned>(OpType::Rx, 0.3, {0});
    CompilationUnit cu(c);
    REQUIRE_THROWS_AS(delay_pass->apply(cu), CircuitInvalidity);
  }
  GIVEN("Measure blocked by classical operation") {
    Circuit c(2, 1);
    add_2qb_gates(c, OpType::Measure, {{0, 0}, {1, 0}});
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    c.add_op<unsigned>(OpType::Measure, {1, 0});
    CompilationUnit cu(c);
    REQUIRE_THROWS_AS(delay_pass->apply(cu), CircuitInvalidity);
  }
  GIVEN("Measure blocked by conditional operation") {
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::CZ, {0, 1});
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    c.add_conditional_gate<unsigned>(OpType::Z, {}, {1}, {0}, 1);
    CompilationUnit cu(c);
    REQUIRE_THROWS_AS(delay_pass->apply(cu), CircuitInvalidity);
  }
  GIVEN("Combined with routing") {
    Circuit test(3, 1);
    add_2qb_gates(test, OpType::CX, {{0, 1}, {1, 2}});
    add_1qb_gates(test, OpType::X, {0, 0});
    test.add_measure(1, 0);
    add_1qb_gates(test, OpType::X, {2, 2});
    test.add_op<unsigned>(OpType::CX, {0, 2});

    Architecture line({{0, 1}, {1, 2}, {2, 3}});
    PlacementPtr pp = std::make_shared<LinePlacement>(line);
    PassPtr route_pass = gen_full_mapping_pass(
        line, pp,
        {std::make_shared<LexiLabellingMethod>(),
         std::make_shared<LexiRouteRoutingMethod>()});
    CompilationUnit cu(test);
    route_pass->apply(cu);
    REQUIRE(delay_pass->apply(cu));
    Command final_command = cu.get_circ_ref().get_commands()[7];
    OpType type = final_command.get_op_ptr()->get_type();
    REQUIRE(type == OpType::Measure);
    REQUIRE(final_command.get_args().front() == Node(3));
  }
}

SCENARIO("RemoveRedundancies and phase") {
  GIVEN("A trivial circuit with nonzero phase") {
    // TKET-679
    Circuit c(1);
    c.add_op<unsigned>(OpType::TK1, {1., 0., 1.}, {0});
    CompilationUnit cu(c);
    REQUIRE(RemoveRedundancies()->apply(cu));
    const Circuit& c1 = cu.get_circ_ref();
    REQUIRE(c1.get_commands().size() == 0);
    REQUIRE(equiv_val(c1.get_phase(), 1.));
  }
  GIVEN("A circuit with a trivial TK2 gate and nonzero phase") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::TK2, {0., 2., 4.}, {0, 1});
    CompilationUnit cu(c);
    REQUIRE(RemoveRedundancies()->apply(cu));
    const Circuit& c1 = cu.get_circ_ref();
    REQUIRE(c1.get_commands().size() == 0);
    REQUIRE(equiv_val(c1.get_phase(), 1.));
  }
  GIVEN("A circuit with a trivial TK2 gate and nonzero phase") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::TK2, {0., 2., 4.}, {0, 1});
    CompilationUnit cu(c);
    REQUIRE(RemoveRedundancies()->apply(cu));
    const Circuit& c1 = cu.get_circ_ref();
    REQUIRE(c1.get_commands().size() == 0);
    REQUIRE(equiv_val(c1.get_phase(), 1.));
  }
}

// Check whether a circuit maps all basis states to basis states.
// All compiler passes should preserve this property.
static bool is_classical_map(const Circuit& c) {
  const Eigen::MatrixXcd u = tket_sim::get_unitary(c);
  return std::all_of(
      u.data(), u.data() + u.size(), [](const std::complex<double>& x) {
        double r = std::abs(x);
        return r < ERR_EPS || r > 1 - ERR_EPS;
      });
}

SCENARIO("CX mapping pass") {
  // TKET-1045
  GIVEN("A device with a linear architecture") {
    Architecture line({{0, 1}, {1, 2}, {2, 3}, {3, 4}});

    // Noise-aware placement and rebase
    PlacementPtr placer = std::make_shared<NoiseAwarePlacement>(line);
    Circuit cx(2);
    cx.add_op<unsigned>(OpType::CX, {0, 1});
    OpTypeSet gateset = all_single_qubit_types();
    gateset.insert(OpType::CX);
    PassPtr rebase = gen_rebase_pass(gateset, cx, CircPool::tk1_to_tk1);

    // Circuit mapping basis states to basis states
    Circuit c(3);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::CCX, {2, 1, 0});
    c.add_op<unsigned>(OpType::CY, {1, 0});
    c.add_op<unsigned>(OpType::CY, {2, 1});
    REQUIRE(is_classical_map(c));

    // Rebase
    CompilationUnit cu_rebase(c);
    REQUIRE(rebase->apply(cu_rebase));
    const Circuit& c_rebased = cu_rebase.get_circ_ref();
    REQUIRE(is_classical_map(c_rebased));

    // Place
    CompilationUnit cu_place(c_rebased);
    gen_placement_pass(placer)->apply(cu_place);
    const Circuit& c_placed = cu_place.get_circ_ref();
    REQUIRE(is_classical_map(c_placed));

    // Route
    LexiRouteRoutingMethod lrrm(50);
    RoutingMethodPtr rmw = std::make_shared<LexiRouteRoutingMethod>(lrrm);
    CompilationUnit cu_route(c_placed);
    gen_routing_pass(
        line, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()})
        ->apply(cu_route);
    const Circuit& c_routed = cu_route.get_circ_ref();

    // Rebase again
    CompilationUnit cu(c_routed);
    rebase->apply(cu);
    const Circuit& c1 = cu.get_circ_ref();
    c1.assert_valid();
    REQUIRE(is_classical_map(c1));
  }
}

SCENARIO("ThreeQubitSquah") {
  GIVEN("A 3-qubit circuit that can be squashed") {
    Circuit c(3);
    for (unsigned i = 0; i < 21; i++) {
      c.add_op<unsigned>(OpType::H, {i % 3});
      c.add_op<unsigned>(OpType::CX, {i % 3, (i + 1) % 3});
      c.add_op<unsigned>(OpType::Rz, 0.25, {(i + 1) % 3});
    }
    CompilationUnit cu(c);
    REQUIRE(ThreeQubitSquash()->apply(cu));
    const Circuit& c1 = cu.get_circ_ref();
    REQUIRE(c1.count_gates(OpType::CX) <= 19);
    REQUIRE(test_statevector_comparison(c, c1));
  }
  GIVEN("A 3-qubit circuit that is non-trivially the identity") {
    Circuit c(3);
    c.add_op<unsigned>(OpType::U3, {1.5, 0, 1.5}, {0});
    c.add_op<unsigned>(OpType::U3, {0.5, 0.75, 1.25}, {1});
    c.add_op<unsigned>(OpType::U3, {0.5, 0, 1}, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::U3, {0.25, 0.25, 1.75}, {1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::U3, {0.5, 0, 0.25}, {1});
    c.add_op<unsigned>(OpType::U3, {3.5, 1.75, 0}, {2});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::U3, {0.5, 1.75, 0}, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::U1, 0.25, {0});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::U1, 0.25, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::U3, {1.5, 1.5, 1.75}, {0});
    c.add_op<unsigned>(OpType::U3, {0.5, 0.75, 1.25}, {1});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::U1, 0.5, {0});
    c.add_op<unsigned>(OpType::U3, {0.5, 0, 0.5}, {2});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::U3, {0.25, 0.25, 1.75}, {1});
    c.add_op<unsigned>(OpType::CX, {1, 2});
    c.add_op<unsigned>(OpType::U3, {0.5, 0, 0.25}, {1});
    c.add_op<unsigned>(OpType::U3, {3.5, 0.25, 0}, {2});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::U3, {0.5, 1.75, 0}, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::U1, 0.25, {0});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_op<unsigned>(OpType::U1, 0.25, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::U1, 0.25, {0});
    c.add_op<unsigned>(OpType::CX, {2, 0});
    c.add_phase(0.25);

    CompilationUnit cu(c);
    REQUIRE(ThreeQubitSquash()->apply(cu));
    const Circuit& c1 = cu.get_circ_ref();
    REQUIRE(c1.get_commands().empty());
  }
}

}  // namespace test_CompilerPass
}  // namespace tket
