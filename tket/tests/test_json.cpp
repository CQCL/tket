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

#include <boost/range/join.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>

#include "Architecture/Architecture.hpp"
#include "Circuit/CircPool.hpp"
#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Circuit/Command.hpp"
#include "CircuitsForTesting.hpp"
#include "Converters/PhasePoly.hpp"
#include "Gate/SymTable.hpp"
#include "Mapping/LexiLabelling.hpp"
#include "Mapping/LexiRoute.hpp"
#include "Mapping/RoutingMethod.hpp"
#include "MeasurementSetup/MeasurementSetup.hpp"
#include "OpType/OpType.hpp"
#include "Ops/OpPtr.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Predicates/PassLibrary.hpp"
#include "Transformations/OptimisationPass.hpp"
#include "Transformations/PauliOptimisation.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/Json.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_json {

template <class T>
bool serialize_deserialize(const T& obj) {
  nlohmann::json j = obj;
  auto new_obj = j.get<T>();
  return obj == new_obj;
}

template <class T>
void check_cases(const std::vector<T>& cases) {
  for (const auto& test : cases) {
    CHECK(serialize_deserialize(test));
  }
}

bool check_circuit(const Circuit& c) {
  nlohmann::json j = c;
  auto new_c = j.get<Circuit>();
  return c.circuit_equality(new_c);
}

SCENARIO("Test Op serialization") {
  GIVEN("OpType") {
    const OpTypeSet metaops = {
        OpType::Input, OpType::Output, OpType::ClInput, OpType::ClOutput,
        OpType::Barrier};
    const OpTypeSet boxes = {OpType::CircBox,      OpType::Unitary1qBox,
                             OpType::Unitary2qBox, OpType::Unitary3qBox,
                             OpType::ExpBox,       OpType::PauliExpBox,
                             OpType::CustomGate,   OpType::CliffBox,
                             OpType::PhasePolyBox, OpType::QControlBox};

    std::set<std::string> type_names;
    for (auto type :
         boost::join(all_gate_types(), boost::join(metaops, boxes))) {
      bool success_insert =
          type_names.insert(optypeinfo().at(type).name).second;
      // check all optype names are unique
      CHECK(success_insert);
      CHECK(serialize_deserialize(type));
    }

    nlohmann::json false_str = "NOTANOPTYPE";
    nlohmann::json correct_str = "Z";
    CHECK(correct_str.get<OpType>() == OpType::Z);
    REQUIRE_THROWS_AS(false_str.get<OpType>(), JsonError);
  }

  GIVEN("Expressions") {
    std::vector<Expr> e_tests = {
        Expr(0.3), Expr("a"), Expr((2 * 3. / 4 - 1)),
        Expr(-0.3 + (3.4 * Expr(SymEngine::sin(Expr("d") - 2.3))))};
    check_cases(e_tests);
  }
}

SCENARIO("Test UnitID serialization") {
  GIVEN("A list of Qubit instances") {
    std::vector<Qubit> test_q = {
        Qubit("test", 1), Qubit(4), Node(3), Qubit("a", {1, 2, 3, 4}),
        Qubit("sdaf", 1, 2)};
    check_cases(test_q);

    std::vector<Bit> test_b = {
        Bit("test", 1), Bit(4), Bit("a", {1, 2, 3, 4}), Bit("sdaf", 1, 2)};
    check_cases(test_b);
  }
}

SCENARIO("Test Command serialization") {
  GIVEN("A test circuit") {
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::Rz, 0.2, {0});
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Measure, {0, 1});
    const qubit_vector_t& q = c.all_qubits();
    const Qubit a = Qubit("a", 1, 2);
    c.add_qubit(a);
    c.add_op<UnitID>(OpType::CnRy, 0.1, {q[0], a, q[1]});
    c.add_barrier({q[0], a});

    check_cases(c.get_commands());
  }
}

SCENARIO("Test Circuit serialization") {
  GIVEN("A simple test circuit") {
    Circuit c(2, 2, "test_circ_1");
    // add standard ops
    c.add_op<unsigned>(OpType::Rz, 0.2, {0});
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Measure, {0, 1});
    // add custom qubit
    const qubit_vector_t& q = c.all_qubits();
    const Qubit a = Qubit("a", 1, 2);
    c.add_qubit(a);
    // variable arity gate
    c.add_op<UnitID>(OpType::CnRy, 0.1, {q[0], a, q[1]});
    // barrier metaop
    c.add_barrier({q[0], a});
    // phase
    c.add_phase(0.3);
    REQUIRE(check_circuit(c));
  }

  GIVEN("An implicit permutation") {
    Circuit circ(3);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 0}, {1, 2}, {2, 1}});

    Transforms::clifford_simp().apply(circ);

    REQUIRE(check_circuit(circ));
  }

  GIVEN("Conditional") {
    Circuit c(2, 3);
    c.add_conditional_gate<unsigned>(OpType::Ry, {-0.75}, {0}, {0, 1}, 1);
    c.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0, 1}, 1);
    c.add_conditional_gate<unsigned>(OpType::Measure, {}, {0, 2}, {0, 1}, 1);

    nlohmann::json j_box = c;
    const Circuit new_c = j_box.get<Circuit>();
    c.circuit_equality(new_c);
    const auto& cond =
        static_cast<const Conditional&>(*new_c.get_commands()[1].get_op_ptr());

    REQUIRE(*cond.get_op() == *get_op_ptr(OpType::CX));
    REQUIRE(cond.get_width() == 2);
    REQUIRE(cond.get_value() == 1);
  }

  GIVEN("A circbox") {
    Circuit c(3, 2, "circbox_base");
    c.add_op<unsigned>(OpType::Rz, 0.2, {0});

    Circuit temp_circ(2, "circbox");
    temp_circ.add_op<unsigned>(OpType::Ry, 0.75, {0});
    temp_circ.add_op<unsigned>(OpType::CX, {0, 1});
    CircBox temp_box(temp_circ);
    c.add_box(temp_box, {0, 1});

    nlohmann::json j_cbox = c;
    const Circuit new_c = j_cbox.get<Circuit>();

    const Command cbox_com = new_c.get_commands()[1];
    const CircBox& c_b = static_cast<const CircBox&>(*cbox_com.get_op_ptr());
    REQUIRE(temp_box == c_b);
    const Circuit& new_temp = *c_b.to_circuit();
    CHECK(new_temp.get_name() == temp_circ.get_name());
    REQUIRE(new_temp == temp_circ);
  }

  GIVEN("Unitary Boxes") {
    Circuit c(3, 2, "unitarybox");
    c.add_op<unsigned>(OpType::Rz, 0.2, {0});

    Circuit setup(1);
    setup.add_op<unsigned>(OpType::TK1, {0.2374, 1.0353, 0.5372}, {0});
    Eigen::Matrix2cd m = get_matrix_from_circ(setup);
    Unitary1qBox mbox(m);
    c.add_box(mbox, {1});

    Eigen::Matrix4cd m2;
    m2 << 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0;
    Unitary2qBox mbox2(m2);
    c.add_box(mbox2, {0, 2});

    Matrix8cd U = Matrix8cd::Zero();
    U(0, 3) = 1;
    U(1, 1) = 1;
    U(2, 7) = 1;
    U(3, 5) = 1;
    U(4, 0) = 1;
    U(5, 4) = 1;
    U(6, 2) = 1;
    U(7, 6) = 1;
    Unitary3qBox mbox3(U);

    c.add_box(mbox3, {0, 1, 2});

    Eigen::Matrix4cd A;
    A << 0., 1., 2., 3., 1., 2., 3. * i_, 4., 2., -3. * i_, 3, 2. - 3. * i_, 3.,
        4., 2. + 3. * i_, 5.;
    ExpBox ebox(A, -0.5);
    c.add_box(ebox, {1, 2});

    nlohmann::json j_mbox = c;
    const Circuit new_c = j_mbox.get<Circuit>();

    const std::vector<Command> coms = new_c.get_commands();
    const auto& m_b = static_cast<const Unitary1qBox&>(*coms[1].get_op_ptr());
    REQUIRE(matrices_are_equal(mbox.get_matrix(), m_b.get_matrix()));
    REQUIRE(mbox == m_b);
    const auto& m2_b = static_cast<const Unitary2qBox&>(*coms[2].get_op_ptr());
    REQUIRE(matrices_are_equal(mbox2.get_matrix(), m2_b.get_matrix()));
    REQUIRE(mbox2 == m2_b);

    const auto& m3_b = static_cast<const Unitary3qBox&>(*coms[3].get_op_ptr());
    REQUIRE(matrices_are_equal(mbox3.get_matrix(), m3_b.get_matrix()));
    REQUIRE(mbox3 == m3_b);

    const auto& exp_b = static_cast<const ExpBox&>(*coms[4].get_op_ptr());

    const auto ebox_m_p = ebox.get_matrix_and_phase();
    const auto exp_b_m_p = exp_b.get_matrix_and_phase();
    REQUIRE(matrices_are_equal(ebox_m_p.first, exp_b_m_p.first));
    REQUIRE(ebox_m_p.second == exp_b_m_p.second);

    REQUIRE(ebox == exp_b);
  }
  GIVEN("Pauli ExpBoxes") {
    Circuit c(4, 2, "paulibox");
    PauliExpBox pbox({Pauli::X, Pauli::Y, Pauli::I, Pauli::Z}, -0.72521);
    c.add_box(pbox, {0, 1, 2, 3});
    nlohmann::json j_pbox = c;
    const Circuit new_c = j_pbox.get<Circuit>();

    const auto& p_b =
        static_cast<const PauliExpBox&>(*new_c.get_commands()[0].get_op_ptr());

    REQUIRE(p_b.get_paulis() == pbox.get_paulis());
    REQUIRE(p_b.get_phase() == pbox.get_phase());
    REQUIRE(p_b == pbox);
  }

  GIVEN("CustomGate") {
    Circuit setup(2);
    Sym a = SymTable::fresh_symbol("a");
    Sym c = SymTable::fresh_symbol("c");
    Expr b(SymTable::fresh_symbol("b"));
    setup.add_op<unsigned>(OpType::Rx, {c}, {0});
    setup.add_op<unsigned>(OpType::CX, {0, 1});
    setup.add_op<unsigned>(OpType::Ry, {a}, {0});
    composite_def_ptr_t def = CompositeGateDef::define_gate("g", setup, {a});
    CustomGate g0(def, {0.2374});
    CustomGate g1(def, {b});

    Circuit circ(3);
    circ.add_box(g0, {0, 1});
    circ.add_box(g1, {1, 2});

    nlohmann::json j_pbox = circ;
    const Circuit new_c = j_pbox.get<Circuit>();

    const std::vector<Command> coms = new_c.get_commands();

    const auto& g_0_new = static_cast<const CustomGate&>(*coms[0].get_op_ptr());

    REQUIRE(g0.get_params() == g_0_new.get_params());
    REQUIRE(*g0.get_gate() == *g_0_new.get_gate());
    REQUIRE(g0 == g_0_new);
    const auto& g_1_new = static_cast<const CustomGate&>(*coms[1].get_op_ptr());

    REQUIRE(g1.get_params() == g_1_new.get_params());
    REQUIRE(*g1.get_gate() == *g_1_new.get_gate());
    REQUIRE(g1 == g_1_new);
  }

  GIVEN("QControlBox") {
    Op_ptr op = get_op_ptr(OpType::Sycamore);
    QControlBox qcbox(op, 2);
    Circuit c(4);
    c.add_box(qcbox, {0, 1, 2, 3});

    nlohmann::json j_box = c;
    const Circuit new_c = j_box.get<Circuit>();

    const auto& qc_b =
        static_cast<const QControlBox&>(*new_c.get_commands()[0].get_op_ptr());

    REQUIRE(qc_b == qcbox);
    REQUIRE(qc_b.get_n_controls() == qcbox.get_n_controls());
    REQUIRE(*qc_b.get_op() == *qc_b.get_op());
  }

  GIVEN("PhasePolyBox") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    PhasePolyBox ppbox(circ);
    Circuit c(3);
    c.add_box(ppbox, {1, 2});
    nlohmann::json j_box = c;

    const Circuit new_c = j_box.get<Circuit>();

    const auto& pp_b =
        static_cast<const PhasePolyBox&>(*new_c.get_commands()[0].get_op_ptr());

    // can use box equality check here as in this case all members are checked
    REQUIRE(pp_b == ppbox);
  }

  GIVEN("Circuits with named operations") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.125, {1}, "foo");
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Circuit circ1(2);
    circ1.add_op<unsigned>(OpType::CX, {0, 1}, "bar");
    circ1.add_op<unsigned>(OpType::Rz, 0.125, {1});
    circ1.add_op<unsigned>(OpType::CX, {0, 1}, "bar");
    REQUIRE(check_circuit(circ));
    REQUIRE(check_circuit(circ1));
    REQUIRE(!(circ == circ1));
  }
}

SCENARIO("Test config serializations") {
  GIVEN("PlacementConfig") {
    PlacementConfig orig(5, 20, 100000, 10, 1);
    nlohmann::json j_config = orig;
    PlacementConfig loaded = j_config.get<PlacementConfig>();
    REQUIRE(orig == loaded);
    nlohmann::json j_loaded = loaded;
    REQUIRE(j_config == j_loaded);
  }
}

SCENARIO("Test device serializations") {
  GIVEN("Architecture") {
    Architecture arc({{0, 1}, {1, 2}});
    nlohmann::json j_arc = arc;
    Architecture loaded_arc = j_arc.get<Architecture>();
    CHECK(arc == loaded_arc);
    nlohmann::json j_loaded_arc = loaded_arc;
    CHECK(j_arc == j_loaded_arc);
    Architecture ring = RingArch(6);
    node_vector_t nodes = ring.get_all_nodes_vec();
    ring.add_connection(nodes.at(0), nodes.at(3), 20);
    nlohmann::json j_ring = ring;
    Architecture loaded_ring = j_ring.get<Architecture>();
    CHECK(ring == loaded_ring);
    nlohmann::json j_loaded_ring = loaded_ring;
    CHECK(j_ring == j_loaded_ring);
  }
  GIVEN("FullyConnected") {
    FullyConnected full(4);
    nlohmann::json j_full = full;
    FullyConnected loaded_full = j_full.get<FullyConnected>();
    CHECK(full == loaded_full);
    nlohmann::json j_loaded_full = loaded_full;
    CHECK(j_full == j_loaded_full);
  }
  GIVEN("DeviceCharacterisation") {
    Architecture ring = RingArch(3);
    node_vector_t nodes = ring.get_all_nodes_vec();
    op_errors_t node_err0{{{OpType::X, 0.3}, {OpType::Y, 0.4}}};
    op_errors_t node_err1{{{OpType::X, 0.2}, {OpType::Y, 0.5}}};
    op_node_errors_t ne{
        {nodes.at(0), node_err0},
        {nodes.at(1), node_err1},
        {nodes.at(2), node_err1}};
    op_errors_t link_err0{{{OpType::CX, 0.1}}};
    op_errors_t link_err1{{{OpType::CX, 0.1}, {OpType::CZ, 0.2}}};
    op_link_errors_t le{
        {{nodes.at(0), nodes.at(1)}, link_err0},
        {{nodes.at(1), nodes.at(2)}, link_err1},
        {{nodes.at(0), nodes.at(2)}, link_err0}};
    avg_readout_errors_t roe{
        {nodes.at(0), 0.02}, {nodes.at(1), 0.01}, {nodes.at(2), 0.98}};
    DeviceCharacterisation op_dc{ne, le, roe};
    nlohmann::json j_op_dc = op_dc;
    DeviceCharacterisation loaded_op_dc = j_op_dc.get<DeviceCharacterisation>();
    CHECK(op_dc == loaded_op_dc);
    nlohmann::json j_loaded_op_dc = loaded_op_dc;
    CHECK(j_op_dc == j_loaded_op_dc);
    avg_node_errors_t avg_ne{
        {nodes.at(0), 0.}, {nodes.at(1), 0.1}, {nodes.at(2), 0.2}};
    avg_link_errors_t avg_le{
        {{nodes.at(0), nodes.at(1)}, 0.},
        {{nodes.at(1), nodes.at(2)}, 0.1},
        {{nodes.at(1), nodes.at(2)}, 0.9}};
    DeviceCharacterisation avg_dc{avg_ne, avg_le, roe};
    nlohmann::json j_avg_dc = avg_dc;
    DeviceCharacterisation loaded_avg_dc =
        j_avg_dc.get<DeviceCharacterisation>();
    CHECK(avg_dc == loaded_avg_dc);
    nlohmann::json j_loaded_avg_dc = loaded_avg_dc;
    CHECK(j_avg_dc == j_loaded_avg_dc);
  }
}

SCENARIO("Test RoutingMethod serializations") {
  RoutingMethod rm;
  nlohmann::json rm_j = rm;
  RoutingMethod loaded_rm_j = rm_j.get<RoutingMethod>();

  Circuit c(2, 2);
  c.add_op<unsigned>(OpType::CX, {0, 1});

  MappingFrontier mf(c);
  MappingFrontier_ptr mf_sp = std::make_shared<MappingFrontier>(mf);
  CHECK(!loaded_rm_j.routing_method(mf_sp, std::make_shared<SquareGrid>(2, 2))
             .first);

  std::vector<RoutingMethodPtr> rmp = {
      std::make_shared<RoutingMethod>(rm),
      std::make_shared<LexiLabellingMethod>(),
      std::make_shared<LexiRouteRoutingMethod>(5)};

  nlohmann::json rmp_j = rmp;
  std::vector<RoutingMethodPtr> loaded_rmp_j =
      rmp_j.get<std::vector<RoutingMethodPtr>>();
  CHECK(!loaded_rmp_j[0]
             ->routing_method(mf_sp, std::make_shared<SquareGrid>(2, 2))
             .first);
  CHECK(loaded_rmp_j[1]
            ->routing_method(mf_sp, std::make_shared<SquareGrid>(2, 2))
            .first);
}

SCENARIO("Test predicate serializations") {
#define BASICPREDJSONTEST(classname)                             \
  GIVEN(#classname) {                                            \
    PredicatePtr pp = std::make_shared<classname>();             \
    nlohmann::json j_pp = pp;                                    \
    PredicatePtr loaded_pp = j_pp.get<PredicatePtr>();           \
    REQUIRE_NOTHROW(dynamic_cast<const classname&>(*loaded_pp)); \
    nlohmann::json j_loaded_pp = loaded_pp;                      \
    REQUIRE(j_pp == j_loaded_pp);                                \
  }
  BASICPREDJSONTEST(NoClassicalControlPredicate)
  BASICPREDJSONTEST(NoFastFeedforwardPredicate)
  BASICPREDJSONTEST(NoClassicalBitsPredicate)
  BASICPREDJSONTEST(NoWireSwapsPredicate)
  BASICPREDJSONTEST(MaxTwoQubitGatesPredicate)
  BASICPREDJSONTEST(CliffordCircuitPredicate)
  BASICPREDJSONTEST(DefaultRegisterPredicate)
  BASICPREDJSONTEST(NoBarriersPredicate)
  BASICPREDJSONTEST(NoMidMeasurePredicate)
  BASICPREDJSONTEST(NoSymbolsPredicate)
  BASICPREDJSONTEST(GlobalPhasedXPredicate)
#undef BASICPREDJSONTEST
  GIVEN("GateSetPredicate") {
    OpTypeSet ops = {OpType::X, OpType::V, OpType::Rz, OpType::ZZMax};
    PredicatePtr gs = std::make_shared<GateSetPredicate>(ops);
    nlohmann::json j_gs = gs;
    PredicatePtr loaded_gs = j_gs.get<PredicatePtr>();
    REQUIRE(
        dynamic_cast<const GateSetPredicate&>(*loaded_gs).get_allowed_types() ==
        ops);
    // Don't check the json equality here since ordering of elements in an
    // OpTypeSet (std::unordered_set<OpType>) is not guaranteed
  }
  GIVEN("PlacementPredicate") {
    node_set_t nodes = {Node(0), Node(14), Node(16)};
    PredicatePtr pl = std::make_shared<PlacementPredicate>(nodes);
    nlohmann::json j_pl = pl;
    PredicatePtr loaded_pl = j_pl.get<PredicatePtr>();
    REQUIRE(
        dynamic_cast<const PlacementPredicate&>(*loaded_pl).get_nodes() ==
        nodes);
    nlohmann::json j_loaded_pl = loaded_pl;
    REQUIRE(j_pl == j_loaded_pl);
  }
  GIVEN("ConnectivityPredicate") {
    Architecture ring = RingArch(3);
    PredicatePtr conn = std::make_shared<ConnectivityPredicate>(ring);
    nlohmann::json j_conn = conn;
    PredicatePtr loaded_conn = j_conn.get<PredicatePtr>();
    REQUIRE(
        dynamic_cast<const ConnectivityPredicate&>(*loaded_conn).get_arch() ==
        ring);
    nlohmann::json j_loaded_conn = loaded_conn;
    REQUIRE(j_conn == j_loaded_conn);
  }
  GIVEN("DirectednessPredicate") {
    Architecture ring = RingArch(3);
    PredicatePtr conn = std::make_shared<DirectednessPredicate>(ring);
    nlohmann::json j_conn = conn;
    PredicatePtr loaded_conn = j_conn.get<PredicatePtr>();
    REQUIRE(
        dynamic_cast<const DirectednessPredicate&>(*loaded_conn).get_arch() ==
        ring);
    nlohmann::json j_loaded_conn = loaded_conn;
    REQUIRE(j_conn == j_loaded_conn);
  }
  GIVEN("MaxNQubitsPredicate") {
    PredicatePtr max = std::make_shared<MaxNQubitsPredicate>(12);
    nlohmann::json j_max = max;
    PredicatePtr loaded_max = j_max.get<PredicatePtr>();
    REQUIRE(
        dynamic_cast<const MaxNQubitsPredicate&>(*loaded_max).get_n_qubits() ==
        12);
    nlohmann::json j_loaded_max = loaded_max;
    REQUIRE(j_max == j_loaded_max);
  }
  GIVEN("UserDefinedPredicate") {
    std::function<bool(const Circuit&)> func = [](const Circuit& c) {
      return false;
    };
    PredicatePtr custom = std::make_shared<UserDefinedPredicate>(func);
    nlohmann::json j_custom = custom;
    REQUIRE_THROWS_AS(j_custom.get<PredicatePtr>(), NotImplemented);
  }
}

SCENARIO("Test compiler pass serializations") {
  Architecture arc = SquareGrid(2, 4, 2);
  RoutingMethodPtr rmp = std::make_shared<LexiRouteRoutingMethod>(80);
  std::vector<RoutingMethodPtr> rcon = {rmp};
  PlacementConfig plcon(5, 20, 100000, 10, 1000);
  PlacementPtr place = std::make_shared<GraphPlacement>(arc, plcon);
  std::map<Qubit, Qubit> qmap = {{Qubit(0), Node(1)}, {Qubit(3), Node(2)}};
  PlacementPtr na_place = std::make_shared<NoiseAwarePlacement>(arc, plcon);
#define COMPPASSJSONTEST(passname, pass)               \
  GIVEN(#passname) {                                   \
    Circuit circ = CircuitsForTesting::get().uccsd;    \
    CompilationUnit cu{circ};                          \
    CompilationUnit copy = cu;                         \
    PassPtr pp = pass;                                 \
    nlohmann::json j_pp = pp;                          \
    PassPtr loaded = j_pp.get<PassPtr>();              \
    pp->apply(cu);                                     \
    loaded->apply(copy);                               \
    REQUIRE(cu.get_circ_ref() == copy.get_circ_ref()); \
    nlohmann::json j_loaded = loaded;                  \
    REQUIRE(j_pp == j_loaded);                         \
  }
  COMPPASSJSONTEST(CommuteThroughMultis, CommuteThroughMultis())
  COMPPASSJSONTEST(
      DecomposeArbitrarilyControlledGates,
      DecomposeArbitrarilyControlledGates())
  COMPPASSJSONTEST(DecomposeBoxes, DecomposeBoxes())
  COMPPASSJSONTEST(DecomposeMultiQubitsCX, DecomposeMultiQubitsCX())
  COMPPASSJSONTEST(DecomposeSingleQubitsTK1, DecomposeSingleQubitsTK1())
  COMPPASSJSONTEST(PeepholeOptimise2Q, PeepholeOptimise2Q())
  COMPPASSJSONTEST(FullPeepholeOptimise, FullPeepholeOptimise())
  COMPPASSJSONTEST(RebaseTket, RebaseTket())
  COMPPASSJSONTEST(RebaseUFR, RebaseUFR())
  COMPPASSJSONTEST(RemoveRedundancies, RemoveRedundancies())
  COMPPASSJSONTEST(SynthesiseHQS, SynthesiseHQS())
  COMPPASSJSONTEST(SynthesiseTK, SynthesiseTK())
  COMPPASSJSONTEST(SynthesiseTket, SynthesiseTket())
  COMPPASSJSONTEST(SynthesiseOQC, SynthesiseOQC())
  COMPPASSJSONTEST(SynthesiseUMD, SynthesiseUMD())
  COMPPASSJSONTEST(SquashTK1, SquashTK1())
  COMPPASSJSONTEST(FlattenRegisters, FlattenRegisters())
  COMPPASSJSONTEST(DelayMeasures, DelayMeasures())
  COMPPASSJSONTEST(RemoveDiscarded, RemoveDiscarded())
  COMPPASSJSONTEST(SimplifyMeasured, SimplifyMeasured())
  COMPPASSJSONTEST(RemoveBarriers, RemoveBarriers())
  COMPPASSJSONTEST(ComposePhasePolyBoxes, ComposePhasePolyBoxes())
  COMPPASSJSONTEST(DecomposeBridges, DecomposeBridges())
  COMPPASSJSONTEST(KAKDecomposition, KAKDecomposition(0.98))
  COMPPASSJSONTEST(ThreeQubitSquash, ThreeQubitSquash(false))
  COMPPASSJSONTEST(
      EulerAngleReduction, gen_euler_pass(OpType::Rx, OpType::Ry, false))
  COMPPASSJSONTEST(RenameQubitsPass, gen_rename_qubits_pass(qmap))
  COMPPASSJSONTEST(CliffordSimp, gen_clifford_simp_pass(true))
  COMPPASSJSONTEST(
      DecomposeSwapsToCXs, gen_decompose_routing_gates_to_cxs_pass(arc, false))
  COMPPASSJSONTEST(
      DecomposeSwapsToCircuit,
      gen_user_defined_swap_decomp_pass(CircPool::SWAP_using_CX_1()))
  COMPPASSJSONTEST(
      OptimisePhaseGadgets, gen_optimise_phase_gadgets(CXConfigType::Star))
  COMPPASSJSONTEST(
      OptimisePairwiseGadgets, gen_pairwise_pauli_gadgets(CXConfigType::Tree))
  COMPPASSJSONTEST(
      PauliSimp, gen_synthesise_pauli_graph(
                     Transforms::PauliSynthStrat::Sets, CXConfigType::Tree))
  COMPPASSJSONTEST(
      GuidedPauliSimp,
      gen_special_UCC_synthesis(
          Transforms::PauliSynthStrat::Pairwise, CXConfigType::Snake))
  COMPPASSJSONTEST(
      SimplifyInitial,
      gen_simplify_initial(
          Transforms::AllowClassical::No, Transforms::CreateAllQubits::Yes,
          std::make_shared<Circuit>(CircPool::X())))
  COMPPASSJSONTEST(PlacementPass, gen_placement_pass(place))
  // TKET-1419
  COMPPASSJSONTEST(NoiseAwarePlacement, gen_placement_pass(na_place))
  COMPPASSJSONTEST(NaivePlacementPass, gen_naive_placement_pass(arc))
#undef COMPPASSJSONTEST
  GIVEN("RoutingPass") {
    // Can only be applied to placed circuits
    Circuit circ = CircuitsForTesting::get().uccsd;
    CompilationUnit cu{circ};
    PassPtr placement = gen_placement_pass(place);
    placement->apply(cu);
    CompilationUnit copy = cu;
    PassPtr pp = gen_routing_pass(arc, rcon);
    nlohmann::json j_pp = pp;
    PassPtr loaded = j_pp.get<PassPtr>();
    pp->apply(cu);
    loaded->apply(copy);
    REQUIRE(cu.get_circ_ref() == copy.get_circ_ref());
    nlohmann::json j_loaded = loaded;
    REQUIRE(j_pp == j_loaded);
  }
  GIVEN("Routing with multiple routing methods") {
    RoutingMethodPtr mrmp =
        std::make_shared<MultiGateReorderRoutingMethod>(60, 80);
    RoutingMethodPtr brmp = std::make_shared<BoxDecompositionRoutingMethod>();
    std::vector<RoutingMethodPtr> mrcon = {mrmp, rmp, brmp};
    Circuit circ = CircuitsForTesting::get().uccsd;
    CompilationUnit cu{circ};
    PassPtr placement = gen_placement_pass(place);
    placement->apply(cu);
    CompilationUnit copy = cu;
    PassPtr pp = gen_routing_pass(arc, mrcon);
    nlohmann::json j_pp = pp;
    PassPtr loaded = j_pp.get<PassPtr>();
    pp->apply(cu);
    loaded->apply(copy);
    REQUIRE(cu.get_circ_ref() == copy.get_circ_ref());
    nlohmann::json j_loaded = loaded;
    REQUIRE(j_pp == j_loaded);
  }
  GIVEN("FullMappingPass") {
    // Sequence pass - deserializable only
    Circuit circ = CircuitsForTesting::get().uccsd;
    CompilationUnit cu{circ};
    CompilationUnit copy = cu;
    PassPtr pp = gen_full_mapping_pass(arc, place, rcon);
    nlohmann::json j_pp;
    j_pp["pass_class"] = "StandardPass";
    j_pp["StandardPass"]["name"] = "FullMappingPass";
    j_pp["StandardPass"]["architecture"] = arc;
    j_pp["StandardPass"]["placement"] = place;

    nlohmann::json config_array;
    for (const auto& con : rcon) {
      config_array.push_back(*con);
    }

    j_pp["StandardPass"]["routing_config"] = config_array;
    PassPtr loaded = j_pp.get<PassPtr>();
    pp->apply(cu);
    loaded->apply(copy);
    REQUIRE(cu.get_circ_ref() == copy.get_circ_ref());
  }
  GIVEN("DefaultMappingPass") {
    // Sequence pass - deserializable only
    Circuit circ = CircuitsForTesting::get().uccsd;
    CompilationUnit cu{circ};
    CompilationUnit copy = cu;
    PassPtr pp = gen_default_mapping_pass(arc, true);
    nlohmann::json j_pp;
    j_pp["pass_class"] = "StandardPass";
    j_pp["StandardPass"]["name"] = "DefaultMappingPass";
    j_pp["StandardPass"]["architecture"] = arc;
    j_pp["StandardPass"]["delay_measures"] = true;
    PassPtr loaded = j_pp.get<PassPtr>();
    pp->apply(cu);
    loaded->apply(copy);
    REQUIRE(cu.get_circ_ref() == copy.get_circ_ref());
  }
  GIVEN("CXMappingPass") {
    // Sequence pass - deserializable only
    Circuit circ = CircuitsForTesting::get().uccsd;
    CompilationUnit cu{circ};
    CompilationUnit copy = cu;
    PassPtr pp = gen_cx_mapping_pass(arc, place, rcon, true, false);
    nlohmann::json j_pp;
    j_pp["pass_class"] = "StandardPass";
    j_pp["StandardPass"]["name"] = "CXMappingPass";
    j_pp["StandardPass"]["architecture"] = arc;
    j_pp["StandardPass"]["placement"] = place;
    nlohmann::json config_array;
    for (const auto& con : rcon) {
      config_array.push_back(*con);
    }
    j_pp["StandardPass"]["routing_config"] = config_array;
    j_pp["StandardPass"]["directed"] = true;
    j_pp["StandardPass"]["delay_measures"] = false;
    PassPtr loaded = j_pp.get<PassPtr>();
    pp->apply(cu);
    loaded->apply(copy);
    REQUIRE(cu.get_circ_ref() == copy.get_circ_ref());
  }
  GIVEN("PauliSquash") {
    // Sequence pass - deserializable only
    Circuit circ = CircuitsForTesting::get().uccsd;
    CompilationUnit cu{circ};
    CompilationUnit copy = cu;
    PassPtr pp =
        PauliSquash(Transforms::PauliSynthStrat::Sets, CXConfigType::Star);
    nlohmann::json j_pp;
    j_pp["pass_class"] = "StandardPass";
    j_pp["StandardPass"]["name"] = "PauliSquash";
    j_pp["StandardPass"]["pauli_synth_strat"] =
        Transforms::PauliSynthStrat::Sets;
    j_pp["StandardPass"]["cx_config"] = CXConfigType::Star;
    PassPtr loaded = j_pp.get<PassPtr>();
    pp->apply(cu);
    loaded->apply(copy);
    REQUIRE(cu.get_circ_ref() == copy.get_circ_ref());
  }
  GIVEN("ContextSimp") {
    // Sequence pass - deserializable only
    Circuit circ = CircuitsForTesting::get().uccsd;
    circ.qubit_create_all();
    CompilationUnit cu{circ};
    CompilationUnit copy = cu;
    PassPtr pp = gen_contextual_pass(
        Transforms::AllowClassical::Yes,
        std::make_shared<Circuit>(CircPool::X()));
    nlohmann::json j_pp;
    j_pp["pass_class"] = "StandardPass";
    j_pp["StandardPass"]["name"] = "ContextSimp";
    j_pp["StandardPass"]["allow_classical"] = true;
    j_pp["StandardPass"]["x_circuit"] = CircPool::X();
    PassPtr loaded = j_pp.get<PassPtr>();
    pp->apply(cu);
    loaded->apply(copy);
    REQUIRE(cu.get_circ_ref() == copy.get_circ_ref());
  }
}

SCENARIO("Test compiler pass combinator serializations") {
  GIVEN("Sequences of passes") {
    Circuit circ = CircuitsForTesting::get().uccsd;
    CompilationUnit cu{circ};
    CompilationUnit copy = cu;
    std::vector<PassPtr> seq_vec = {
        gen_synthesise_pauli_graph(), gen_clifford_simp_pass()};
    PassPtr seq = std::make_shared<SequencePass>(seq_vec);
    nlohmann::json j_seq = seq;
    PassPtr loaded_seq = j_seq.get<PassPtr>();
    seq->apply(cu);
    loaded_seq->apply(copy);
    REQUIRE(cu.get_circ_ref() == copy.get_circ_ref());
    nlohmann::json j_loaded_seq = loaded_seq;
    REQUIRE(j_seq == j_loaded_seq);
  }
  GIVEN("A complex pass with multiple combinators") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::Z, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::Z, {0});
    circ.add_op<unsigned>(OpType::Z, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    CompilationUnit cu{circ};
    CompilationUnit copy = cu;
    PredicatePtr gate_set =
        std::make_shared<GateSetPredicate>(OpTypeSet{OpType::Z});
    PassPtr seq = std::make_shared<SequencePass>(
        std::vector<PassPtr>{RemoveRedundancies(), CommuteThroughMultis()});
    PassPtr rep = std::make_shared<RepeatUntilSatisfiedPass>(seq, gate_set);
    PassPtr comb =
        std::make_shared<SequencePass>(std::vector<PassPtr>{rep, RebaseTket()});
    nlohmann::json j_comb = comb;
    PassPtr loaded_comb = j_comb.get<PassPtr>();
    comb->apply(cu);
    loaded_comb->apply(copy);
    REQUIRE(cu.get_circ_ref() == copy.get_circ_ref());
    nlohmann::json j_loaded_comb = loaded_comb;
    REQUIRE(j_comb == j_loaded_comb);
  }
}

SCENARIO("Test QubitPauliString serialization") {
  QubitPauliString qps(
      {{Qubit(2), Pauli::X}, {Qubit(7), Pauli::Y}, {Qubit(0), Pauli::I}});

  nlohmann::json j_qps = qps;
  QubitPauliString new_qps = j_qps.get<QubitPauliString>();

  REQUIRE(qps == new_qps);
}

SCENARIO("Test MeasurementSetup serializations") {
  GIVEN("MeasurementBitMap") {
    MeasurementSetup::MeasurementBitMap map(0, {0, 1}, 1);
    nlohmann::json j_map = map;
    nlohmann::json j_correct_map = {
        {"circ_index", 0}, {"bits", {0, 1}}, {"invert", true}};
    REQUIRE(j_map == j_correct_map);
    MeasurementSetup::MeasurementBitMap map_loaded =
        j_map.get<MeasurementSetup::MeasurementBitMap>();
    nlohmann::json j_loaded_map = map_loaded;
    REQUIRE(j_loaded_map == j_map);
  }
  GIVEN("MeasurementBitMap with default constructor") {
    MeasurementSetup::MeasurementBitMap map;
    nlohmann::json j_map = map;
    nlohmann::json j_correct_map = {
        {"circ_index", 0},
        {"bits", nlohmann::json::array()},
        {"invert", false}};
    REQUIRE(j_map == j_correct_map);
    MeasurementSetup::MeasurementBitMap map_loaded =
        j_map.get<MeasurementSetup::MeasurementBitMap>();
    nlohmann::json j_loaded_map = map_loaded;
    REQUIRE(j_loaded_map == j_map);
  }
  GIVEN("MeasurementSetup") {
    MeasurementSetup ms;
    Circuit mc(2, 2);
    mc.add_measure(0, 0);
    mc.add_measure(1, 1);
    Circuit mc2(2, 2);
    mc2.add_measure(0, 1);
    mc2.add_measure(1, 0);
    ms.add_measurement_circuit(mc);
    ms.add_measurement_circuit(mc2);
    Qubit q0(q_default_reg(), 0);
    Qubit q1(q_default_reg(), 1);
    QubitPauliString ii;
    QubitPauliString zi({{q0, Pauli::Z}});
    QubitPauliString iz({{q1, Pauli::Z}});
    QubitPauliString zz({{q0, Pauli::Z}, {q1, Pauli::Z}});
    QubitPauliString xx({{q0, Pauli::X}, {q1, Pauli::X}});
    QubitPauliString yy({{q0, Pauli::Y}, {q1, Pauli::Y}});
    ms.add_result_for_term(ii, {0, {}, false});
    ms.add_result_for_term(zi, {0, {0}, false});
    ms.add_result_for_term(iz, {0, {1}, false});
    ms.add_result_for_term(zz, {0, {0, 1}, false});
    ms.add_result_for_term(zi, {1, {0}, true});
    ms.add_result_for_term(xx, {1, {0, 1}, true});
    ms.add_result_for_term(yy, {1, {0, 1}, true});
    nlohmann::json j_ms = ms;
    nlohmann::json j_circs = {mc, mc2};
    nlohmann::json j_result_map = {
        {ii,
         {{{"circ_index", 0},
           {"bits", nlohmann::json::array()},
           {"invert", false}}}},
        {iz, {{{"circ_index", 0}, {"bits", {1}}, {"invert", false}}}},
        {xx, {{{"circ_index", 1}, {"bits", {0, 1}}, {"invert", true}}}},
        {yy, {{{"circ_index", 1}, {"bits", {0, 1}}, {"invert", true}}}},
        {zi,
         {{{"circ_index", 0}, {"bits", {0}}, {"invert", false}},
          {{"circ_index", 1}, {"bits", {0}}, {"invert", true}}}},
        {zz, {{{"circ_index", 0}, {"bits", {0, 1}}, {"invert", false}}}},
    };
    REQUIRE(j_ms["circs"] == j_circs);
    REQUIRE(j_ms["result_map"] == j_result_map);
    MeasurementSetup ms_loaded = j_ms.get<MeasurementSetup>();
    nlohmann::json j_loaded_ms = ms_loaded;
    REQUIRE(j_loaded_ms == j_ms);
  }
  GIVEN("Empty MeasurementSetup") {
    MeasurementSetup ms;
    nlohmann::json j_ms = ms;
    nlohmann::json j_correct_ms = {
        {"circs", nlohmann::json::array()},
        {"result_map", nlohmann::json::array()}};
    REQUIRE(j_ms == j_correct_ms);
    MeasurementSetup ms_loaded = j_ms.get<MeasurementSetup>();
    nlohmann::json j_loaded_ms = ms_loaded;
    REQUIRE(j_loaded_ms == j_ms);
  }
}

}  // namespace test_json
}  // namespace tket
