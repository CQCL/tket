#include <algorithm>
#include <catch2/catch_test_macros.hpp>

#include "Circuit/Circuit.hpp"
#include "Mapping/AASLabelling.hpp"
#include "Mapping/AASRoute.hpp"
#include "Mapping/LexiLabelling.hpp"
#include "Mapping/LexiRoute.hpp"
#include "Mapping/MappingManager.hpp"
#include "OpType/OpType.hpp"
#include "OpType/OpTypeFunctions.hpp"
#include "Predicates/CompilationUnit.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Predicates/PassLibrary.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Transformations/ContextualReduction.hpp"
#include "testutil.hpp"

namespace tket {
SCENARIO("Test aas route in RV3") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3), Node("node_test", 4), Node("node_test", 5),
      Node("test_node", 6), Node("node_test", 7), Node("node_test", 8),
      Node("node_test", 9), Node("node_test", 10)};

  Architecture architecture(
      {{nodes[0], nodes[1]},
       {nodes[1], nodes[2]},
       {nodes[2], nodes[3]},
       {nodes[3], nodes[4]},
       {nodes[2], nodes[5]},
       {nodes[5], nodes[6]},
       {nodes[4], nodes[7]},
       {nodes[7], nodes[8]},
       {nodes[8], nodes[9]},
       {nodes[9], nodes[10]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);

  GIVEN(
      "AASRoute - test AASRouteRoutingMethod routing_method placed and gates") {
    Circuit circ(11);
    std::vector<Qubit> qubits = circ.all_qubits();

    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});

    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
    circ.add_op<UnitID>(OpType::CX, {qubits[10], qubits[0]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[0]});

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]},  {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]},  {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]},  {qubits[8], nodes[8]},
        {qubits[9], nodes[9]}, {qubits[10], nodes[10]}};

    circ.rename_units(rename_map);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    AASRouteRoutingMethod aasrm(1, aas::CNotSynthType::Rec);

    // this will fail because the cx ae in fron of the ppb
    REQUIRE(!aasrm.routing_method(mf, shared_arc).first);
  }
  GIVEN("AASRoute - test AASRouteRoutingMethod routing_method placed") {
    Circuit circ(11);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]},  {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]},  {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]},  {qubits[8], nodes[8]},
        {qubits[9], nodes[9]}, {qubits[10], nodes[10]}};

    Circuit ppb_circ(11);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    circ.rename_units(rename_map);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    AASRouteRoutingMethod aasrm(1, aas::CNotSynthType::Rec);

    REQUIRE(aasrm.routing_method(mf, shared_arc).first);
  }
  GIVEN("AASRoute - test AASRouteRoutingMethod routing_method unplaced") {
    Circuit circ(11);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]},  {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]},  {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]},  {qubits[8], nodes[8]},
        {qubits[9], nodes[9]}, {qubits[10], nodes[10]}};

    Circuit ppb_circ(11);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    AASRouteRoutingMethod aasrm(1, aas::CNotSynthType::Rec);

    // this will fail because of the unplaced qubits
    REQUIRE(!aasrm.routing_method(mf, shared_arc).first);
  }
  GIVEN("AASRouteRoutingMethod - test routing_method I") {
    std::vector<Node> nodes_mixed = {
        Node("node_test", 0), Node("test_node", 1), Node("node_test", 2)};

    Architecture architecture_mixed(
        {{nodes_mixed[0], nodes_mixed[1]}, {nodes_mixed[1], nodes_mixed[2]}});
    ArchitecturePtr shared_arc_mixed =
        std::make_shared<Architecture>(architecture_mixed);

    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes_mixed[0]},
        {qubits[1], nodes_mixed[1]},
        {qubits[2], nodes_mixed[2]}};

    Circuit ppb_circ(3);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});

    PhasePolyBox ppbox(ppb_circ);

    circ.add_box(ppbox, qubits);
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});

    // testing this without interacting with the lexilabelling or aas labelling
    circ.rename_units(rename_map);

    /*std::vector<Qubit> qubits_after_rename = circ.all_qubits();

    for (const Command& com : circ) {
      OpType ot = com.get_op_ptr()->get_type();
      unit_vector_t qbs = com.get_args();
      switch (ot) {
        case OpType::PhasePolyBox: {
          Op_ptr op = com.get_op_ptr();
          const PhasePolyBox& ppb = static_cast<const PhasePolyBox&>(*op);
          Circuit circuit_ppb_place(*ppb.to_circuit());
        }
      }
    }

    Circuit circ_flatt_copy(circ);

    circ_flatt_copy.flatten_registers();

    std::vector<Qubit> qubits_after_flatten = circ_flatt_copy.all_qubits();

    for (const Command& com : circ_flatt_copy) {
      OpType ot = com.get_op_ptr()->get_type();
      unit_vector_t qbs = com.get_args();
      switch (ot) {
        case OpType::PhasePolyBox: {
          Op_ptr op = com.get_op_ptr();
          const PhasePolyBox& ppb = static_cast<const PhasePolyBox&>(*op);
          Circuit circuit_ppb_place(*ppb.to_circuit());
        }
      }
    }*/

    Circuit circ_copy(circ);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    AASRouteRoutingMethod aasrm(1, aas::CNotSynthType::Rec);

    REQUIRE(aasrm.routing_method(mf, shared_arc_mixed).first);
    // aasrm.routing_method(mf, shared_arc_mixed);

    REQUIRE(test_unitary_comparison(mf->circuit_, circ_copy));

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture_mixed);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu0(mf->circuit_, preds);
    REQUIRE(cu0.check_all_predicates());
  }
  GIVEN("AASRouteRoutingMethod - test routing_method II") {
    std::vector<Node> nodes_mixed = {
        Node("node_test", 0), Node("test_node", 1), Node("node_test", 2)};

    Architecture architecture_mixed(
        {{nodes_mixed[0], nodes_mixed[1]}, {nodes_mixed[1], nodes_mixed[2]}});

    ArchitecturePtr shared_arc_mixed =
        std::make_shared<Architecture>(architecture_mixed);

    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes_mixed[0]},
        {qubits[1], nodes_mixed[1]},
        {qubits[2], nodes_mixed[2]}};

    Circuit ppb_circ(3);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});

    ppb_circ.add_op<UnitID>(OpType::Rz, 0.22, {qubits[0]});
    ppb_circ.add_op<UnitID>(OpType::Rz, 0.33, {qubits[1]});
    ppb_circ.add_op<UnitID>(OpType::Rz, 0.55, {qubits[2]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});

    circ.rename_units(rename_map);

    /*std::vector<Qubit> qubits_after_rename = circ.all_qubits();

    for (const Command& com : circ) {
      OpType ot = com.get_op_ptr()->get_type();
      unit_vector_t qbs = com.get_args();
      switch (ot) {
        case OpType::PhasePolyBox: {
          Op_ptr op = com.get_op_ptr();
          const PhasePolyBox& ppb = static_cast<const PhasePolyBox&>(*op);
          Circuit circuit_ppb_place(*ppb.to_circuit());
        }
      }
    }

    Circuit circ_flatt_copy(circ);

    circ_flatt_copy.flatten_registers();

    std::vector<Qubit> qubits_after_flatten = circ_flatt_copy.all_qubits();

    for (const Command& com : circ_flatt_copy) {
      OpType ot = com.get_op_ptr()->get_type();
      unit_vector_t qbs = com.get_args();
      switch (ot) {
        case OpType::PhasePolyBox: {
          Op_ptr op = com.get_op_ptr();
          const PhasePolyBox& ppb = static_cast<const PhasePolyBox&>(*op);
          Circuit circuit_ppb_place(*ppb.to_circuit());
        }
      }
    }*/

    Circuit circ_copy(circ);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    AASRouteRoutingMethod aasrm(1, aas::CNotSynthType::Rec);

    REQUIRE(aasrm.routing_method(mf, shared_arc_mixed).first);
    // aasrm.routing_method(mf, shared_arc_mixed);

    REQUIRE(test_unitary_comparison(mf->circuit_, circ_copy));

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture_mixed);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu0(mf->circuit_, preds);
    REQUIRE(cu0.check_all_predicates());
  }
  GIVEN("AASRouteRoutingMethod  and LexiRouteRoutingMethod I") {
    std::vector<Node> nodes_mixed = {
        Node("node_test", 0), Node("test_node", 1), Node("node_test", 2)};

    Architecture architecture_mixed(
        {{nodes_mixed[0], nodes_mixed[1]}, {nodes_mixed[1], nodes_mixed[2]}});

    ArchitecturePtr shared_arc_mixed =
        std::make_shared<Architecture>(architecture_mixed);

    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes_mixed[0]},
        {qubits[1], nodes_mixed[1]},
        {qubits[2], nodes_mixed[2]}};

    Circuit ppb_circ(3);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    circ.rename_units(rename_map);

    Circuit circ_copy(circ);

    MappingManager mm(shared_arc_mixed);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<AASRouteRoutingMethod>(1, aas::CNotSynthType::Rec),
        std::make_shared<LexiRouteRoutingMethod>(100),
    };

    mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture_mixed);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);
    REQUIRE(cu.check_all_predicates());
    REQUIRE(test_unitary_comparison(circ, circ_copy));
    REQUIRE(circ.n_gates() == 4);
  }

  GIVEN("AASRouteRoutingMethod  and LexiRouteRoutingMethod II") {
    std::vector<Node> nodes_mixed = {
        Node("node_test", 0), Node("test_node", 1), Node("node_test", 2)};

    Architecture architecture_mixed(
        {{nodes_mixed[0], nodes_mixed[1]}, {nodes_mixed[1], nodes_mixed[2]}});

    ArchitecturePtr shared_arc_mixed =
        std::make_shared<Architecture>(architecture_mixed);

    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes_mixed[0]},
        {qubits[1], nodes_mixed[1]},
        {qubits[2], nodes_mixed[2]}};

    Circuit ppb_circ(3);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});

    circ.rename_units(rename_map);

    Circuit circ_copy(circ);

    MappingManager mm(shared_arc_mixed);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<AASRouteRoutingMethod>(1),
        std::make_shared<LexiRouteRoutingMethod>(100),
    };

    mm.route_circuit(circ, vrm);
  }
  GIVEN("AASRouteRoutingMethod  and LexiRouteRoutingMethod III") {
    std::vector<Node> nodes_mixed = {
        Node("node_test", 0), Node("test_node", 1), Node("node_test", 2)};

    Architecture architecture_mixed(
        {{nodes_mixed[0], nodes_mixed[1]}, {nodes_mixed[1], nodes_mixed[2]}});

    ArchitecturePtr shared_arc_mixed =
        std::make_shared<Architecture>(architecture_mixed);

    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes_mixed[0]},
        {qubits[1], nodes_mixed[1]},
        {qubits[2], nodes_mixed[2]}};

    Circuit ppb_circ(3);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    Circuit circ_copy(circ);

    MappingManager mm(shared_arc_mixed);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<AASRouteRoutingMethod>(1, aas::CNotSynthType::Rec),
        std::make_shared<AASLabellingMethod>(),
    };

    mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture_mixed);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);
    REQUIRE(cu.check_all_predicates());
    REQUIRE(test_unitary_comparison(circ, circ_copy));
    REQUIRE(circ.n_gates() == 1);
  }
  GIVEN("AASRouteRoutingMethod  and LexiRouteRoutingMethod IV") {
    std::vector<Node> nodes_mixed = {
        Node("node_test", 0), Node("test_node", 1), Node("node_test", 2)};

    Architecture architecture_mixed(
        {{nodes_mixed[0], nodes_mixed[1]}, {nodes_mixed[1], nodes_mixed[2]}});

    ArchitecturePtr shared_arc_mixed =
        std::make_shared<Architecture>(architecture_mixed);

    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes_mixed[0]},
        {qubits[1], nodes_mixed[1]},
        {qubits[2], nodes_mixed[2]}};

    Circuit ppb_circ(3);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    Circuit circ_copy(circ);

    MappingManager mm(shared_arc_mixed);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<AASRouteRoutingMethod>(1, aas::CNotSynthType::Rec),
        std::make_shared<AASLabellingMethod>(),
    };

    mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture_mixed);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);
    REQUIRE(cu.check_all_predicates());
    REQUIRE(test_unitary_comparison(circ, circ_copy));
    REQUIRE(circ.n_gates() == 4);
  }
  GIVEN("AASRouteRoutingMethod  and LexiRouteRoutingMethod V") {
    std::vector<Node> nodes_mixed = {
        Node("node_test", 0), Node("test_node", 1), Node("node_test", 2)};

    Architecture architecture_mixed(
        {{nodes_mixed[0], nodes_mixed[1]}, {nodes_mixed[1], nodes_mixed[2]}});

    ArchitecturePtr shared_arc_mixed =
        std::make_shared<Architecture>(architecture_mixed);

    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes_mixed[0]},
        {qubits[1], nodes_mixed[1]},
        {qubits[2], nodes_mixed[2]}};

    Circuit ppb_circ(3);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    circ.rename_units(rename_map);

    Circuit circ_copy(circ);

    MappingManager mm(shared_arc_mixed);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<AASRouteRoutingMethod>(1, aas::CNotSynthType::Rec),
        std::make_shared<LexiLabellingMethod>(),
        std::make_shared<LexiRouteRoutingMethod>(100),
    };

    mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture_mixed);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);
    REQUIRE(cu.check_all_predicates());
    REQUIRE(test_unitary_comparison(circ, circ_copy));
    REQUIRE(circ.n_gates() == 4);
  }
  GIVEN("AASRouteRoutingMethod  and LexiRouteRoutingMethod VI") {
    std::vector<Node> nodes_mixed = {
        Node("node_test", 0), Node("test_node", 1), Node("node_test", 2)};

    Architecture architecture_mixed(
        {{nodes_mixed[0], nodes_mixed[1]}, {nodes_mixed[1], nodes_mixed[2]}});

    ArchitecturePtr shared_arc_mixed =
        std::make_shared<Architecture>(architecture_mixed);

    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes_mixed[0]},
        {qubits[1], nodes_mixed[1]},
        {qubits[2], nodes_mixed[2]}};

    Circuit ppb_circ(3);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    circ.rename_units(rename_map);

    Circuit circ_copy(circ);

    MappingManager mm(shared_arc_mixed);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<AASRouteRoutingMethod>(1, aas::CNotSynthType::Rec),
        std::make_shared<LexiLabellingMethod>(),
        std::make_shared<LexiRouteRoutingMethod>(100),
        std::make_shared<AASLabellingMethod>(),
    };

    mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture_mixed);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);
    REQUIRE(cu.check_all_predicates());
    REQUIRE(test_unitary_comparison(circ, circ_copy));
    REQUIRE(circ.n_gates() == 4);
  }

  GIVEN("AASRouteRoutingMethod  and LexiRouteRoutingMethod VII") {
    Circuit circ(11);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]},  {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]},  {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]},  {qubits[8], nodes[8]},
        {qubits[9], nodes[9]}, {qubits[10], nodes[10]}};

    Circuit ppb_circ(11);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[5]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[5], qubits[4]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});

    circ.rename_units(rename_map);

    Circuit circ_copy(circ);

    MappingManager mm(shared_arc);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<AASLabellingMethod>(),
        std::make_shared<AASRouteRoutingMethod>(1, aas::CNotSynthType::Rec),
        std::make_shared<LexiLabellingMethod>(),
        std::make_shared<LexiRouteRoutingMethod>(100),
    };

    mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);
    REQUIRE(cu.check_all_predicates());

    // REQUIRE(test_unitary_comparison(circ, circ_copy));
    // will fail because of the inserted swaps
    REQUIRE(circ.n_gates() == 18);
    REQUIRE(circ.count_gates(OpType::CX) == 15);
    REQUIRE(circ.count_gates(OpType::SWAP) == 3);
  }
  GIVEN("AASRouteRoutingMethod  and LexiRouteRoutingMethod, only route") {
    Circuit circ(11);
    std::vector<Qubit> qubits = circ.all_qubits();

    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});

    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
    circ.add_op<UnitID>(OpType::CX, {qubits[10], qubits[0]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[0]});

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]},  {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]},  {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]},  {qubits[8], nodes[8]},
        {qubits[9], nodes[9]}, {qubits[10], nodes[10]}};

    Circuit ppb_circ(11);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[5]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[5], qubits[4]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});

    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
    circ.add_op<UnitID>(OpType::CX, {qubits[10], qubits[0]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[0]});

    circ.rename_units(rename_map);

    Circuit circ_copy(circ);

    MappingManager mm(shared_arc);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<AASLabellingMethod>(),
        std::make_shared<AASRouteRoutingMethod>(1, aas::CNotSynthType::Rec),
        std::make_shared<LexiLabellingMethod>(),
        std::make_shared<LexiRouteRoutingMethod>(100),
    };

    mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);
    REQUIRE(cu.check_all_predicates());

    REQUIRE(circ.n_gates() == 61);
    REQUIRE(circ.count_gates(OpType::CX) == 37);
    REQUIRE(circ.count_gates(OpType::SWAP) == 21);
  }
  GIVEN(
      "AASRouteRoutingMethod  and LexiRouteRoutingMethod only route two "
      "boxes") {
    Circuit circ(11);
    std::vector<Qubit> qubits = circ.all_qubits();

    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});

    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
    circ.add_op<UnitID>(OpType::CX, {qubits[10], qubits[0]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[0]});

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]},  {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]},  {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]},  {qubits[8], nodes[8]},
        {qubits[9], nodes[9]}, {qubits[10], nodes[10]}};

    Circuit ppb_circ(11);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[5]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[7], qubits[4]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[4]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[9], qubits[4]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});

    Circuit ppb_circ_2(11);

    ppb_circ_2.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    ppb_circ_2.add_op<UnitID>(OpType::CX, {qubits[2], qubits[5]});
    ppb_circ_2.add_op<UnitID>(OpType::CX, {qubits[5], qubits[4]});

    PhasePolyBox ppbox2(ppb_circ_2);
    circ.add_box(ppbox2, qubits);

    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
    circ.add_op<UnitID>(OpType::CX, {qubits[10], qubits[0]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[0]});

    circ.rename_units(rename_map);

    Circuit circ_copy(circ);

    MappingManager mm(shared_arc);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<AASLabellingMethod>(),
        std::make_shared<AASRouteRoutingMethod>(1, aas::CNotSynthType::Rec),
        std::make_shared<LexiLabellingMethod>(),
        std::make_shared<LexiRouteRoutingMethod>(100),
    };

    mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);
    REQUIRE(cu.check_all_predicates());

    REQUIRE(circ.n_gates() == 115);
    REQUIRE(circ.count_gates(OpType::CX) == 89);
    REQUIRE(circ.count_gates(OpType::SWAP) == 24);
  }
  GIVEN("AAS + Lexi, Label and Route") {
    Circuit circ(11);
    std::vector<Qubit> qubits = circ.all_qubits();

    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});

    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
    circ.add_op<UnitID>(OpType::CX, {qubits[10], qubits[0]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[0]});

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]},  {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]},  {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]},  {qubits[8], nodes[8]},
        {qubits[9], nodes[9]}, {qubits[10], nodes[10]}};

    Circuit ppb_circ(11);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[5]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[7], qubits[4]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[4]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[9], qubits[4]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});

    Circuit ppb_circ_2(11);

    ppb_circ_2.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    ppb_circ_2.add_op<UnitID>(OpType::CX, {qubits[2], qubits[5]});
    ppb_circ_2.add_op<UnitID>(OpType::CX, {qubits[5], qubits[4]});

    PhasePolyBox ppbox2(ppb_circ_2);
    circ.add_box(ppbox2, qubits);

    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
    circ.add_op<UnitID>(OpType::CX, {qubits[10], qubits[0]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[0]});

    Circuit circ_copy(circ);

    MappingManager mm(shared_arc);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<AASLabellingMethod>(),
        std::make_shared<AASRouteRoutingMethod>(1, aas::CNotSynthType::Rec),
        std::make_shared<LexiLabellingMethod>(),
        std::make_shared<LexiRouteRoutingMethod>(100),
    };

    mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);
    REQUIRE(cu.check_all_predicates());

    REQUIRE(circ.n_gates() == 275);
    REQUIRE(circ.count_gates(OpType::CX) == 244);
    REQUIRE(circ.count_gates(OpType::SWAP) == 28);
  }
}

SCENARIO("Test aas route in RV3 - long test", "[.long]") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3), Node("node_test", 4), Node("node_test", 5),
      Node("test_node", 6), Node("node_test", 7), Node("node_test", 8),
      Node("node_test", 9), Node("node_test", 10)};

  Architecture architecture(
      {{nodes[0], nodes[1]},
       {nodes[1], nodes[2]},
       {nodes[2], nodes[3]},
       {nodes[3], nodes[4]},
       {nodes[2], nodes[5]},
       {nodes[5], nodes[6]},
       {nodes[4], nodes[7]},
       {nodes[7], nodes[8]},
       {nodes[8], nodes[9]},
       {nodes[9], nodes[10]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  GIVEN("AASRouteRoutingMethod - test routing_method III") {
    Circuit circ(11);
    std::vector<Qubit> qubits = circ.all_qubits();

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]},  {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]},  {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]},  {qubits[8], nodes[8]},
        {qubits[9], nodes[9]}, {qubits[10], nodes[10]}};

    Circuit ppb_circ(11);

    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    ppb_circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});

    PhasePolyBox ppbox(ppb_circ);
    circ.add_box(ppbox, qubits);

    circ.rename_units(rename_map);

    Circuit circ_copy(circ);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    AASRouteRoutingMethod aasrm(1, aas::CNotSynthType::Rec);

    REQUIRE(aasrm.routing_method(mf, shared_arc).first);

    // aasrm.routing_method(mf, shared_arc);

    REQUIRE(test_unitary_comparison(mf->circuit_, circ_copy));

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu0(mf->circuit_, preds);
    REQUIRE(cu0.check_all_predicates());
  }
}

}  // namespace tket
