#include <algorithm>
#include <catch2/catch.hpp>

#include "Circuit/Circuit.hpp"
#include "Mapping/AASRoute.hpp"
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

  GIVEN("test check_method in AASRouteRoutingMethod I") {
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

    AASRouteRoutingMethod aasrm(aas::CNotSynthType::Rec, 1);

    REQUIRE(!aasrm.check_method(mf, shared_arc));
  }
  GIVEN("test check_method in AASRouteRoutingMethod II") {
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

    AASRouteRoutingMethod aasrm(aas::CNotSynthType::Rec, 1);

    REQUIRE(aasrm.check_method(mf, shared_arc));
  }
  GIVEN("test check_method in AASRouteRoutingMethod III") {
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

    AASRouteRoutingMethod aasrm(aas::CNotSynthType::Rec, 1);

    REQUIRE(!aasrm.check_method(mf, shared_arc));
  }
  GIVEN("test routing_method in AASRouteRoutingMethod I") {
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

    std::vector<Qubit> qubits_after_rename = circ.all_qubits();

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
    }

    Circuit circ_copy(circ);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    AASRouteRoutingMethod aasrm(aas::CNotSynthType::Rec, 1);

    REQUIRE(aasrm.check_method(mf, shared_arc_mixed));
    aasrm.routing_method(mf, shared_arc_mixed);

    REQUIRE(test_unitary_comparison(mf->circuit_, circ_copy));

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture_mixed);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu0(mf->circuit_, preds);
    REQUIRE(cu0.check_all_predicates());
  }
  GIVEN("test routing_method in AASRouteRoutingMethod II") {
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

    std::vector<Qubit> qubits_after_rename = circ.all_qubits();

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
    }

    Circuit circ_copy(circ);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    AASRouteRoutingMethod aasrm(aas::CNotSynthType::Rec, 1);

    REQUIRE(aasrm.check_method(mf, shared_arc_mixed));
    aasrm.routing_method(mf, shared_arc_mixed);

    REQUIRE(test_unitary_comparison(mf->circuit_, circ_copy));

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture_mixed);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu0(mf->circuit_, preds);
    REQUIRE(cu0.check_all_predicates());
  }
  GIVEN("test routing_method in AASRouteRoutingMethod III") {
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

    AASRouteRoutingMethod aasrm(aas::CNotSynthType::Rec, 1);

    REQUIRE(aasrm.check_method(mf, shared_arc));

    aasrm.routing_method(mf, shared_arc);

    REQUIRE(test_unitary_comparison(mf->circuit_, circ_copy));

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu0(mf->circuit_, preds);
    REQUIRE(cu0.check_all_predicates());
  }
  GIVEN("test routing_method in AASRouteRoutingMethod with lexi route I") {
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
        std::make_shared<AASRouteRoutingMethod>(aas::CNotSynthType::Rec, 1),
        std::make_shared<LexiRouteRoutingMethod>(100),
    };

    // this will not use the aas route, which is not clear to me
    mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture_mixed);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);
    REQUIRE(!cu.check_all_predicates());  // TODO

    REQUIRE(test_unitary_comparison(circ, circ_copy));

    REQUIRE(circ.n_gates() == 1);  // TODO
  }

  GIVEN("test routing_method in AASRouteRoutingMethod with lexi route II") {
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

    std::cout << "\n\nprint out circuit:\n";

    for (auto g : circ_copy) {
      std::cout << g << std::endl;
    }

    std::cout << "print end circuit\n\n";

    MappingManager mm(shared_arc_mixed);

    // AASRouteRoutingMethod aasrm(aas::CNotSynthType::Rec, 1);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<AASRouteRoutingMethod>(aas::CNotSynthType::Rec, 1),
        std::make_shared<LexiRouteRoutingMethod>(100),
    };

    mm.route_circuit(circ, vrm);
  }
}
}  // namespace tket