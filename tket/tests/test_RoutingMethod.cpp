#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "Mapping/MappingFrontier.hpp"
#include "Mapping/RoutingMethodCircuit.hpp"
#include "Placement/Placement.hpp"

namespace tket {

SCENARIO("Test RoutingMethod default methods.") {
  RoutingMethod rm;
  Architecture arc(
      {{Node("t", 1), Node("t", 0)}, {Node("t", 2), Node("t", 1)}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
  Circuit circ(3);
  MappingFrontier_ptr mf = std::make_shared<MappingFrontier>(circ);
  unit_map_t empty;
  std::pair<bool, unit_map_t> rm_return = rm.routing_method(mf, shared_arc);
  REQUIRE(!rm_return.first);
  REQUIRE(rm_return.second == empty);
}

std::tuple<bool, Circuit, unit_map_t, unit_map_t>
test_routing_method_mf_simple_relabel(
    const Circuit& c, const ArchitecturePtr& a) {
  Circuit copy(c);
  std::vector<Qubit> qs = copy.all_qubits();
  std::vector<Node> ns = a->get_all_nodes_vec();
  // enforce in tests that ns >= qs, this is testing purposes only so fine...
  unit_map_t rename_map, final_map;
  for (unsigned i = 0; i < qs.size(); i++) {
    rename_map.insert({qs[i], ns[i]});
    final_map.insert({ns[i], ns[i]});
  }
  copy.rename_units(rename_map);
  return {true, copy, rename_map, final_map};
}

std::tuple<bool, Circuit, unit_map_t, unit_map_t>
test_routing_method_mf_swap_perm(const Circuit& c, const ArchitecturePtr& a) {
  if (c.n_qubits() > 2 && a->n_nodes() > 2) {
    Circuit copy(c);
    std::vector<Qubit> qs = copy.all_qubits();
    std::vector<Node> ns = a->get_all_nodes_vec();
    // enforce in tests that ns >= qs, this is testing purposes only so fine...
    unit_map_t rename_map, final_map;
    for (unsigned i = 0; i < qs.size(); i++) {
      rename_map.insert({qs[i], ns[i]});
      final_map.insert({ns[i], ns[i]});
    }
    copy.rename_units(rename_map);
    MappingFrontier mf(copy);
    // n.b. add_swap permutes out edge of both boundaries,
    mf.add_swap(Node("t", 0), Node("t", 1));

    return {true, copy, rename_map, final_map};
  } else {
    return {false, Circuit(), {}, {}};
  }
}

std::tuple<bool, Circuit, unit_map_t, unit_map_t>
test_routing_method_mf_swap_no_perm(
    const Circuit& c, const ArchitecturePtr& a) {
  if (c.n_qubits() > 2 && a->n_nodes() > 2) {
    Circuit copy(c);
    std::vector<Qubit> qs = copy.all_qubits();
    std::vector<Node> ns = a->get_all_nodes_vec();
    // enforce in tests that ns >= qs, this is testing purposes only so fine...
    unit_map_t rename_map, final_map;
    for (unsigned i = 0; i < qs.size(); i++) {
      rename_map.insert({qs[i], ns[i]});
      final_map.insert({ns[i], ns[i]});
    }
    copy.rename_units(rename_map);
    MappingFrontier mf(copy);
    // n.b. add_swap permutes out edge of both boundaries,
    mf.add_swap(Node("t", 0), Node("t", 1));
    final_map[Node("t", 0)] = Node("t", 1);
    final_map[Node("t", 1)] = Node("t", 0);

    return {true, copy, rename_map, final_map};
  } else {
    return {false, Circuit(), {}, {}};
  }
}

std::tuple<bool, Circuit, unit_map_t, unit_map_t>
test_routing_method_circuit_no_perm(
    const Circuit& c, const ArchitecturePtr& a) {
  if (c.n_qubits() > 2 && a->n_nodes() > 2) {
    Circuit copy(c.n_qubits());
    copy.add_op<unsigned>(OpType::SWAP, {0, 1});
    copy.add_op<unsigned>(OpType::CX, {1, 0});
    copy.add_op<unsigned>(OpType::CX, {1, 0});

    std::vector<Qubit> qs = copy.all_qubits();
    std::vector<Node> ns = a->get_all_nodes_vec();
    // enforce in tests that ns >= qs, this is testing purposes only so fine...
    unit_map_t rename_map, final_map;
    for (unsigned i = 0; i < qs.size(); i++) {
      rename_map.insert({qs[i], ns[i]});
      final_map.insert({ns[i], ns[i]});
    }
    copy.rename_units(rename_map);
    MappingFrontier mf(copy);
    final_map[Node("t", 0)] = Node("t", 1);
    final_map[Node("t", 1)] = Node("t", 0);
    return {true, copy, rename_map, final_map};
  } else {
    return {false, Circuit(), {}, {}};
  }
}

SCENARIO("Test RoutingMethodCircuit checking criteria") {
  RoutingMethodCircuit rmc(test_routing_method_mf_swap_no_perm, 5, 5);
  Circuit c(2), circ3(3);
  c.add_op<unsigned>(OpType::CX, {0, 1});
  circ3.add_op<unsigned>(OpType::CX, {0, 2});
  circ3.add_op<unsigned>(OpType::CX, {2, 1});
  MappingFrontier_ptr mf2 = std::make_shared<MappingFrontier>(c);
  MappingFrontier_ptr mf3 = std::make_shared<MappingFrontier>(circ3);

  Architecture arc(
      {{Node("t", 1), Node("t", 0)}, {Node("t", 2), Node("t", 1)}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);

  std::pair<bool, unit_map_t> res0 = rmc.routing_method(mf2, shared_arc);
  REQUIRE(!res0.first);
  std::pair<bool, unit_map_t> res1 = rmc.routing_method(mf3, shared_arc);
  REQUIRE(res1.first);
}
SCENARIO("Test RoutingMethodCircuit::routing_method") {
  Circuit comp(3);
  comp.add_op<unsigned>(OpType::SWAP, {0, 1});
  comp.add_op<unsigned>(OpType::CX, {1, 0});
  comp.add_op<unsigned>(OpType::CX, {1, 0});
  comp.add_op<unsigned>(OpType::CX, {1, 0});
  comp.add_op<unsigned>(OpType::CX, {1, 0});
  auto qbs = comp.all_qubits();
  unit_map_t rename_map = {
      {qbs[0], Node("t", 0)}, {qbs[1], Node("t", 1)}, {qbs[2], Node("t", 2)}};
  comp.rename_units(rename_map);
  qubit_map_t permutation = {
      {Node("t", 0), Node("t", 1)}, {Node("t", 1), Node("t", 0)}};
  comp.permute_boundary_output(permutation);

  GIVEN("Non-implicit Permutation method, using MappingFrontier::add_swap") {
    RoutingMethodCircuit rmc(test_routing_method_mf_swap_no_perm, 2, 2);
    Circuit c(3);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});

    MappingFrontier_ptr mf = std::make_shared<MappingFrontier>(c);
    Architecture arc(
        {{Node("t", 1), Node("t", 0)}, {Node("t", 2), Node("t", 1)}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
    std::pair<bool, unit_map_t> output = rmc.routing_method(mf, shared_arc);
    unit_map_t empty;
    REQUIRE(output.first);
    REQUIRE(output.second == empty);
    REQUIRE(c == comp);
  }
  GIVEN("Non-implicit Permutation method, using circuit replacement") {
    RoutingMethodCircuit rmc(test_routing_method_circuit_no_perm, 2, 2);
    Circuit c(3);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});

    MappingFrontier_ptr mf = std::make_shared<MappingFrontier>(c);
    Architecture arc(
        {{Node("t", 1), Node("t", 0)}, {Node("t", 2), Node("t", 1)}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
    std::pair<bool, unit_map_t> output = rmc.routing_method(mf, shared_arc);
    unit_map_t empty;
    REQUIRE(output.first);
    REQUIRE(output.second == empty);
    REQUIRE(c == comp);
  }
  GIVEN("Implicit Permutation method, using MappingFrontier::add_swap") {
    RoutingMethodCircuit rmc(test_routing_method_mf_swap_perm, 2, 2);
    Circuit c(3);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});

    MappingFrontier_ptr mf = std::make_shared<MappingFrontier>(c);
    Architecture arc(
        {{Node("t", 1), Node("t", 0)}, {Node("t", 2), Node("t", 1)}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
    std::pair<bool, unit_map_t> output = rmc.routing_method(mf, shared_arc);
    unit_map_t empty;
    REQUIRE(output.first);
    REQUIRE(output.second == empty);

    Circuit comp1(3);
    comp1.add_op<unsigned>(OpType::SWAP, {0, 1});
    comp1.add_op<unsigned>(OpType::CX, {1, 0});
    comp1.add_op<unsigned>(OpType::CX, {1, 0});
    comp1.add_op<unsigned>(OpType::CX, {0, 1});
    comp1.add_op<unsigned>(OpType::CX, {0, 1});
    qbs = comp1.all_qubits();
    rename_map = {
        {qbs[0], Node("t", 0)}, {qbs[1], Node("t", 1)}, {qbs[2], Node("t", 2)}};
    comp1.rename_units(rename_map);

    REQUIRE(c == comp1);
  }
}

SCENARIO("Test RoutingMethodCircuit produces correct map") {
  RoutingMethodCircuit rmc(test_routing_method_mf_simple_relabel, 5, 5);
  Architecture arc({{Node(0), Node(1)}, {Node(1), Node(2)}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
  Circuit c(3);
  c.add_op<unsigned>(OpType::CX, {0, 1});
  c.add_op<unsigned>(OpType::CX, {1, 2});
  std::shared_ptr<unit_bimaps_t> maps = std::make_shared<unit_bimaps_t>();
  // Initialise the maps by the same way it's done with CompilationUnit
  for (const UnitID& u : c.all_units()) {
    maps->initial.insert({u, u});
    maps->final.insert({u, u});
  }
  Placement pl(arc);
  qubit_mapping_t partial_map;
  partial_map.insert({Qubit(0), Node(0)});
  partial_map.insert({Qubit(1), Node(1)});
  // We leave q[2] unplaced
  pl.place_with_map(c, partial_map, maps);
  MappingFrontier_ptr mf = std::make_shared<MappingFrontier>(c, maps);

  std::pair<bool, unit_map_t> res = rmc.routing_method(mf, shared_arc);

  for (const Qubit& q : c.all_qubits()) {
    REQUIRE(maps->initial.right.find(q) != maps->initial.right.end());
    REQUIRE(maps->final.right.find(q) != maps->final.right.end());
  }
}

}  // namespace tket
