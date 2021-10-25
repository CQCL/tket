#include <catch2/catch.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "Mapping/MappingFrontier.hpp"
#include "Mapping/RoutingMethodCircuit.hpp"

namespace tket {

SCENARIO("Test RoutingMethod default methods.") {
  RoutingMethod rm;
  Architecture arc(
      {{Node("t", 1), Node("t", 0)}, {Node("t", 2), Node("t", 1)}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
  Circuit circ(3);
  std::shared_ptr<MappingFrontier> mf = std::make_shared<MappingFrontier>(circ);
  REQUIRE(!rm.check_method(mf, shared_arc));
  unit_map_t empty;
  REQUIRE(rm.routing_method(mf, shared_arc) == empty);
}

// These two method are not completely reflective of what is necessary for
// routing Their design is to minimally test the required features of the
// methods, not to actually succesfully route a circuit
bool test_check_method(const Circuit& c, const ArchitecturePtr& a) {
  if (c.n_qubits() > 2 && a->n_uids() > 2) {
    return true;
  } else {
    return false;
  }
}

std::tuple<Circuit, unit_map_t, unit_map_t> test_routing_method_mf_swap_perm(
    const Circuit& c, const ArchitecturePtr& a) {
  Circuit copy(c);
  std::vector<Qubit> qs = copy.all_qubits();
  std::vector<Node> ns = a->get_all_uids_vec();
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

  return std::make_tuple(copy, rename_map, final_map);
}

std::tuple<Circuit, unit_map_t, unit_map_t> test_routing_method_mf_swap_no_perm(
    const Circuit& c, const ArchitecturePtr& a) {
  Circuit copy(c);
  std::vector<Qubit> qs = copy.all_qubits();
  std::vector<Node> ns = a->get_all_uids_vec();
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

  return std::make_tuple(copy, rename_map, final_map);
}

std::tuple<Circuit, unit_map_t, unit_map_t> test_routing_method_circuit_no_perm(
    const Circuit& c, const ArchitecturePtr& a) {
  Circuit copy(c.n_qubits());
  copy.add_op<unsigned>(OpType::SWAP, {0, 1});
  copy.add_op<unsigned>(OpType::CX, {1, 0});
  copy.add_op<unsigned>(OpType::CX, {1, 0});

  std::vector<Qubit> qs = copy.all_qubits();
  std::vector<Node> ns = a->get_all_uids_vec();
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
  return std::make_tuple(copy, rename_map, final_map);
}

SCENARIO("Test RoutingMethodCircuit::check_method") {
  RoutingMethodCircuit rmc(
      test_routing_method_mf_swap_no_perm, test_check_method, 5, 5);
  Circuit c(2), circ3(3);
  c.add_op<unsigned>(OpType::CX, {0, 1});
  circ3.add_op<unsigned>(OpType::CX, {0, 2});
  circ3.add_op<unsigned>(OpType::CX, {2, 1});
  std::shared_ptr<MappingFrontier> mf2 = std::make_shared<MappingFrontier>(c);
  std::shared_ptr<MappingFrontier> mf3 =
      std::make_shared<MappingFrontier>(circ3);

  Architecture arc(
      {{Node("t", 1), Node("t", 0)}, {Node("t", 2), Node("t", 1)}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);

  REQUIRE(!rmc.check_method(mf2, shared_arc));
  REQUIRE(rmc.check_method(mf3, shared_arc));
}
SCENARIO("Test RoutingMethodCircuit::route_method") {
  Circuit comp(2);
  comp.add_op<unsigned>(OpType::SWAP, {0, 1});
  comp.add_op<unsigned>(OpType::CX, {1, 0});
  comp.add_op<unsigned>(OpType::CX, {1, 0});
  comp.add_op<unsigned>(OpType::CX, {1, 0});
  comp.add_op<unsigned>(OpType::CX, {1, 0});
  auto qbs = comp.all_qubits();
  unit_map_t rename_map = {{qbs[0], Node("t", 0)}, {qbs[1], Node("t", 1)}};
  comp.rename_units(rename_map);
  qubit_map_t permutation = {
      {Node("t", 0), Node("t", 1)}, {Node("t", 1), Node("t", 0)}};
  comp.permute_boundary_output(permutation);

  GIVEN("Non-implicit Permutation method, using MappingFrontier::add_swap") {
    RoutingMethodCircuit rmc(
        test_routing_method_mf_swap_no_perm, test_check_method, 2, 2);
    Circuit c(2);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});

    std::shared_ptr<MappingFrontier> mf = std::make_shared<MappingFrontier>(c);
    Architecture arc(
        {{Node("t", 1), Node("t", 0)}, {Node("t", 2), Node("t", 1)}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
    unit_map_t output = rmc.routing_method(mf, shared_arc), empty;
    REQUIRE(output == empty);
    REQUIRE(c == comp);
  }
  GIVEN("Non-implicit Permutation method, using circuit replacement") {
    RoutingMethodCircuit rmc(
        test_routing_method_circuit_no_perm, test_check_method, 2, 2);
    Circuit c(2);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});

    std::shared_ptr<MappingFrontier> mf = std::make_shared<MappingFrontier>(c);
    Architecture arc(
        {{Node("t", 1), Node("t", 0)}, {Node("t", 2), Node("t", 1)}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
    unit_map_t output = rmc.routing_method(mf, shared_arc), empty;
    REQUIRE(output == empty);
    REQUIRE(c == comp);
  }
  GIVEN("Implicit Permutation method, using MappingFrontier::add_swap") {
    RoutingMethodCircuit rmc(
        test_routing_method_mf_swap_perm, test_check_method, 2, 2);
    Circuit c(2);
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::CX, {0, 1});

    std::shared_ptr<MappingFrontier> mf = std::make_shared<MappingFrontier>(c);
    Architecture arc(
        {{Node("t", 1), Node("t", 0)}, {Node("t", 2), Node("t", 1)}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
    unit_map_t output = rmc.routing_method(mf, shared_arc), empty;
    REQUIRE(output == empty);

    Circuit comp1(2);
    comp1.add_op<unsigned>(OpType::SWAP, {0, 1});
    comp1.add_op<unsigned>(OpType::CX, {1, 0});
    comp1.add_op<unsigned>(OpType::CX, {1, 0});
    comp1.add_op<unsigned>(OpType::CX, {0, 1});
    comp1.add_op<unsigned>(OpType::CX, {0, 1});
    qbs = comp1.all_qubits();
    rename_map = {{qbs[0], Node("t", 0)}, {qbs[1], Node("t", 1)}};
    comp1.rename_units(rename_map);

    REQUIRE(c == comp1);
  }
}
}  // namespace tket