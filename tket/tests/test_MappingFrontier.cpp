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

#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "Circuit/ClassicalExpBox.hpp"
#include "Mapping/MappingManager.hpp"

namespace tket {

SCENARIO("Test MappingFrontier initialisation, advance_frontier_boundary.") {
  GIVEN("A typical Circuit and Architecture with uninitialised boundary") {
    Circuit circ;
    circ.add_q_register("test_nodes", 4);

    std::vector<Qubit> qubits = circ.all_qubits();

    Vertex v1 = circ.add_op<UnitID>(OpType::X, {qubits[0]});
    Vertex v8 = circ.add_op<UnitID>(OpType::S, {qubits[3]});
    Vertex v9 = circ.add_op<UnitID>(OpType::T, {qubits[3]});
    Vertex v2 = circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
    Vertex v3 = circ.add_op<UnitID>(OpType::CY, {qubits[2], qubits[3]});
    Vertex v4 = circ.add_op<UnitID>(OpType::H, {qubits[0]});
    Vertex v5 = circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[2]});
    Vertex v6 = circ.add_op<UnitID>(OpType::Y, {qubits[0]});
    Vertex v7 = circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[1]});

    std::vector<Node> nodes = {Node(0), Node(1), Node(2), Node(3)};

    Architecture arc(
        {{nodes[0], nodes[1]}, {nodes[1], nodes[3]}, {nodes[2], nodes[1]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);

    MappingFrontier m(circ);
    MappingFrontier mf(m);
    mf.advance_frontier_boundary(shared_arc);

    VertPort vp0 = mf.linear_boundary->get<TagKey>().find(nodes[0])->second;
    VertPort vp1 = mf.linear_boundary->get<TagKey>().find(nodes[1])->second;
    VertPort vp2 = mf.linear_boundary->get<TagKey>().find(nodes[2])->second;
    VertPort vp3 = mf.linear_boundary->get<TagKey>().find(nodes[3])->second;

    Edge e0 = mf.circuit_.get_nth_out_edge(vp0.first, vp0.second);
    Edge e1 = mf.circuit_.get_nth_out_edge(vp1.first, vp1.second);
    Edge e2 = mf.circuit_.get_nth_out_edge(vp2.first, vp2.second);
    Edge e3 = mf.circuit_.get_nth_out_edge(vp3.first, vp3.second);

    REQUIRE(mf.circuit_.source(e0) == v4);
    REQUIRE(mf.circuit_.target(e0) == v5);
    REQUIRE(mf.circuit_.source(e1) == v2);
    REQUIRE(mf.circuit_.target(e1) == v7);
    REQUIRE(
        mf.circuit_.get_OpType_from_Vertex(mf.circuit_.source(e2)) ==
        OpType::Input);
    REQUIRE(mf.circuit_.target(e2) == v3);
    REQUIRE(mf.circuit_.source(e3) == v9);
    REQUIRE(mf.circuit_.target(e3) == v3);

    mf.advance_frontier_boundary(shared_arc);
    // check that advance_frontier_boundary doesn't incorrectly move boundary
    // forwards
    vp0 = mf.linear_boundary->get<TagKey>().find(nodes[0])->second;
    vp1 = mf.linear_boundary->get<TagKey>().find(nodes[1])->second;
    vp2 = mf.linear_boundary->get<TagKey>().find(nodes[2])->second;
    vp3 = mf.linear_boundary->get<TagKey>().find(nodes[3])->second;

    e0 = mf.circuit_.get_nth_out_edge(vp0.first, vp0.second);
    e1 = mf.circuit_.get_nth_out_edge(vp1.first, vp1.second);
    e2 = mf.circuit_.get_nth_out_edge(vp2.first, vp2.second);
    e3 = mf.circuit_.get_nth_out_edge(vp3.first, vp3.second);

    REQUIRE(mf.circuit_.source(e0) == v4);
    REQUIRE(mf.circuit_.target(e0) == v5);
    REQUIRE(mf.circuit_.source(e1) == v2);
    REQUIRE(mf.circuit_.target(e1) == v7);
    REQUIRE(
        mf.circuit_.get_OpType_from_Vertex(mf.circuit_.source(e2)) ==
        OpType::Input);
    REQUIRE(mf.circuit_.target(e2) == v3);
    REQUIRE(mf.circuit_.source(e3) == v9);
    REQUIRE(mf.circuit_.target(e3) == v3);
  }

  GIVEN("A circuit with measurements and classically controlled operations") {
    Circuit circ(3, 1);
    std::vector<Qubit> qubits = circ.all_qubits();
    // All gates are physically permitted
    Vertex v0 = circ.add_op<unsigned>(OpType::Measure, {0, 0});
    Vertex v1 =
        circ.add_conditional_gate<unsigned>(OpType::Rx, {0.6}, {0}, {0}, 1);
    Vertex v2 =
        circ.add_conditional_gate<unsigned>(OpType::Rz, {0.6}, {1}, {0}, 1);
    Vertex v3 = circ.add_op<unsigned>(OpType::X, {2});
    std::vector<Node> nodes = {Node(0), Node(1), Node(2)};

    Architecture arc({{nodes[0], nodes[1]}, {nodes[1], nodes[2]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]}};
    circ.rename_units(rename_map);
    MappingFrontier mf(circ);
    mf.advance_frontier_boundary(shared_arc);
    VertPort vp0 = mf.linear_boundary->get<TagKey>().find(nodes[0])->second;
    VertPort vp1 = mf.linear_boundary->get<TagKey>().find(nodes[1])->second;
    VertPort vp2 = mf.linear_boundary->get<TagKey>().find(nodes[2])->second;
    Op_ptr op = circ.get_Op_ptr_from_Vertex(vp0.first);
    Op_ptr op2 = circ.get_Op_ptr_from_Vertex(vp1.first);
    Op_ptr op3 = circ.get_Op_ptr_from_Vertex(vp2.first);
    REQUIRE(vp0.first == v1);
    REQUIRE(vp1.first == v2);
    REQUIRE(vp2.first == v3);
  }
  GIVEN(
      "A circuit with multi edge bundles of booleans, conditional gates with "
      "multiple inputs, conditional 2-qubit gates.") {
    Circuit circ(4, 4);

    Vertex v0 =
        circ.add_conditional_gate<unsigned>(OpType::X, {}, {0}, {0, 1}, 1);
    Vertex v1 = circ.add_conditional_gate<unsigned>(OpType::Y, {}, {1}, {1}, 0);
    Vertex v2 = circ.add_op<unsigned>(OpType::CX, {1, 2});
    Vertex v3 = circ.add_measure(2, 2);
    Vertex v4 = circ.add_op<unsigned>(OpType::CX, {3, 2});
    Vertex v5 = circ.add_measure(3, 3);
    Vertex v6 =
        circ.add_conditional_gate<unsigned>(OpType::Z, {}, {3}, {1, 2}, 0);
    Vertex v7 = circ.add_measure(3, 3);
    Vertex v8 = circ.add_barrier(
        {Qubit(0), Qubit(1), Qubit(2), Qubit(3), Bit(1), Bit(2), Bit(3)});
    Vertex v9 =
        circ.add_conditional_gate<unsigned>(OpType::Z, {}, {3}, {1, 2}, 0);

    std::vector<Node> nodes = {Node(0), Node(1), Node(2), Node(3)};
    Architecture arc(
        {{nodes[0], nodes[1]}, {nodes[1], nodes[2]}, {nodes[2], nodes[3]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
    std::vector<Qubit> qubits = circ.all_qubits();
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};

    circ.rename_units(rename_map);
    std::vector<Bit> bits = circ.all_bits();
    MappingFrontier mf(circ);

    REQUIRE(
        mf.boolean_boundary->get<TagKey>().find(bits[0]) !=
        mf.boolean_boundary->get<TagKey>().end());
    REQUIRE(
        mf.boolean_boundary->get<TagKey>().find(bits[1]) !=
        mf.boolean_boundary->get<TagKey>().end());
    REQUIRE(
        mf.boolean_boundary->get<TagKey>().find(bits[2]) ==
        mf.boolean_boundary->get<TagKey>().end());
    REQUIRE(
        mf.boolean_boundary->get<TagKey>().find(bits[3]) ==
        mf.boolean_boundary->get<TagKey>().end());

    mf.advance_frontier_boundary(shared_arc);

    VertPort vp_q_0 = mf.linear_boundary->get<TagKey>().find(nodes[0])->second;
    VertPort vp_q_1 = mf.linear_boundary->get<TagKey>().find(nodes[1])->second;
    VertPort vp_q_2 = mf.linear_boundary->get<TagKey>().find(nodes[2])->second;
    VertPort vp_q_3 = mf.linear_boundary->get<TagKey>().find(nodes[3])->second;
    // note c[0] and c[1] not linear_boundary as they are immediately boolean
    VertPort vp_b_2 = mf.linear_boundary->get<TagKey>().find(bits[2])->second;
    VertPort vp_b_3 = mf.linear_boundary->get<TagKey>().find(bits[3])->second;

    REQUIRE(
        circ.get_OpType_from_Vertex(circ.target(circ.get_nth_out_edge(
            vp_q_0.first, vp_q_0.second))) == OpType::Output);
    REQUIRE(
        circ.get_OpType_from_Vertex(circ.target(circ.get_nth_out_edge(
            vp_q_1.first, vp_q_1.second))) == OpType::Output);
    REQUIRE(
        circ.get_OpType_from_Vertex(circ.target(circ.get_nth_out_edge(
            vp_q_2.first, vp_q_2.second))) == OpType::Output);
    REQUIRE(
        circ.get_OpType_from_Vertex(circ.target(circ.get_nth_out_edge(
            vp_q_3.first, vp_q_3.second))) == OpType::Output);
    REQUIRE(
        circ.get_OpType_from_Vertex(circ.target(circ.get_nth_out_edge(
            vp_b_2.first, vp_b_2.second))) == OpType::ClOutput);
    REQUIRE(
        circ.get_OpType_from_Vertex(circ.target(circ.get_nth_out_edge(
            vp_b_3.first, vp_b_3.second))) == OpType::ClOutput);

    // in and then removed from boolean boundary
    REQUIRE(
        mf.boolean_boundary->get<TagKey>().find(bits[2]) ==
        mf.boolean_boundary->get<TagKey>().end());
    // not in boolean boundary because bool not used in condition
    REQUIRE(
        mf.boolean_boundary->get<TagKey>().find(bits[3]) ==
        mf.boolean_boundary->get<TagKey>().end());
  }
}

SCENARIO("Test MappingFrontier get_default_to_linear_boundary_unit_map") {
  Circuit circ;
  circ.add_q_register("test_nodes", 4);
  std::vector<Qubit> qubits = circ.all_qubits();
  MappingFrontier mf(circ);
  unit_map_t d_2_q = mf.get_default_to_linear_boundary_unit_map();
  REQUIRE(d_2_q[Qubit(0)] == qubits[0]);
  REQUIRE(d_2_q[Qubit(1)] == qubits[1]);
  REQUIRE(d_2_q[Qubit(2)] == qubits[2]);
  REQUIRE(d_2_q[Qubit(3)] == qubits[3]);
}

SCENARIO("Test MappingFrontier get_frontier_subcircuit.") {
  GIVEN(
      "A typical circuit, MappingFrontier with depth 1 and depth 3 "
      "subcircuit returns, no renaming units.") {
    Circuit circ;
    circ.add_q_register("test_nodes", 4);
    std::vector<Qubit> qubits = circ.all_qubits();

    Vertex v1 = circ.add_op<UnitID>(OpType::X, {qubits[0]});
    Vertex v8 = circ.add_op<UnitID>(OpType::S, {qubits[3]});
    Vertex v9 = circ.add_op<UnitID>(OpType::T, {qubits[3]});
    Vertex v2 = circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
    Vertex v3 = circ.add_op<UnitID>(OpType::CY, {qubits[2], qubits[3]});
    Vertex v4 = circ.add_op<UnitID>(OpType::H, {qubits[0]});
    Vertex v5 = circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[2]});
    Vertex v6 = circ.add_op<UnitID>(OpType::Y, {qubits[0]});
    Vertex v7 = circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[1]});

    std::vector<Node> nodes = {Node(0), Node(1), Node(2), Node(3)};

    Architecture arc(
        {{nodes[0], nodes[1]}, {nodes[1], nodes[3]}, {nodes[2], nodes[1]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);

    MappingFrontier mf_1(circ);
    MappingFrontier mf_3(circ);

    mf_1.advance_frontier_boundary(shared_arc);
    Subcircuit sc1 = mf_1.get_frontier_subcircuit(1, 7);
    mf_3.advance_frontier_boundary(shared_arc);
    Subcircuit sc3 = mf_3.get_frontier_subcircuit(3, 7);

    Circuit frontier_circuit_1 = mf_1.circuit_.subcircuit(sc1);

    Circuit comparison_circuit(4);
    comparison_circuit.add_op<unsigned>(OpType::CY, {2, 3});
    REQUIRE(frontier_circuit_1 == comparison_circuit);

    Circuit frontier_circuit_3 = mf_3.circuit_.subcircuit(sc3);
    comparison_circuit.add_op<unsigned>(OpType::CZ, {0, 2});
    comparison_circuit.add_op<unsigned>(OpType::Y, {0});
    comparison_circuit.add_op<unsigned>(OpType::CX, {3, 1});
    REQUIRE(frontier_circuit_3 == comparison_circuit);
  }

  GIVEN(
      "A typical circuit but with non-contiguous Qubit Labelling. "
      "MappingFrontier with depth 1 and depth 3 "
      "subcircuit returns, no renaming units.") {
    Circuit circ(4);
    Qubit q0("label_0", 1);
    Qubit q1("label_1", 3);
    Qubit q2("label_2", 0);
    Qubit q3("label_3", 2);
    std::vector<Qubit> qubits = {q0, q1, q2, q3};
    std::map<UnitID, UnitID> new_units = {
        {Qubit(0), q0}, {Qubit(1), q1}, {Qubit(2), q2}, {Qubit(3), q3}};
    circ.rename_units(new_units);

    Vertex v1 = circ.add_op<UnitID>(OpType::X, {qubits[0]});
    Vertex v8 = circ.add_op<UnitID>(OpType::S, {qubits[3]});
    Vertex v9 = circ.add_op<UnitID>(OpType::T, {qubits[3]});
    Vertex v2 = circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
    Vertex v3 = circ.add_op<UnitID>(OpType::CY, {qubits[2], qubits[3]});
    Vertex v4 = circ.add_op<UnitID>(OpType::H, {qubits[0]});
    Vertex v5 = circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[2]});
    Vertex v6 = circ.add_op<UnitID>(OpType::Y, {qubits[0]});
    Vertex v7 = circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[1]});

    std::vector<Node> nodes = {Node(0), Node(1), Node(2), Node(3)};

    Architecture arc(
        {{nodes[0], nodes[1]}, {nodes[1], nodes[3]}, {nodes[2], nodes[1]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);

    std::map<UnitID, UnitID> rename_map = {
        {q0, nodes[0]}, {q1, nodes[1]}, {q2, nodes[2]}, {q3, nodes[3]}};

    circ.rename_units(rename_map);

    MappingFrontier mf_1(circ);
    MappingFrontier mf_3(circ);

    mf_1.advance_frontier_boundary(shared_arc);
    Subcircuit sc1 = mf_1.get_frontier_subcircuit(1, 7);
    mf_3.advance_frontier_boundary(shared_arc);
    Subcircuit sc3 = mf_3.get_frontier_subcircuit(3, 7);

    Circuit frontier_circuit_1 = mf_1.circuit_.subcircuit(sc1);

    frontier_circuit_1.rename_units(
        mf_1.get_default_to_linear_boundary_unit_map());
    Circuit comparison_circuit(4);
    std::map<UnitID, UnitID> rename_map_default = {
        {Qubit(0), nodes[0]},
        {Qubit(1), nodes[1]},
        {Qubit(2), nodes[2]},
        {Qubit(3), nodes[3]}};
    comparison_circuit.rename_units(rename_map_default);
    comparison_circuit.add_op<UnitID>(OpType::CY, {nodes[2], nodes[3]});
    REQUIRE(frontier_circuit_1 == comparison_circuit);
    Circuit frontier_circuit_3 = mf_3.circuit_.subcircuit(sc3);
    frontier_circuit_3.rename_units(
        mf_3.get_default_to_linear_boundary_unit_map());

    comparison_circuit.add_op<UnitID>(OpType::CZ, {nodes[0], nodes[2]});
    comparison_circuit.add_op<UnitID>(OpType::Y, {nodes[0]});
    comparison_circuit.add_op<UnitID>(OpType::CX, {nodes[3], nodes[1]});
    REQUIRE(frontier_circuit_3 == comparison_circuit);
  }
}

SCENARIO("Test update_linear_boundary_uids.") {
  Circuit circ(10);
  std::vector<Qubit> qbs = circ.all_qubits();
  MappingFrontier mf(circ);
  GIVEN("Empty relabelling.") { mf.update_linear_boundary_uids({}); }
  GIVEN("Relabel some qubits to same qubit.") {
    mf.update_linear_boundary_uids(
        {{qbs[0], qbs[0]}, {qbs[2], qbs[2]}, {qbs[7], qbs[7]}});
    REQUIRE(mf.linear_boundary->get<TagKey>().find(qbs[0])->first == qbs[0]);
    REQUIRE(mf.linear_boundary->get<TagKey>().find(qbs[2])->first == qbs[2]);
    REQUIRE(mf.linear_boundary->get<TagKey>().find(qbs[7])->first == qbs[7]);
  }
  GIVEN("Relabel to already present qubit, check boundary has qubit removed.") {
    mf.update_linear_boundary_uids({{qbs[0], qbs[1]}});
    REQUIRE(mf.linear_boundary->get<TagKey>().size() == 9);
  }
  GIVEN("Relabel to new UnitID.") {
    mf.update_linear_boundary_uids({{qbs[0], Node("tn", 6)}});
    REQUIRE(
        mf.linear_boundary->get<TagKey>().find(qbs[0]) ==
        mf.linear_boundary->get<TagKey>().end());
  }
}

SCENARIO("Test permute_subcircuit_q_out_hole.") {
  GIVEN("Quantum Boundary and Permutation have size mismatch.") {
    Circuit circ(0);
    circ.add_q_register("test_nodes", 4);
    Qubit q0("test_nodes", 0);
    Qubit q1("test_nodes", 1);
    Qubit q2("test_nodes", 2);
    Qubit q3("test_nodes", 3);

    circ.add_op<UnitID>(OpType::X, {q0});
    circ.add_op<UnitID>(OpType::CX, {q0, q1});
    circ.add_op<UnitID>(OpType::CY, {q2, q3});
    circ.add_op<UnitID>(OpType::CZ, {q0, q2});
    circ.add_op<UnitID>(OpType::CX, {q3, q1});

    std::vector<Node> nodes = {Node(0), Node(1), Node(2), Node(3)};

    Architecture arc(
        {{nodes[0], nodes[1]}, {nodes[1], nodes[3]}, {nodes[2], nodes[1]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);

    std::map<UnitID, UnitID> rename_map = {
        {q0, nodes[0]}, {q1, nodes[1]}, {q2, nodes[2]}, {q3, nodes[3]}};
    circ.rename_units(rename_map);

    MappingFrontier mf(circ);

    mf.advance_frontier_boundary(shared_arc);
    Subcircuit sc = mf.get_frontier_subcircuit(2, 5);
    unit_map_t permutation = {{nodes[0], nodes[1]}};

    REQUIRE_THROWS_AS(
        mf.permute_subcircuit_q_out_hole(permutation, sc),
        MappingFrontierError);
  }
  GIVEN(
      "Quantum Boundary and permutation have same size, but UnitID don't "
      "match.") {
    Circuit circ(0);
    circ.add_q_register("test_nodes", 4);
    Qubit q0("test_nodes", 0);
    Qubit q1("test_nodes", 1);
    Qubit q2("test_nodes", 2);
    Qubit q3("test_nodes", 3);

    circ.add_op<UnitID>(OpType::X, {q0});
    circ.add_op<UnitID>(OpType::CX, {q0, q1});
    circ.add_op<UnitID>(OpType::CY, {q2, q3});
    circ.add_op<UnitID>(OpType::CZ, {q0, q2});
    circ.add_op<UnitID>(OpType::CX, {q3, q1});

    std::vector<Node> nodes = {Node(0), Node(1), Node(2), Node(3)};

    Architecture arc(
        {{nodes[0], nodes[1]}, {nodes[1], nodes[3]}, {nodes[2], nodes[1]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);

    std::map<UnitID, UnitID> rename_map = {
        {q0, nodes[0]}, {q1, nodes[1]}, {q2, nodes[2]}, {q3, nodes[3]}};
    circ.rename_units(rename_map);

    MappingFrontier mf(circ);

    mf.advance_frontier_boundary(shared_arc);
    Subcircuit sc = mf.get_frontier_subcircuit(2, 5);
    unit_map_t permutation = {
        {nodes[0], nodes[1]},
        {nodes[1], nodes[2]},
        {nodes[2], nodes[3]},
        {Node(4), nodes[0]}};

    REQUIRE_THROWS_AS(
        mf.permute_subcircuit_q_out_hole(permutation, sc),
        MappingFrontierError);
  }
  GIVEN("A four qubit subcircuit where every qubit is permuted by given map.") {
    Circuit circ(0);
    circ.add_q_register("test_nodes", 4);
    Qubit q0("test_nodes", 0);
    Qubit q1("test_nodes", 1);
    Qubit q2("test_nodes", 2);
    Qubit q3("test_nodes", 3);

    circ.add_op<UnitID>(OpType::X, {q0});
    circ.add_op<UnitID>(OpType::CX, {q0, q1});
    circ.add_op<UnitID>(OpType::CY, {q2, q3});
    circ.add_op<UnitID>(OpType::CZ, {q0, q2});
    circ.add_op<UnitID>(OpType::CX, {q3, q1});

    std::vector<Node> nodes = {Node(0), Node(1), Node(2), Node(3)};

    Architecture arc(
        {{nodes[0], nodes[1]}, {nodes[1], nodes[3]}, {nodes[2], nodes[1]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);

    std::map<UnitID, UnitID> rename_map = {
        {q0, nodes[0]}, {q1, nodes[1]}, {q2, nodes[2]}, {q3, nodes[3]}};
    circ.rename_units(rename_map);

    MappingFrontier mf(circ);

    mf.advance_frontier_boundary(shared_arc);
    Subcircuit sc = mf.get_frontier_subcircuit(2, 5);
    // assume only 1 subcircuit
    EdgeVec original_q_out = sc.q_out_hole;

    unit_map_t permutation = {
        {nodes[0], nodes[1]},
        {nodes[1], nodes[2]},
        {nodes[2], nodes[3]},
        {nodes[3], nodes[0]}};
    mf.permute_subcircuit_q_out_hole(permutation, sc);

    EdgeVec permuted_q_out = sc.q_out_hole;

    REQUIRE(original_q_out[1] == permuted_q_out[0]);
    REQUIRE(original_q_out[2] == permuted_q_out[1]);
    REQUIRE(original_q_out[3] == permuted_q_out[2]);
    REQUIRE(original_q_out[0] == permuted_q_out[3]);
  }
  GIVEN("A four qubit subcircuit with a partial permutation.") {
    Circuit circ(0);
    circ.add_q_register("test_nodes", 4);
    Qubit q0("test_nodes", 0);
    Qubit q1("test_nodes", 1);
    Qubit q2("test_nodes", 2);
    Qubit q3("test_nodes", 3);

    Vertex v1 = circ.add_op<UnitID>(OpType::X, {q0});
    Vertex v2 = circ.add_op<UnitID>(OpType::CX, {q0, q1});
    Vertex v3 = circ.add_op<UnitID>(OpType::CY, {q2, q3});
    Vertex v5 = circ.add_op<UnitID>(OpType::CZ, {q0, q2});
    Vertex v7 = circ.add_op<UnitID>(OpType::CX, {q3, q1});

    std::vector<Node> nodes = {Node(0), Node(1), Node(2), Node(3)};

    Architecture arc(
        {{nodes[0], nodes[1]}, {nodes[1], nodes[3]}, {nodes[2], nodes[1]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);

    std::map<UnitID, UnitID> rename_map = {
        {q0, nodes[0]}, {q1, nodes[1]}, {q2, nodes[2]}, {q3, nodes[3]}};
    circ.rename_units(rename_map);

    MappingFrontier mf(circ);

    mf.advance_frontier_boundary(shared_arc);
    Subcircuit sc = mf.get_frontier_subcircuit(2, 5);
    // assume only 1 subcircuit
    EdgeVec original_q_out = sc.q_out_hole;

    unit_map_t permutation = {
        {nodes[0], nodes[1]},
        {nodes[1], nodes[0]},
        {nodes[2], nodes[2]},
        {nodes[3], nodes[3]}};
    mf.permute_subcircuit_q_out_hole(permutation, sc);

    EdgeVec permuted_q_out = sc.q_out_hole;

    REQUIRE(original_q_out[1] == permuted_q_out[0]);
    REQUIRE(original_q_out[0] == permuted_q_out[1]);
    REQUIRE(original_q_out[2] == permuted_q_out[2]);
    REQUIRE(original_q_out[3] == permuted_q_out[3]);
  }
}
SCENARIO("Test MappingFrontier::advance_next_2qb_slice") {
  std::vector<Node> nodes = {Node("test_node", 0), Node("test_node", 1),
                             Node("test_node", 2), Node("node_test", 3),
                             Node("node_test", 4), Node("node_test", 5),
                             Node("test_node", 6), Node("node_test", 7)};
  // n0 -- n1 -- n2 -- n3 -- n4
  //             |     |
  //             n5    n7
  //             |
  //             n6
  Architecture architecture(
      {{nodes[0], nodes[1]},
       {nodes[1], nodes[2]},
       {nodes[2], nodes[3]},
       {nodes[3], nodes[4]},
       {nodes[2], nodes[5]},
       {nodes[5], nodes[6]},
       {nodes[3], nodes[7]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  GIVEN("One CX to find in next slice.") {
    Circuit circ(8);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::X, {qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[7]});
    //                   n7
    //                   |
    // n0 -- n1 -- n2 -- n3 -- n4
    //             |
    //             n5
    //             |
    //             n6
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]}, {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]}};
    circ.rename_units(rename_map);
    MappingFrontier mf(circ);
    // gets to first two cx
    mf.advance_frontier_boundary(shared_arc);

    VertPort vp0 = mf.linear_boundary->get<TagKey>().find(nodes[0])->second;
    VertPort vp4 = mf.linear_boundary->get<TagKey>().find(nodes[4])->second;
    VertPort vp6 = mf.linear_boundary->get<TagKey>().find(nodes[6])->second;
    VertPort vp7 = mf.linear_boundary->get<TagKey>().find(nodes[7])->second;

    Edge e0 = mf.circuit_.get_nth_out_edge(vp0.first, vp0.second);
    Edge e4 = mf.circuit_.get_nth_out_edge(vp4.first, vp4.second);
    Edge e6 = mf.circuit_.get_nth_out_edge(vp6.first, vp6.second);
    Edge e7 = mf.circuit_.get_nth_out_edge(vp7.first, vp7.second);

    Vertex v0 = mf.circuit_.target(e0);
    Vertex v4 = mf.circuit_.target(e4);
    Vertex v6 = mf.circuit_.target(e6);
    Vertex v7 = mf.circuit_.target(e7);

    REQUIRE(v0 == v4);
    REQUIRE(v6 == v7);

    mf.advance_next_2qb_slice(5);
    VertPort vp2 = mf.linear_boundary->get<TagKey>().find(nodes[2])->second;
    vp7 = mf.linear_boundary->get<TagKey>().find(nodes[7])->second;

    Edge e2 = mf.circuit_.get_nth_out_edge(vp2.first, vp2.second);
    e7 = mf.circuit_.get_nth_out_edge(vp7.first, vp7.second);

    Vertex v2 = mf.circuit_.target(e2);
    v7 = mf.circuit_.target(e7);

    REQUIRE(v2 == v7);
  }
  GIVEN(
      "Three CX to find in next slice 1, Two CX and one CZ in next slice 2. ") {
    Circuit circ(8);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[5]});
    circ.add_op<UnitID>(OpType::X, {qubits[0]});
    circ.add_op<UnitID>(OpType::CX, {qubits[4], qubits[1]});
    circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[0]});
    circ.add_op<UnitID>(OpType::X, {qubits[1]});
    circ.add_op<UnitID>(OpType::CX, {qubits[4], qubits[1]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[3], qubits[7]});
    //                   n7
    //                   |
    // n0 -- n1 -- n2 -- n3 -- n4
    //             |
    //             n5
    //             |
    //             n6
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]}, {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]}};
    circ.rename_units(rename_map);
    MappingFrontier mf(circ);
    // gets to first two cx
    mf.advance_frontier_boundary(shared_arc);

    VertPort vp0 = mf.linear_boundary->get<TagKey>().find(nodes[0])->second;
    VertPort vp4 = mf.linear_boundary->get<TagKey>().find(nodes[4])->second;
    VertPort vp6 = mf.linear_boundary->get<TagKey>().find(nodes[6])->second;
    VertPort vp7 = mf.linear_boundary->get<TagKey>().find(nodes[7])->second;

    Edge e0 = mf.circuit_.get_nth_out_edge(vp0.first, vp0.second);
    Edge e4 = mf.circuit_.get_nth_out_edge(vp4.first, vp4.second);
    Edge e6 = mf.circuit_.get_nth_out_edge(vp6.first, vp6.second);
    Edge e7 = mf.circuit_.get_nth_out_edge(vp7.first, vp7.second);

    Vertex v0 = mf.circuit_.target(e0);
    Vertex v4 = mf.circuit_.target(e4);
    Vertex v6 = mf.circuit_.target(e6);
    Vertex v7 = mf.circuit_.target(e7);

    // get edges
    // then get target...
    REQUIRE(v0 == v4);
    REQUIRE(v6 == v7);

    mf.advance_next_2qb_slice(1);
    vp0 = mf.linear_boundary->get<TagKey>().find(nodes[0])->second;
    VertPort vp1 = mf.linear_boundary->get<TagKey>().find(nodes[1])->second;
    VertPort vp2 = mf.linear_boundary->get<TagKey>().find(nodes[2])->second;
    vp4 = mf.linear_boundary->get<TagKey>().find(nodes[4])->second;
    VertPort vp5 = mf.linear_boundary->get<TagKey>().find(nodes[5])->second;
    vp7 = mf.linear_boundary->get<TagKey>().find(nodes[7])->second;

    e0 = mf.circuit_.get_nth_out_edge(vp0.first, vp0.second);
    Edge e1 = mf.circuit_.get_nth_out_edge(vp1.first, vp1.second);
    Edge e2 = mf.circuit_.get_nth_out_edge(vp2.first, vp2.second);
    e4 = mf.circuit_.get_nth_out_edge(vp4.first, vp4.second);
    Edge e5 = mf.circuit_.get_nth_out_edge(vp5.first, vp5.second);
    e7 = mf.circuit_.get_nth_out_edge(vp7.first, vp7.second);

    v0 = mf.circuit_.target(e0);
    Vertex v1 = mf.circuit_.target(e1);
    Vertex v2 = mf.circuit_.target(e2);
    v4 = mf.circuit_.target(e4);
    Vertex v5 = mf.circuit_.target(e5);
    v7 = mf.circuit_.target(e7);

    REQUIRE(v1 == v4);
    REQUIRE(v0 == v5);
    REQUIRE(v2 == v7);

    mf.advance_next_2qb_slice(1);
    vp0 = mf.linear_boundary->get<TagKey>().find(nodes[0])->second;
    vp1 = mf.linear_boundary->get<TagKey>().find(nodes[1])->second;
    vp2 = mf.linear_boundary->get<TagKey>().find(nodes[2])->second;
    VertPort vp3 = mf.linear_boundary->get<TagKey>().find(nodes[3])->second;
    vp4 = mf.linear_boundary->get<TagKey>().find(nodes[4])->second;
    vp7 = mf.linear_boundary->get<TagKey>().find(nodes[7])->second;

    e0 = mf.circuit_.get_nth_out_edge(vp0.first, vp0.second);
    e1 = mf.circuit_.get_nth_out_edge(vp1.first, vp1.second);
    e2 = mf.circuit_.get_nth_out_edge(vp2.first, vp2.second);
    Edge e3 = mf.circuit_.get_nth_out_edge(vp3.first, vp3.second);
    e4 = mf.circuit_.get_nth_out_edge(vp4.first, vp4.second);
    e7 = mf.circuit_.get_nth_out_edge(vp7.first, vp7.second);

    v0 = mf.circuit_.target(e0);
    v1 = mf.circuit_.target(e1);
    v2 = mf.circuit_.target(e2);
    Vertex v3 = mf.circuit_.target(e3);
    v4 = mf.circuit_.target(e4);
    v7 = mf.circuit_.target(e7);

    REQUIRE(v0 == v2);
    REQUIRE(v1 == v4);
    REQUIRE(v3 == v7);
  }
}
SCENARIO("Test MappingFrontier::add_qubit") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3)};
  Circuit circ(3);
  std::vector<Qubit> qubits = circ.all_qubits();
  circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
  circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[2]});
  std::map<UnitID, UnitID> rename_map = {
      {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]}};
  circ.rename_units(rename_map);

  MappingFrontier mf(circ);
  mf.add_ancilla(nodes[3]);

  REQUIRE(circ.all_qubits().size() == 4);
  REQUIRE(mf.circuit_.all_qubits().size() == 4);
  REQUIRE(mf.linear_boundary->size() == 4);
  REQUIRE(mf.linear_boundary->find(nodes[3]) != mf.linear_boundary->end());
}

SCENARIO("Test MappingFrontier::add_swap") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3)};
  Circuit circ(4);
  std::vector<Qubit> qubits = circ.all_qubits();
  circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
  circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[2]});
  circ.add_op<UnitID>(OpType::CZ, {qubits[1], qubits[3]});

  std::map<UnitID, UnitID> rename_map = {
      {qubits[0], nodes[0]},
      {qubits[1], nodes[1]},
      {qubits[2], nodes[2]},
      {qubits[3], nodes[3]}};
  circ.rename_units(rename_map);
  MappingFrontier mf(circ);
  REQUIRE(mf.add_swap(nodes[0], nodes[1]));

  std::vector<Command> commands = mf.circuit_.get_commands();
  REQUIRE(commands.size() == 4);
  Command swap_c = commands[0];
  unit_vector_t uids = {nodes[0], nodes[1]};
  REQUIRE(swap_c.get_args() == uids);
  REQUIRE(*swap_c.get_op_ptr() == *get_op_ptr(OpType::SWAP));

  Command cx_c = commands[1];
  uids = {nodes[1], nodes[0]};
  REQUIRE(cx_c.get_args() == uids);
  REQUIRE(*cx_c.get_op_ptr() == *get_op_ptr(OpType::CX));

  cx_c = commands[2];
  uids = {nodes[0], nodes[2]};
  REQUIRE(cx_c.get_args() == uids);
  REQUIRE(*cx_c.get_op_ptr() == *get_op_ptr(OpType::CX));

  cx_c = commands[3];
  uids = {nodes[0], nodes[3]};
  REQUIRE(cx_c.get_args() == uids);
  REQUIRE(*cx_c.get_op_ptr() == *get_op_ptr(OpType::CZ));

  Node new_node("new_node", 8);
  REQUIRE(mf.add_swap(nodes[0], new_node));

  commands = mf.circuit_.get_commands();
  REQUIRE(commands.size() == 5);
  swap_c = commands[0];
  uids = {nodes[0], nodes[1]};

  REQUIRE(swap_c.get_args() == uids);
  REQUIRE(*swap_c.get_op_ptr() == *get_op_ptr(OpType::SWAP));

  swap_c = commands[1];
  uids = {nodes[0], new_node};
  REQUIRE(swap_c.get_args() == uids);
  REQUIRE(*swap_c.get_op_ptr() == *get_op_ptr(OpType::SWAP));

  cx_c = commands[2];
  uids = {nodes[1], new_node};
  REQUIRE(cx_c.get_args() == uids);
  REQUIRE(*cx_c.get_op_ptr() == *get_op_ptr(OpType::CX));

  cx_c = commands[3];
  uids = {new_node, nodes[2]};
  REQUIRE(cx_c.get_args() == uids);
  REQUIRE(*cx_c.get_op_ptr() == *get_op_ptr(OpType::CX));

  cx_c = commands[4];
  uids = {new_node, nodes[3]};
  REQUIRE(cx_c.get_args() == uids);
  REQUIRE(*cx_c.get_op_ptr() == *get_op_ptr(OpType::CZ));

  // swap on same pair of nodes returns false
  REQUIRE(!mf.add_swap(nodes[0], new_node));
}

SCENARIO("Test MappingFrontier::add_swap, classical wires edge case") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3)};
  Circuit circ(4, 3);
  std::vector<Qubit> qubits = circ.all_qubits();
  std::vector<Bit> bits = circ.all_bits();
  circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
  circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[3]});
  circ.add_measure(3, 0);
  circ.add_conditional_gate<UnitID>(
      OpType::Y, {}, {qubits[2]}, {bits[0], bits[1], bits[2]}, 3);
  circ.add_conditional_gate<UnitID>(OpType::X, {}, {qubits[1]}, {bits[2]}, 1);
  circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[0]});
  circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[0]});

  Architecture architecture(
      {{nodes[0], nodes[1]}, {nodes[0], nodes[2]}, {nodes[0], nodes[3]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  MappingFrontier mf(circ);
  mf.advance_frontier_boundary(shared_arc);
  REQUIRE(mf.add_swap(qubits[0], qubits[2]));
}
SCENARIO("Test MappingFrontier::add_bridge") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3)};
  Circuit circ(4);
  std::vector<Qubit> qubits = circ.all_qubits();
  circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
  circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[2]});
  circ.add_op<UnitID>(OpType::CZ, {qubits[1], qubits[3]});
  std::map<UnitID, UnitID> rename_map = {
      {qubits[0], nodes[0]},
      {qubits[1], nodes[1]},
      {qubits[2], nodes[2]},
      {qubits[3], nodes[3]}};
  circ.rename_units(rename_map);
  MappingFrontier mf(circ);
  mf.add_bridge(nodes[0], nodes[2], nodes[1]);

  std::vector<Command> commands = mf.circuit_.get_commands();
  REQUIRE(commands.size() == 3);
  Command bridge_c = commands[0];
  unit_vector_t uids = {nodes[0], nodes[2], nodes[1]};
  REQUIRE(bridge_c.get_args() == uids);
  REQUIRE(*bridge_c.get_op_ptr() == *get_op_ptr(OpType::BRIDGE));

  Command cx_c = commands[1];
  uids = {nodes[1], nodes[2]};
  REQUIRE(cx_c.get_args() == uids);
  REQUIRE(*cx_c.get_op_ptr() == *get_op_ptr(OpType::CX));

  cx_c = commands[2];
  uids = {nodes[1], nodes[3]};
  REQUIRE(cx_c.get_args() == uids);
  REQUIRE(*cx_c.get_op_ptr() == *get_op_ptr(OpType::CZ));
}
SCENARIO("Test MappingFrontier set_linear_boundary") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3)};
  Architecture architecture(
      {{nodes[0], nodes[1]}, {nodes[1], nodes[2]}, {nodes[2], nodes[3]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  Circuit circ(4);
  std::vector<Qubit> qubits = circ.all_qubits();
  circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
  circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[2]});
  circ.add_op<UnitID>(OpType::CZ, {qubits[2], qubits[3]});
  std::map<UnitID, UnitID> rename_map = {
      {qubits[0], nodes[0]},
      {qubits[1], nodes[1]},
      {qubits[2], nodes[2]},
      {qubits[3], nodes[3]}};
  circ.rename_units(rename_map);
  MappingFrontier mf(circ);

  unit_vertport_frontier_t copy;
  for (const std::pair<UnitID, VertPort>& pair :
       mf.linear_boundary->get<TagKey>()) {
    copy.insert({pair.first, pair.second});
  }

  VertPort vp0_c = copy.get<TagKey>().find(nodes[0])->second;
  VertPort vp1_c = copy.get<TagKey>().find(nodes[1])->second;
  VertPort vp2_c = copy.get<TagKey>().find(nodes[2])->second;
  VertPort vp3_c = copy.get<TagKey>().find(nodes[3])->second;

  mf.advance_frontier_boundary(shared_arc);

  VertPort vp0_in = mf.linear_boundary->get<TagKey>().find(nodes[0])->second;
  VertPort vp1_in = mf.linear_boundary->get<TagKey>().find(nodes[1])->second;
  VertPort vp2_in = mf.linear_boundary->get<TagKey>().find(nodes[2])->second;
  VertPort vp3_in = mf.linear_boundary->get<TagKey>().find(nodes[3])->second;

  REQUIRE(vp0_in.first != vp0_c.first);
  REQUIRE(vp1_in.first != vp1_c.first);
  REQUIRE(vp2_in.first != vp2_c.first);
  REQUIRE(vp3_in.first != vp3_c.first);

  mf.set_linear_boundary(copy);

  vp0_in = mf.linear_boundary->get<TagKey>().find(nodes[0])->second;
  vp1_in = mf.linear_boundary->get<TagKey>().find(nodes[1])->second;
  vp2_in = mf.linear_boundary->get<TagKey>().find(nodes[2])->second;
  vp3_in = mf.linear_boundary->get<TagKey>().find(nodes[3])->second;

  REQUIRE(vp0_in.first == vp0_c.first);
  REQUIRE(vp1_in.first == vp1_c.first);
  REQUIRE(vp2_in.first == vp2_c.first);
  REQUIRE(vp3_in.first == vp3_c.first);
}

SCENARIO("Test MappingFrontier maps checking") {
  Circuit circ(3);
  GIVEN("Valid maps") {
    std::shared_ptr<unit_bimaps_t> maps = std::make_shared<unit_bimaps_t>();
    maps->initial.insert({Qubit(0), Qubit(0)});
    maps->final.insert({Qubit(0), Qubit(0)});
    maps->initial.insert({Qubit(1), Qubit(1)});
    maps->final.insert({Qubit(1), Qubit(1)});
    maps->initial.insert({Qubit(2), Qubit(2)});
    maps->final.insert({Qubit(2), Qubit(2)});
    MappingFrontier mf(circ, maps);
  }
  GIVEN("Maps with wrong size") {
    std::shared_ptr<unit_bimaps_t> maps = std::make_shared<unit_bimaps_t>();
    maps->initial.insert({Qubit(0), Qubit(0)});
    maps->final.insert({Qubit(0), Qubit(0)});
    maps->initial.insert({Qubit(1), Qubit(1)});
    maps->final.insert({Qubit(1), Qubit(1)});
    REQUIRE_THROWS_AS(MappingFrontier(circ, maps), MappingFrontierError);
  }
  GIVEN("Uids not found in map") {
    std::shared_ptr<unit_bimaps_t> maps = std::make_shared<unit_bimaps_t>();
    maps->initial.insert({Qubit(0), Node(0)});
    maps->final.insert({Qubit(0), Qubit(0)});
    maps->initial.insert({Qubit(1), Qubit(1)});
    maps->final.insert({Qubit(1), Qubit(1)});
    maps->initial.insert({Qubit(2), Qubit(2)});
    maps->final.insert({Qubit(2), Qubit(2)});

    REQUIRE_THROWS_AS(MappingFrontier(circ, maps), MappingFrontierError);
  }
}

}  // namespace tket
