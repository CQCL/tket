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

#include "ArchAwareSynth/SteinerTree.hpp"

namespace tket {

SCENARIO("Build some basic trees") {
  GIVEN("1-vertex tree in 2-vertex graph") {
    MatrixXb connectivity(2, 2);
    connectivity << 0, 1, 1, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{0};
    aas::SteinerTree st(handler, nodes_to_add, 0);
    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf, aas::SteinerNodeType::OutOfTree};
    REQUIRE(st.node_types == correct_types);
    MatrixXb correct_neighbours(2, 2);
    correct_neighbours << 0, 0, 0, 0;
    REQUIRE(nodes_to_add.size() == 0);
  }
  GIVEN("2-vertex tree in 2-vertex graph") {
    MatrixXb connectivity(2, 2);
    connectivity << 0, 1, 1, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{0, 1};
    aas::SteinerTree st(handler, nodes_to_add, 0);
    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf, aas::SteinerNodeType::Leaf};
    REQUIRE(st.node_types == correct_types);
    REQUIRE(nodes_to_add.size() == 0);
  }
  GIVEN("2-vertex tree in 3-vertex graph") {
    MatrixXb connectivity(3, 3);
    connectivity << 0, 1, 0, 1, 0, 1, 0, 1, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{0, 2};
    aas::SteinerTree st(handler, nodes_to_add, 0);
    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf, aas::SteinerNodeType::ZeroInTree,
        aas::SteinerNodeType::Leaf};
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("Different 2-vertex tree in 3-vertex graph") {
    MatrixXb connectivity(3, 3);
    connectivity << 0, 1, 0, 1, 0, 1, 0, 1, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{0, 1};
    aas::SteinerTree st1(handler, nodes_to_add, 0);
    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf, aas::SteinerNodeType::Leaf,
        aas::SteinerNodeType::OutOfTree};
    REQUIRE(st1.node_types == correct_types);
    MatrixXb correct_neighbours(3, 3);
    correct_neighbours << 0, 1, 0, 1, 0, 0, 0, 0, 0;

    nodes_to_add = {0, 1};
    aas::SteinerTree st2(handler, nodes_to_add, 0);
    REQUIRE(st2.node_types == correct_types);
  }
  GIVEN("3-vertex tree in 4-vertex graph") {
    MatrixXb connectivity(4, 4);
    connectivity << 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{0, 2, 3};
    aas::SteinerTree st1(handler, nodes_to_add, 0);
    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf, aas::SteinerNodeType::ZeroInTree,
        aas::SteinerNodeType::OneInTree, aas::SteinerNodeType::Leaf};
    REQUIRE(st1.node_types == correct_types);
  }
  GIVEN("2-vertex tree in 6-vertex graph") {
    MatrixXb connectivity(6, 6);
    connectivity << 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{0, 2};
    aas::SteinerTree st(handler, nodes_to_add, 0);
    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,      aas::SteinerNodeType::OutOfTree,
        aas::SteinerNodeType::Leaf,      aas::SteinerNodeType::OutOfTree,
        aas::SteinerNodeType::OutOfTree, aas::SteinerNodeType::OutOfTree};
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("2-vertex tree in 6-vertex graph, check Zero in Tree") {
    MatrixXb connectivity(6, 6);
    connectivity << 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{1, 3};
    aas::SteinerTree st(handler, nodes_to_add, 0);
    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::ZeroInTree, aas::SteinerNodeType::Leaf,
        aas::SteinerNodeType::OutOfTree,  aas::SteinerNodeType::Leaf,
        aas::SteinerNodeType::OutOfTree,  aas::SteinerNodeType::OutOfTree};
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("2-vertex tree in 6-vertex graph, check OneInTree") {
    MatrixXb connectivity(6, 6);
    connectivity << 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{0, 1, 3};
    aas::SteinerTree st(handler, nodes_to_add, 0);
    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::OneInTree, aas::SteinerNodeType::Leaf,
        aas::SteinerNodeType::OutOfTree, aas::SteinerNodeType::Leaf,
        aas::SteinerNodeType::OutOfTree, aas::SteinerNodeType::OutOfTree};
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("2-vertex tree in 6-vertex graph, check OneInTree, unsorted") {
    MatrixXb connectivity(6, 6);
    connectivity << 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{3, 0, 1};
    aas::SteinerTree st(handler, nodes_to_add, 0);
    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::OneInTree, aas::SteinerNodeType::Leaf,
        aas::SteinerNodeType::OutOfTree, aas::SteinerNodeType::Leaf,
        aas::SteinerNodeType::OutOfTree, aas::SteinerNodeType::OutOfTree};
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("2-vertex tree in 6-vertex graph, cyclic") {
    MatrixXb connectivity(6, 6);
    connectivity << 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0,
        0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{3, 0, 1};
    aas::SteinerTree st(handler, nodes_to_add, 0);
    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::OneInTree, aas::SteinerNodeType::Leaf,
        aas::SteinerNodeType::OutOfTree, aas::SteinerNodeType::Leaf,
        aas::SteinerNodeType::OutOfTree, aas::SteinerNodeType::OutOfTree};
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("2-vertex tree in 6-vertex graph, acyclic") {
    MatrixXb connectivity(2, 2);
    connectivity << 0, 1, 1, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{0, 1};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf, aas::SteinerNodeType::Leaf};
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("2-vertex tree in 6-vertex graph, acyclic") {
    MatrixXb connectivity(2, 2);
    connectivity << 0, 1, 1, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{0, 1};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf, aas::SteinerNodeType::Leaf};
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("2-vertex tree in 6-vertex graph, acyclic") {
    MatrixXb connectivity(3, 3);
    connectivity << 0, 1, 1, 1, 0, 0, 1, 0, 0;
    aas::PathHandler handler(connectivity);
    std::list<unsigned> nodes_to_add{0, 1, 2};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::OneInTree, aas::SteinerNodeType::Leaf,
        aas::SteinerNodeType::Leaf};
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("complex architecture") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(2)},
         {Node(2), Node(3)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 2};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::OneInTree, aas::SteinerNodeType::Leaf,
        aas::SteinerNodeType::Leaf, aas::SteinerNodeType::OutOfTree};
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("complex architecture II") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(0), Node(4)},
         {Node(1), Node(5)},
         {Node(5), Node(6)},
         {Node(6), Node(7)},
         {Node(7), Node(8)},
         {Node(8), Node(9)},
         {Node(9), Node(10)},
         {Node(10), Node(11)},
         {Node(11), Node(2)},
         {Node(2), Node(3)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 2, 11, 8};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::OneInTree,   // 0
        aas::SteinerNodeType::Leaf,        // 1
        aas::SteinerNodeType::OneInTree,   // 2
        aas::SteinerNodeType::OutOfTree,   // 3
        aas::SteinerNodeType::OutOfTree,   // 4
        aas::SteinerNodeType::OutOfTree,   // 5
        aas::SteinerNodeType::OutOfTree,   // 6
        aas::SteinerNodeType::OutOfTree,   // 7
        aas::SteinerNodeType::Leaf,        // 8
        aas::SteinerNodeType::ZeroInTree,  // 9
        aas::SteinerNodeType::ZeroInTree,  // 10
        aas::SteinerNodeType::OneInTree};  // 11
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("complex architecture III") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(0), Node(4)},
         {Node(1), Node(5)},
         {Node(5), Node(6)},
         {Node(6), Node(7)},
         {Node(7), Node(8)},
         {Node(0), Node(8)},
         {Node(8), Node(9)},
         {Node(9), Node(10)},
         {Node(10), Node(11)},
         {Node(11), Node(2)},
         {Node(2), Node(3)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 2, 11, 8};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::OneInTree,  // 0
        aas::SteinerNodeType::Leaf,       // 1
        aas::SteinerNodeType::OneInTree,  // 2
        aas::SteinerNodeType::OutOfTree,  // 3
        aas::SteinerNodeType::OutOfTree,  // 4
        aas::SteinerNodeType::OutOfTree,  // 5
        aas::SteinerNodeType::OutOfTree,  // 6
        aas::SteinerNodeType::OutOfTree,  // 7
        aas::SteinerNodeType::Leaf,       // 8
        aas::SteinerNodeType::OutOfTree,  // 9
        aas::SteinerNodeType::OutOfTree,  // 10
        aas::SteinerNodeType::Leaf};      // 11
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("binar tree") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(2), Node(6)},
         {Node(3), Node(7)},
         {Node(3), Node(8)},
         {Node(4), Node(9)},
         {Node(4), Node(10)},
         {Node(5), Node(11)},
         {Node(5), Node(12)},
         {Node(6), Node(13)},
         {Node(6), Node(14)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 2, 11, 8};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::OneInTree,   // 0
        aas::SteinerNodeType::OneInTree,   // 1
        aas::SteinerNodeType::OneInTree,   // 2
        aas::SteinerNodeType::ZeroInTree,  // 3
        aas::SteinerNodeType::OutOfTree,   // 4
        aas::SteinerNodeType::ZeroInTree,  // 5
        aas::SteinerNodeType::OutOfTree,   // 6
        aas::SteinerNodeType::OutOfTree,   // 7
        aas::SteinerNodeType::Leaf,        // 8
        aas::SteinerNodeType::OutOfTree,   // 9
        aas::SteinerNodeType::OutOfTree,   // 10
        aas::SteinerNodeType::Leaf,        // 11
        aas::SteinerNodeType::OutOfTree,   // 12
        aas::SteinerNodeType::OutOfTree,   // 13
        aas::SteinerNodeType::OutOfTree};  // 14
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("binar tree II") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(2), Node(6)},
         {Node(3), Node(7)},
         {Node(3), Node(8)},
         {Node(4), Node(9)},
         {Node(4), Node(10)},
         {Node(5), Node(11)},
         {Node(5), Node(12)},
         {Node(6), Node(13)},
         {Node(6), Node(14)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 2, 11, 8, 14};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::OneInTree,   // 0
        aas::SteinerNodeType::OneInTree,   // 1
        aas::SteinerNodeType::OneInTree,   // 2
        aas::SteinerNodeType::ZeroInTree,  // 3
        aas::SteinerNodeType::OutOfTree,   // 4
        aas::SteinerNodeType::ZeroInTree,  // 5
        aas::SteinerNodeType::ZeroInTree,  // 6
        aas::SteinerNodeType::OutOfTree,   // 7
        aas::SteinerNodeType::Leaf,        // 8
        aas::SteinerNodeType::OutOfTree,   // 9
        aas::SteinerNodeType::OutOfTree,   // 10
        aas::SteinerNodeType::Leaf,        // 11
        aas::SteinerNodeType::OutOfTree,   // 12
        aas::SteinerNodeType::OutOfTree,   // 13
        aas::SteinerNodeType::Leaf};       // 14
    REQUIRE(st.node_types == correct_types);
  }

  GIVEN("operations available") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(2), Node(6)},
         {Node(3), Node(7)},
         {Node(3), Node(8)},
         {Node(4), Node(9)},
         {Node(4), Node(10)},
         {Node(5), Node(11)},
         {Node(5), Node(12)},
         {Node(6), Node(13)},
         {Node(6), Node(14)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,        // 0
        aas::SteinerNodeType::Leaf,        // 1
        aas::SteinerNodeType::OutOfTree,   // 2
        aas::SteinerNodeType::OutOfTree,   // 3
        aas::SteinerNodeType::OutOfTree,   // 4
        aas::SteinerNodeType::OutOfTree,   // 5
        aas::SteinerNodeType::OutOfTree,   // 6
        aas::SteinerNodeType::OutOfTree,   // 7
        aas::SteinerNodeType::OutOfTree,   // 8
        aas::SteinerNodeType::OutOfTree,   // 9
        aas::SteinerNodeType::OutOfTree,   // 10
        aas::SteinerNodeType::OutOfTree,   // 11
        aas::SteinerNodeType::OutOfTree,   // 12
        aas::SteinerNodeType::OutOfTree,   // 13
        aas::SteinerNodeType::OutOfTree};  // 14
    REQUIRE(st.node_types == correct_types);

    aas::PathHandler handler2(archi);

    aas::OperationList res2 = st.operations_available(handler2);

    aas::OperationList res = {std::pair(0, 1), std::pair(1, 0)};

    REQUIRE(res2 == res);
  }

  GIVEN("cost of operation") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(2), Node(6)},
         {Node(3), Node(7)},
         {Node(3), Node(8)},
         {Node(4), Node(9)},
         {Node(4), Node(10)},
         {Node(5), Node(11)},
         {Node(5), Node(12)},
         {Node(6), Node(13)},
         {Node(6), Node(14)},
         {Node(7), Node(15)},
         {Node(7), Node(16)},
         {Node(8), Node(17)},
         {Node(8), Node(18)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 11, 8, 14, 5};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::OneInTree,   // 0
        aas::SteinerNodeType::OneInTree,   // 1
        aas::SteinerNodeType::ZeroInTree,  // 2
        aas::SteinerNodeType::ZeroInTree,  // 3
        aas::SteinerNodeType::OutOfTree,   // 4
        aas::SteinerNodeType::OneInTree,   // 5
        aas::SteinerNodeType::ZeroInTree,  // 6
        aas::SteinerNodeType::OutOfTree,   // 7
        aas::SteinerNodeType::Leaf,        // 8
        aas::SteinerNodeType::OutOfTree,   // 9
        aas::SteinerNodeType::OutOfTree,   // 10
        aas::SteinerNodeType::Leaf,        // 11
        aas::SteinerNodeType::OutOfTree,   // 12
        aas::SteinerNodeType::OutOfTree,   // 13
        aas::SteinerNodeType::Leaf,        // 14
        aas::SteinerNodeType::OutOfTree,   // 15
        aas::SteinerNodeType::OutOfTree,   // 16
        aas::SteinerNodeType::OutOfTree,   // 17
        aas::SteinerNodeType::OutOfTree};  // 18
    REQUIRE(st.node_types == correct_types);
    REQUIRE(st.cost_of_operation(2, 3) == 0);
    REQUIRE(st.cost_of_operation(3, 2) == 0);
    REQUIRE(st.cost_of_operation(3, 7) == 0);
    REQUIRE(st.node_types == correct_types);

    REQUIRE(st.cost_of_operation(1, 3) == -1);
    REQUIRE(st.cost_of_operation(0, 1) == 1);
    REQUIRE(st.cost_of_operation(1, 4) == 1);
    REQUIRE(st.cost_of_operation(5, 11) == -1);
    REQUIRE(st.node_types == correct_types);

    REQUIRE(st.cost_of_operation(11, 5) == 1);
    REQUIRE(st.cost_of_operation(8, 3) == -1);
    REQUIRE(st.cost_of_operation(14, 6) == -1);
    REQUIRE(st.cost_of_operation(8, 17) == 1);
    REQUIRE(st.cost_of_operation(8, 18) == 1);
    REQUIRE(st.cost_of_operation(8, 14) == -1);  // not connected
    REQUIRE(st.node_types == correct_types);

    REQUIRE(st.cost_of_operation(4, 9) == 0);
    REQUIRE(st.cost_of_operation(4, 10) == 0);
    REQUIRE(st.cost_of_operation(4, 1) == 0);
    REQUIRE(st.cost_of_operation(7, 3) == 0);
    REQUIRE(st.cost_of_operation(7, 15) == 0);
    REQUIRE(st.cost_of_operation(7, 16) == 0);
    REQUIRE(st.cost_of_operation(18, 8) == 0);
    REQUIRE(st.node_types == correct_types);
  }

  GIVEN("cost of operation II") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 4};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,       // 0
        aas::SteinerNodeType::OneInTree,  // 1
        aas::SteinerNodeType::OutOfTree,  // 2
        aas::SteinerNodeType::OutOfTree,  // 3
        aas::SteinerNodeType::Leaf};      // 4
    REQUIRE(st.node_types == correct_types);
    REQUIRE(st.cost_of_operation(0, 4) == -1);
    REQUIRE(st.cost_of_operation(4, 0) == -1);
    REQUIRE(st.cost_of_operation(0, 1) == 1);
    REQUIRE(st.cost_of_operation(1, 0) == -1);
    REQUIRE(st.cost_of_operation(0, 2) == 1);  // not connected
    REQUIRE(st.cost_of_operation(2, 0) == 0);  // not connected
    REQUIRE(st.node_types == correct_types);
  }
  GIVEN("add row") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(2), Node(6)},
         {Node(3), Node(7)},
         {Node(3), Node(8)},
         {Node(4), Node(9)},
         {Node(4), Node(10)},
         {Node(5), Node(11)},
         {Node(5), Node(12)},
         {Node(6), Node(13)},
         {Node(6), Node(14)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 2, 11, 8, 14};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::OneInTree,   // 0
        aas::SteinerNodeType::OneInTree,   // 1
        aas::SteinerNodeType::OneInTree,   // 2
        aas::SteinerNodeType::ZeroInTree,  // 3
        aas::SteinerNodeType::OutOfTree,   // 4
        aas::SteinerNodeType::ZeroInTree,  // 5
        aas::SteinerNodeType::ZeroInTree,  // 6
        aas::SteinerNodeType::OutOfTree,   // 7
        aas::SteinerNodeType::Leaf,        // 8
        aas::SteinerNodeType::OutOfTree,   // 9
        aas::SteinerNodeType::OutOfTree,   // 10
        aas::SteinerNodeType::Leaf,        // 11
        aas::SteinerNodeType::OutOfTree,   // 12
        aas::SteinerNodeType::OutOfTree,   // 13
        aas::SteinerNodeType::Leaf};       // 14
    REQUIRE(st.node_types == correct_types);
    st.add_row(5, 6);
    st.add_row(0, 1);
    st.add_row(0, 3);
    st.add_row(0, 14);
    st.add_row(0, 6);
    st.add_row(7, 8);
    st.add_row(9, 10);
    st.add_row(1, 2);
    st.add_row(13, 14);
    st.add_row(14, 13);
    st.add_row(10, 6);
    st.add_row(10, 1);
    std::vector<aas::SteinerNodeType> correct_types_2{
        aas::SteinerNodeType::Leaf,        // 0
        aas::SteinerNodeType::ZeroInTree,  // 1
        aas::SteinerNodeType::OneInTree,   // 2
        aas::SteinerNodeType::OneInTree,   // 3
        aas::SteinerNodeType::OutOfTree,   // 4
        aas::SteinerNodeType::ZeroInTree,  // 5
        aas::SteinerNodeType::OneInTree,   // 6
        aas::SteinerNodeType::OutOfTree,   // 7
        aas::SteinerNodeType::Leaf,        // 8
        aas::SteinerNodeType::OutOfTree,   // 9
        aas::SteinerNodeType::OutOfTree,   // 10
        aas::SteinerNodeType::Leaf,        // 11
        aas::SteinerNodeType::OutOfTree,   // 12
        aas::SteinerNodeType::OutOfTree,   // 13
        aas::SteinerNodeType::OutOfTree};  // 14
    REQUIRE(st.node_types == correct_types_2);
  }
  GIVEN("add row II") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 4};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,       // 0
        aas::SteinerNodeType::OneInTree,  // 1
        aas::SteinerNodeType::OutOfTree,  // 2
        aas::SteinerNodeType::OutOfTree,  // 3
        aas::SteinerNodeType::Leaf};      // 4
    REQUIRE(st.node_types == correct_types);
    st.add_row(0, 1);
    st.add_row(0, 2);
    st.add_row(0, 3);
    st.add_row(0, 4);
  }
  GIVEN("add row III") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 4};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,       // 0
        aas::SteinerNodeType::OneInTree,  // 1
        aas::SteinerNodeType::OutOfTree,  // 2
        aas::SteinerNodeType::OutOfTree,  // 3
        aas::SteinerNodeType::Leaf};      // 4
    REQUIRE(st.node_types == correct_types);
    st.add_row(2, 0);
    REQUIRE(st.node_types == correct_types);
    st.add_row(1, 2);
    st.add_row(2, 1);
  }
  GIVEN("add row IV") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 4};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,       // 0
        aas::SteinerNodeType::OneInTree,  // 1
        aas::SteinerNodeType::OutOfTree,  // 2
        aas::SteinerNodeType::OutOfTree,  // 3
        aas::SteinerNodeType::Leaf};      // 4
    REQUIRE(st.node_types == correct_types);
    st.add_row(1, 4);
    std::vector<aas::SteinerNodeType> correct_types2{
        aas::SteinerNodeType::Leaf,        // 0
        aas::SteinerNodeType::Leaf,        // 1
        aas::SteinerNodeType::OutOfTree,   // 2
        aas::SteinerNodeType::OutOfTree,   // 3
        aas::SteinerNodeType::OutOfTree};  // 4
    REQUIRE(st.node_types == correct_types2);
  }
  GIVEN("add row IVb") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 4};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,       // 0
        aas::SteinerNodeType::OneInTree,  // 1
        aas::SteinerNodeType::OutOfTree,  // 2
        aas::SteinerNodeType::OutOfTree,  // 3
        aas::SteinerNodeType::Leaf};      // 4
    REQUIRE(st.node_types == correct_types);
    st.add_row(4, 1);
    std::vector<aas::SteinerNodeType> correct_types_2{
        aas::SteinerNodeType::Leaf,        // 0
        aas::SteinerNodeType::ZeroInTree,  // 1
        aas::SteinerNodeType::OutOfTree,   // 2
        aas::SteinerNodeType::OutOfTree,   // 3
        aas::SteinerNodeType::Leaf};       // 4
    REQUIRE(st.node_types == correct_types_2);
  }
  GIVEN("add row V") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 4};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,       // 0
        aas::SteinerNodeType::OneInTree,  // 1
        aas::SteinerNodeType::OutOfTree,  // 2
        aas::SteinerNodeType::OutOfTree,  // 3
        aas::SteinerNodeType::Leaf};      // 4
    REQUIRE(st.node_types == correct_types);
    st.add_row(1, 0);
    std::vector<aas::SteinerNodeType> correct_types_2{
        aas::SteinerNodeType::OutOfTree,  // 0
        aas::SteinerNodeType::Leaf,       // 1
        aas::SteinerNodeType::OutOfTree,  // 2
        aas::SteinerNodeType::OutOfTree,  // 3
        aas::SteinerNodeType::Leaf};      // 4
    REQUIRE(st.node_types == correct_types_2);
  }
  GIVEN("add row VI") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 4};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,       // 0
        aas::SteinerNodeType::OneInTree,  // 1
        aas::SteinerNodeType::OutOfTree,  // 2
        aas::SteinerNodeType::OutOfTree,  // 3
        aas::SteinerNodeType::Leaf};      // 4
    REQUIRE(st.node_types == correct_types);
    st.add_row(0, 1);
    std::vector<aas::SteinerNodeType> correct_types_2{
        aas::SteinerNodeType::Leaf,        // 0
        aas::SteinerNodeType::ZeroInTree,  // 1
        aas::SteinerNodeType::OutOfTree,   // 2
        aas::SteinerNodeType::OutOfTree,   // 3
        aas::SteinerNodeType::Leaf};       // 4
    REQUIRE(st.node_types == correct_types_2);
  }
  GIVEN("add row VII") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1, 4};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,       // 0
        aas::SteinerNodeType::OneInTree,  // 1
        aas::SteinerNodeType::OutOfTree,  // 2
        aas::SteinerNodeType::OutOfTree,  // 3
        aas::SteinerNodeType::Leaf};      // 4
    REQUIRE(st.node_types == correct_types);
    st.add_row(0, 1);
    std::vector<aas::SteinerNodeType> correct_types_2{
        aas::SteinerNodeType::Leaf,        // 0
        aas::SteinerNodeType::ZeroInTree,  // 1
        aas::SteinerNodeType::OutOfTree,   // 2
        aas::SteinerNodeType::OutOfTree,   // 3
        aas::SteinerNodeType::Leaf};       // 4
    REQUIRE(st.node_types == correct_types_2);
  }
  GIVEN("add row VIII") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 1};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,        // 0
        aas::SteinerNodeType::Leaf,        // 1
        aas::SteinerNodeType::OutOfTree,   // 2
        aas::SteinerNodeType::OutOfTree,   // 3
        aas::SteinerNodeType::OutOfTree};  // 4
    REQUIRE(st.node_types == correct_types);
    st.add_row(1, 0);
    std::vector<aas::SteinerNodeType> correct_types_2{
        aas::SteinerNodeType::OutOfTree,   // 0
        aas::SteinerNodeType::OutOfTree,   // 1
        aas::SteinerNodeType::OutOfTree,   // 2
        aas::SteinerNodeType::OutOfTree,   // 3
        aas::SteinerNodeType::OutOfTree};  // 4
    REQUIRE(st.node_types == correct_types_2);
  }
  GIVEN("nodes") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 2};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,        // 0
        aas::SteinerNodeType::ZeroInTree,  // 1
        aas::SteinerNodeType::Leaf,        // 2
        aas::SteinerNodeType::OutOfTree,   // 3
        aas::SteinerNodeType::OutOfTree};  // 4
    REQUIRE(st.node_types == correct_types);
    std::vector<unsigned> all_nodes = st.nodes();

    std::vector<unsigned> expected_all_nodes = {0, 1, 2};

    REQUIRE(all_nodes == expected_all_nodes);
  }
  GIVEN("nodes") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(1), Node(3)},
         {Node(1), Node(4)}});

    aas::PathHandler handler(archi);
    std::list<unsigned> nodes_to_add{0, 3};
    aas::SteinerTree st(handler, nodes_to_add, 0);

    std::vector<aas::SteinerNodeType> correct_types{
        aas::SteinerNodeType::Leaf,        // 0
        aas::SteinerNodeType::ZeroInTree,  // 1
        aas::SteinerNodeType::OutOfTree,   // 2
        aas::SteinerNodeType::Leaf,        // 3
        aas::SteinerNodeType::OutOfTree};  // 4
    REQUIRE(st.node_types == correct_types);
    std::vector<unsigned> all_nodes = st.nodes();

    std::vector<unsigned> expected_all_nodes = {0, 1, 3};

    REQUIRE(all_nodes == expected_all_nodes);
  }
  GIVEN("test swap cnot synth - 1") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(3), Node(6)},
         {Node(4), Node(7)},
         {Node(5), Node(8)},
         {Node(6), Node(9)}});

    aas::PathHandler handler(archi);

    MatrixXb matrix = MatrixXb::Identity(10, 10);
    DiagMatrix CNOT_matrix = DiagMatrix(matrix);

    aas::CNotSwapSynth cnot(handler, CNOT_matrix);
    Circuit result = cnot.get_circuit();
    REQUIRE(cnot.valid_result());
  }
  GIVEN("test swap cnot synth - 2") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(3), Node(6)},
         {Node(4), Node(7)},
         {Node(5), Node(8)},
         {Node(6), Node(9)}});

    aas::PathHandler handler(archi);

    MatrixXb mat = MatrixXb::Identity(10, 10);
    mat << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  // 0
        0, 1, 1, 1, 1, 1, 1, 1, 1, 1,     // 1
        0, 0, 1, 1, 1, 1, 1, 1, 1, 1,     // 2
        0, 0, 0, 1, 1, 1, 1, 1, 1, 1,     // 3
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1,     // 4
        0, 0, 0, 0, 0, 1, 1, 1, 1, 1,     // 5
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1,     // 6
        0, 0, 0, 0, 0, 0, 0, 1, 1, 1,     // 7
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1,     // 8
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1;     // 9

    DiagMatrix CNOT_matrix = DiagMatrix(mat);

    aas::CNotSwapSynth cnot(handler, CNOT_matrix);
    Circuit result = cnot.get_circuit();
    REQUIRE(cnot.valid_result());
  }
}
}  // namespace tket
