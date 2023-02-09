// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "ArchAwareSynth/Path.hpp"
#include "ArchAwareSynth/SteinerTree.hpp"

namespace tket {

SCENARIO("Check graph construction is correct") {
  GIVEN("2-vertex graph") {
    MatrixXb connectivity(2, 2);
    connectivity << 0, 1,  // 0
        1, 0;              // 1
    aas::PathHandler handler1(connectivity);
    aas::MatrixXu correct_distance_matrix(2, 2);
    correct_distance_matrix << 0, 1, 1, 0;
    REQUIRE(handler1.get_distance_matrix() == correct_distance_matrix);

    aas::MatrixXu correct_path_matrix(2, 2);
    correct_path_matrix << 0, 1, 0, 1;
    REQUIRE(handler1.get_path_matrix() == correct_path_matrix);

    std::list<unsigned> path1 = handler1.find_path(0, 1);
    REQUIRE(path1.size() == 2);
    REQUIRE(*++path1.begin() == 1);
    std::list<unsigned> path2 = handler1.find_path(1, 0);
    REQUIRE(path2.size() == 2);
    REQUIRE(*++path2.begin() == 0);
  }
  GIVEN("3-vertex graph") {
    MatrixXb connectivity(3, 3);
    connectivity << 0, 1, 0,  // 0
        1, 0, 1,              // 1
        0, 1, 0;              // 2
    aas::PathHandler handler2(connectivity);
    aas::MatrixXu correct_distance_matrix(3, 3);
    correct_distance_matrix << 0, 1, 2, 1, 0, 1, 2, 1, 0;
    REQUIRE(handler2.get_distance_matrix() == correct_distance_matrix);

    aas::MatrixXu correct_path_matrix(3, 3);
    correct_path_matrix << 0, 1, 1, 0, 1, 2, 1, 1, 2;
    REQUIRE(handler2.get_path_matrix() == correct_path_matrix);

    std::list<unsigned> path1 = handler2.find_path(0, 2);
    std::list<unsigned> correct_path1{0, 1, 2};
    REQUIRE(path1 == correct_path1);
  }
  GIVEN("4-vertex graph") {
    MatrixXb connectivity(4, 4);
    connectivity << 0, 1, 1, 0,  // 0
        1, 0, 1, 0,              // 1
        1, 1, 0, 1,              // 2
        0, 0, 1, 0;              // 3
    aas::PathHandler handler3(connectivity);
    aas::MatrixXu correct_distance_matrix(4, 4);
    correct_distance_matrix << 0, 1, 1, 2,  // 0
        1, 0, 1, 2,                         // 1
        1, 1, 0, 1,                         // 2
        2, 2, 1, 0;
    REQUIRE(handler3.get_distance_matrix() == correct_distance_matrix);
    aas::MatrixXu correct_path_matrix(4, 4);
    correct_path_matrix << 0, 1, 2, 2,  // 0
        0, 1, 2, 2,                     // 1
        0, 1, 2, 3,                     // 2
        2, 2, 2, 3;                     // 3
    REQUIRE(handler3.get_path_matrix() == correct_path_matrix);
  }
  GIVEN("3-vertex graph, unreachable node") {
    MatrixXb connectivity(3, 3);
    connectivity << 0, 1, 0,  // 0
        1, 0, 0,              // 1
        0, 0, 0;              // 2

    aas::PathHandler handler;
    handler = aas::PathHandler(connectivity);

    aas::MatrixXu correct_distance_matrix(3, 3);
    REQUIRE(handler.get_distance_matrix()(0, 0) == 0);
    REQUIRE(handler.get_distance_matrix()(0, 1) == 1);
    REQUIRE(handler.get_distance_matrix()(0, 2) >= 3);
    REQUIRE(handler.get_distance_matrix()(1, 0) == 1);
    REQUIRE(handler.get_distance_matrix()(1, 1) == 0);
    REQUIRE(handler.get_distance_matrix()(1, 2) >= 3);
    REQUIRE(handler.get_distance_matrix()(2, 0) >= 3);
    REQUIRE(handler.get_distance_matrix()(2, 1) >= 3);
    REQUIRE(handler.get_distance_matrix()(2, 2) == 0);

    aas::MatrixXu correct_path_matrix(3, 3);
    correct_path_matrix << 0, 1, 3,  // 0
        0, 1, 3,                     // 1
        3, 3, 2;                     // 2
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);
  }
  GIVEN("4-vertex graph") {
    MatrixXb connectivity(4, 4);
    connectivity << 0, 1, 1, 0,  // 0
        1, 0, 1, 0,              // 1
        1, 1, 0, 1,              // 2
        0, 0, 1, 0;              // 3
    aas::PathHandler handler(connectivity);
    aas::MatrixXu correct_distance_matrix(4, 4);
    correct_distance_matrix << 0, 1, 1, 2,  // 0
        1, 0, 1, 2,                         // 1
        1, 1, 0, 1,                         // 2
        2, 2, 1, 0;                         // 3
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);
    aas::MatrixXu correct_path_matrix(4, 4);
    correct_path_matrix << 0, 1, 2, 2,  // 0
        0, 1, 2, 2,                     // 1
        0, 1, 2, 3,                     // 2
        2, 2, 2, 3;                     // 3
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);
  }
  GIVEN("wrong 4-vertex graph") {
    MatrixXb connectivity(4, 4);
    connectivity << 1, 1, 1, 0,  // 0
        0, 0, 1, 0,              // 1
        1, 1, 0, 1,              // 2
        0, 0, 1, 0;              // 3
    aas::PathHandler handler(connectivity);

    aas::MatrixXu correct_distance_matrix(4, 4);
    correct_distance_matrix << 0, 1, 1, 2,  // 0
        2, 0, 1, 2,                         // 1
        1, 1, 0, 1,                         // 2
        2, 2, 1, 0;                         // 3
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);
    aas::MatrixXu correct_path_matrix(4, 4);
    correct_path_matrix << 0, 1, 2, 2,  // 0
        2, 1, 2, 2,                     // 1
        0, 1, 2, 3,                     // 2
        2, 2, 2, 3;                     // 3
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);
  }
  GIVEN("wrong 4-vertex graph complet graph") {
    MatrixXb connectivity(4, 4);
    connectivity << 0, 1, 1, 1,  // 0
        1, 0, 1, 1,              // 1
        1, 1, 0, 1,              // 2
        1, 1, 1, 0;              // 3
    aas::PathHandler handler(connectivity);

    aas::MatrixXu correct_distance_matrix(4, 4);
    correct_distance_matrix << 0, 1, 1, 1,  // 0
        1, 0, 1, 1,                         // 1
        1, 1, 0, 1,                         // 2
        1, 1, 1, 0;                         // 3
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);
    aas::MatrixXu correct_path_matrix(4, 4);
    correct_path_matrix << 0, 1, 2, 3,  // 0
        0, 1, 2, 3,                     // 1
        0, 1, 2, 3,                     // 2
        0, 1, 2, 3;                     // 3
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);
  }
  GIVEN("6-vertex graph") {
    MatrixXb connectivity(6, 6);
    connectivity << 0, 1, 1, 1, 0, 1,  // 0
        1, 0, 1, 1, 1, 1,              // 1
        1, 1, 0, 0, 1, 1,              // 2
        1, 1, 0, 0, 0, 1,              // 3
        0, 1, 1, 0, 0, 1,              // 4
        1, 1, 1, 1, 1, 0;              // 5

    aas::PathHandler handler(connectivity);
    aas::MatrixXu correct_distance_matrix(6, 6);
    correct_distance_matrix << 0, 1, 1, 1, 2, 1,  // 0
        1, 0, 1, 1, 1, 1,                         // 1
        1, 1, 0, 2, 1, 1,                         // 2
        1, 1, 2, 0, 2, 1,                         // 3
        2, 1, 1, 2, 0, 1,                         // 4
        1, 1, 1, 1, 1, 0;                         // 5
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);
    aas::MatrixXu correct_path_matrix(6, 6);
    correct_path_matrix << 0, 1, 2, 3, 1, 5,  // 0
        0, 1, 2, 3, 4, 5,                     // 1
        0, 1, 2, 0, 4, 5,                     // 2
        0, 1, 0, 3, 1, 5,                     // 3
        1, 1, 2, 1, 4, 5,                     // 4
        0, 1, 2, 3, 4, 5;                     // 5
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);
  }
  GIVEN("6-vertex spare graph") {
    MatrixXb connectivity(6, 6);
    connectivity << 0, 1, 1, 1, 1, 1,  // 0
        1, 0, 0, 0, 0, 0,              // 1
        1, 0, 0, 0, 0, 0,              // 2
        1, 0, 0, 0, 0, 0,              // 3
        1, 0, 0, 0, 0, 0,              // 4
        1, 0, 0, 0, 0, 0;              // 5

    aas::PathHandler handler(connectivity);
    aas::MatrixXu correct_distance_matrix(6, 6);
    correct_distance_matrix << 0, 1, 1, 1, 1, 1,  // 0
        1, 0, 2, 2, 2, 2,                         // 1
        1, 2, 0, 2, 2, 2,                         // 2
        1, 2, 2, 0, 2, 2,                         // 3
        1, 2, 2, 2, 0, 2,                         // 4
        1, 2, 2, 2, 2, 0;                         // 5
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);
    aas::MatrixXu correct_path_matrix(6, 6);
    correct_path_matrix << 0, 1, 2, 3, 4, 5,  // 0
        0, 1, 0, 0, 0, 0,                     // 1
        0, 0, 2, 0, 0, 0,                     // 2
        0, 0, 0, 3, 0, 0,                     // 3
        0, 0, 0, 0, 4, 0,                     // 4
        0, 0, 0, 0, 0, 5;                     // 5
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);
  }
  GIVEN("4-vertex graph given as architecture") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(2)},
         {Node(2), Node(3)}});

    aas::PathHandler handler(archi);
    aas::MatrixXu correct_distance_matrix(4, 4);
    correct_distance_matrix << 0, 1, 1, 2,  // 0
        1, 0, 1, 2,                         // 1
        1, 1, 0, 1,                         // 2
        2, 2, 1, 0;                         // 3
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);
    aas::MatrixXu correct_path_matrix(4, 4);
    correct_path_matrix << 0, 1, 2, 2,  // 0
        0, 1, 2, 2,                     // 1
        0, 1, 2, 3,                     // 2
        2, 2, 2, 3;                     // 3
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);
  }

  GIVEN("acyclic path generation") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(1), Node(2)},
         {Node(2), Node(3)}});

    aas::PathHandler handler(archi);
    aas::MatrixXu correct_distance_matrix(4, 4);
    correct_distance_matrix << 0, 1, 1, 2,  // 0
        1, 0, 1, 2,                         // 1
        1, 1, 0, 1,                         // 2
        2, 2, 1, 0;                         // 3
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);
    aas::MatrixXu correct_path_matrix(4, 4);
    correct_path_matrix << 0, 1, 2, 2,  // 0
        0, 1, 2, 2,                     // 1
        0, 1, 2, 3,                     // 2
        2, 2, 2, 3;                     // 3
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);

    aas::PathHandler handi = handler.construct_acyclic_handler();

    aas::MatrixXu correct_distance_matrix_2(4, 4);
    correct_distance_matrix_2 << 0, 2, 1, 2,  // 0
        2, 0, 1, 2,                           // 1
        1, 1, 0, 1,                           // 2
        2, 2, 1, 0;                           // 3
    REQUIRE(handi.get_distance_matrix() == correct_distance_matrix_2);
    aas::MatrixXu correct_path_matrix_2(4, 4);
    correct_path_matrix_2 << 0, 2, 2, 2,  // 0
        2, 1, 2, 2,                       // 1
        0, 1, 2, 3,                       // 2
        2, 2, 2, 3;                       // 3
    REQUIRE(handi.get_path_matrix() == correct_path_matrix_2);
  }
  GIVEN("acyclic path generation II") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(0)}});

    aas::PathHandler handler(archi);

    aas::MatrixXu correct_distance_matrix(4, 4);
    correct_distance_matrix << 0, 1, 2, 1,  // 0
        1, 0, 1, 2,                         // 1
        2, 1, 0, 1,                         // 2
        1, 2, 1, 0;                         // 3
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);

    aas::MatrixXu correct_path_matrix(4, 4);
    correct_path_matrix << 0, 1, 1, 3,  // 0
        0, 1, 2, 0,                     // 1
        1, 1, 2, 3,                     // 2
        0, 0, 2, 3;                     // 3
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);

    aas::PathHandler handi = handler.construct_acyclic_handler();

    aas::MatrixXu correct_distance_matrix_2(4, 4);
    correct_distance_matrix_2 << 0, 1, 2, 1,  // 0
        1, 0, 1, 2,                           // 1
        2, 1, 0, 3,                           // 2
        1, 2, 3, 0;                           // 3
    REQUIRE(handi.get_distance_matrix() == correct_distance_matrix_2);

    aas::MatrixXu correct_path_matrix_2(4, 4);
    correct_path_matrix_2 << 0, 1, 1, 3,  // 0
        0, 1, 2, 0,                       // 1
        1, 1, 2, 1,                       // 2
        0, 0, 0, 3;                       // 3
    REQUIRE(handi.get_path_matrix() == correct_path_matrix_2);

    aas::PathHandler handi2 = handi.construct_acyclic_handler();

    aas::MatrixXu correct_distance_matrix_3(4, 4);
    correct_distance_matrix_3 << 0, 1, 2, 1,  // 0
        1, 0, 1, 2,                           // 1
        2, 1, 0, 3,                           // 2
        1, 2, 3, 0;                           // 3
    REQUIRE(handi2.get_distance_matrix() == correct_distance_matrix_3);

    aas::MatrixXu correct_path_matrix_3(4, 4);
    correct_path_matrix_3 << 0, 1, 1, 3,  // 0
        0, 1, 2, 0,                       // 1
        1, 1, 2, 1,                       // 2
        0, 0, 0, 3;                       // 3
    REQUIRE(handi2.get_path_matrix() == correct_path_matrix_3);
  }
  GIVEN("acyclic path generation III") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(0)}});

    aas::PathHandler handler(archi);

    aas::MatrixXu correct_distance_matrix(4, 4);
    correct_distance_matrix << 0, 1, 2, 1,  // 0
        1, 0, 1, 2,                         // 1
        2, 1, 0, 1,                         // 2
        1, 2, 1, 0;                         // 3
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);

    aas::MatrixXu correct_path_matrix(4, 4);
    correct_path_matrix << 0, 1, 1, 3,  // 0
        0, 1, 2, 0,                     // 1
        1, 1, 2, 3,                     // 2
        0, 0, 2, 3;                     // 3
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);

    aas::PathHandler handi = handler.construct_acyclic_handler();

    aas::MatrixXu correct_distance_matrix_2(4, 4);
    correct_distance_matrix_2 << 0, 1, 2, 1,  // 0
        1, 0, 1, 2,                           // 1
        2, 1, 0, 3,                           // 2
        1, 2, 3, 0;                           // 3
    REQUIRE(handi.get_distance_matrix() == correct_distance_matrix_2);

    aas::MatrixXu correct_path_matrix_2(4, 4);
    correct_path_matrix_2 << 0, 1, 1, 3,  // 0
        0, 1, 2, 0,                       // 1
        1, 1, 2, 1,                       // 2
        0, 0, 0, 3;                       // 3
    REQUIRE(handi.get_path_matrix() == correct_path_matrix_2);
  }
  GIVEN("acyclic path generation IV") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)},
         {Node(4), Node(0)}});

    aas::PathHandler handler(archi);

    aas::MatrixXu correct_distance_matrix(5, 5);
    correct_distance_matrix << 0, 1, 2, 2, 1,  // 0
        1, 0, 1, 2, 2,                         // 1
        2, 1, 0, 1, 2,                         // 2
        2, 2, 1, 0, 1,                         // 3
        1, 2, 2, 1, 0;                         // 4
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);

    aas::MatrixXu correct_path_matrix(5, 5);
    correct_path_matrix << 0, 1, 1, 4, 4,  // 0
        0, 1, 2, 2, 0,                     // 1
        1, 1, 2, 3, 3,                     // 2
        4, 2, 2, 3, 4,                     // 3
        0, 0, 3, 3, 4;                     // 4
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);

    aas::PathHandler handi = handler.construct_acyclic_handler();

    aas::MatrixXu correct_distance_matrix_2(5, 5);
    correct_distance_matrix_2 << 0, 1, 2, 2, 1,  // 0
        1, 0, 1, 3, 2,                           // 1
        2, 1, 0, 4, 3,                           // 2
        2, 3, 4, 0, 1,                           // 3
        1, 2, 3, 1, 0;                           // 4
    REQUIRE(handi.get_distance_matrix() == correct_distance_matrix_2);

    aas::MatrixXu correct_path_matrix_2(5, 5);
    correct_path_matrix_2 << 0, 1, 1, 4, 4,  // 0
        0, 1, 2, 0, 0,                       // 1
        1, 1, 2, 1, 1,                       // 2
        4, 4, 4, 3, 4,                       // 3
        0, 0, 0, 3, 4;                       // 4
    REQUIRE(handi.get_path_matrix() == correct_path_matrix_2);
  }
  GIVEN("acyclic path generation V") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)},
         {Node(4), Node(0)},
         {Node(5), Node(0)},
         {Node(5), Node(1)},
         {Node(5), Node(2)},
         {Node(5), Node(3)},
         {Node(5), Node(4)}});

    aas::PathHandler handler(archi);

    aas::MatrixXu correct_distance_matrix(6, 6);
    correct_distance_matrix << 0, 1, 2, 2, 1, 1,  // 0
        1, 0, 1, 2, 2, 1,                         // 1
        2, 1, 0, 1, 2, 1,                         // 2
        2, 2, 1, 0, 1, 1,                         // 3
        1, 2, 2, 1, 0, 1,                         // 4
        1, 1, 1, 1, 1, 0;                         // 5
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);

    aas::MatrixXu correct_path_matrix(6, 6);
    correct_path_matrix << 0, 1, 1, 4, 4, 5,  // 0
        0, 1, 2, 2, 0, 5,                     // 1
        1, 1, 2, 3, 3, 5,                     // 2
        4, 2, 2, 3, 4, 5,                     // 3
        0, 0, 3, 3, 4, 5,                     // 4
        0, 1, 2, 3, 4, 5;                     // 5
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);

    aas::PathHandler handi = handler.construct_acyclic_handler();

    aas::MatrixXu correct_distance_matrix_2(6, 6);
    correct_distance_matrix_2 << 0, 2, 2, 2, 2, 1,  // 0
        2, 0, 2, 2, 2, 1,                           // 1
        2, 2, 0, 2, 2, 1,                           // 2
        2, 2, 2, 0, 2, 1,                           // 3
        2, 2, 2, 2, 0, 1,                           // 4
        1, 1, 1, 1, 1, 0;                           // 5
    REQUIRE(handi.get_distance_matrix() == correct_distance_matrix_2);

    aas::MatrixXu correct_path_matrix_2(6, 6);
    correct_path_matrix_2 << 0, 5, 5, 5, 5, 5,  // 0
        5, 1, 5, 5, 5, 5,                       // 1
        5, 5, 2, 5, 5, 5,                       // 2
        5, 5, 5, 3, 5, 5,                       // 3
        5, 5, 5, 5, 4, 5,                       // 4
        0, 1, 2, 3, 4, 5;                       // 5
    REQUIRE(handi.get_path_matrix() == correct_path_matrix_2);
  }
  GIVEN("acyclic path generation VI") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(5), Node(6)},
         {Node(6), Node(7)},
         {Node(7), Node(0)},
         {Node(2), Node(0)},
         {Node(3), Node(0)},
         {Node(4), Node(0)},
         {Node(5), Node(0)},
         {Node(6), Node(0)}});

    aas::PathHandler handler(archi);

    aas::MatrixXu correct_distance_matrix(8, 8);
    correct_distance_matrix << 0, 1, 1, 1, 1, 1, 1, 1,  // 0
        1, 0, 1, 2, 2, 2, 2, 2,                         // 1
        1, 1, 0, 1, 2, 2, 2, 2,                         // 2
        1, 2, 1, 0, 1, 2, 2, 2,                         // 3
        1, 2, 2, 1, 0, 1, 2, 2,                         // 4
        1, 2, 2, 2, 1, 0, 1, 2,                         // 5
        1, 2, 2, 2, 2, 1, 0, 1,                         // 6
        1, 2, 2, 2, 2, 2, 1, 0;                         // 7
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);

    aas::MatrixXu correct_path_matrix(8, 8);
    correct_path_matrix << 0, 1, 2, 3, 4, 5, 6, 7,  // 0
        0, 1, 2, 0, 0, 0, 0, 0,                     // 1
        0, 1, 2, 3, 0, 0, 0, 0,                     // 2
        0, 0, 2, 3, 4, 0, 0, 0,                     // 3
        0, 0, 0, 3, 4, 5, 0, 0,                     // 4
        0, 0, 0, 0, 4, 5, 6, 0,                     // 5
        0, 0, 0, 0, 0, 5, 6, 7,                     // 6
        0, 0, 0, 0, 0, 0, 6, 7;                     // 7
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);

    aas::PathHandler handi = handler.construct_acyclic_handler();

    aas::MatrixXu correct_distance_matrix_2(8, 8);
    correct_distance_matrix_2 << 0, 1, 1, 1, 1, 1, 1, 1,  // 0
        1, 0, 2, 2, 2, 2, 2, 2,                           // 1
        1, 2, 0, 2, 2, 2, 2, 2,                           // 2
        1, 2, 2, 0, 2, 2, 2, 2,                           // 3
        1, 2, 2, 2, 0, 2, 2, 2,                           // 4
        1, 2, 2, 2, 2, 0, 2, 2,                           // 5
        1, 2, 2, 2, 2, 2, 0, 2,                           // 6
        1, 2, 2, 2, 2, 2, 2, 0;                           // 7
    REQUIRE(handi.get_distance_matrix() == correct_distance_matrix_2);

    aas::MatrixXu correct_path_matrix_2(8, 8);
    correct_path_matrix_2 << 0, 1, 2, 3, 4, 5, 6, 7,  // 0
        0, 1, 0, 0, 0, 0, 0, 0,                       // 1
        0, 0, 2, 0, 0, 0, 0, 0,                       // 2
        0, 0, 0, 3, 0, 0, 0, 0,                       // 3
        0, 0, 0, 0, 4, 0, 0, 0,                       // 4
        0, 0, 0, 0, 0, 5, 0, 0,                       // 5
        0, 0, 0, 0, 0, 0, 6, 0,                       // 6
        0, 0, 0, 0, 0, 0, 0, 7;                       // 7
    REQUIRE(handi.get_path_matrix() == correct_path_matrix_2);
  }
  GIVEN("acyclic path generation VII") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(3)},
         {Node(1), Node(2)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(3), Node(4)},
         {Node(4), Node(5)}});

    aas::PathHandler handler(archi);

    aas::MatrixXu correct_distance_matrix(6, 6);
    correct_distance_matrix << 0, 1, 2, 1, 2, 3,  // 0
        1, 0, 1, 2, 1, 2,                         // 1
        2, 1, 0, 3, 2, 1,                         // 2
        1, 2, 3, 0, 1, 2,                         // 3
        2, 1, 2, 1, 0, 1,                         // 4
        3, 2, 1, 2, 1, 0;                         // 5
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);

    aas::MatrixXu correct_path_matrix(6, 6);
    correct_path_matrix << 0, 1, 1, 3, 1, 1,  // 0
        0, 1, 2, 0, 4, 2,                     // 1
        1, 1, 2, 1, 1, 5,                     // 2
        0, 0, 0, 3, 4, 4,                     // 3
        1, 1, 1, 3, 4, 5,                     // 4
        2, 2, 2, 4, 4, 5;                     // 5
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);

    aas::PathHandler handi = handler.construct_acyclic_handler();

    aas::MatrixXu correct_distance_matrix_2(6, 6);
    correct_distance_matrix_2 << 0, 1, 2, 3, 2, 3,  // 0
        1, 0, 1, 2, 1, 2,                           // 1
        2, 1, 0, 3, 2, 3,                           // 2
        3, 2, 3, 0, 1, 2,                           // 3
        2, 1, 2, 1, 0, 1,                           // 4
        3, 2, 3, 2, 1, 0;                           // 5
    REQUIRE(handi.get_distance_matrix() == correct_distance_matrix_2);

    aas::MatrixXu correct_path_matrix_2(6, 6);
    correct_path_matrix_2 << 0, 1, 1, 1, 1, 1,  // 0
        0, 1, 2, 4, 4, 4,                       // 1
        1, 1, 2, 1, 1, 1,                       // 2
        4, 4, 4, 3, 4, 4,                       // 3
        1, 1, 1, 3, 4, 5,                       // 4
        4, 4, 4, 4, 4, 5;                       // 5
    REQUIRE(handi.get_path_matrix() == correct_path_matrix_2);
  }
  GIVEN("acyclic path generation VIII") {
    const Architecture archi(
        {{Node(0), Node(1)},
         {Node(0), Node(3)},
         {Node(1), Node(2)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(3), Node(4)},
         {Node(3), Node(6)},
         {Node(4), Node(5)},
         {Node(4), Node(7)},
         {Node(5), Node(8)},
         {Node(6), Node(7)},
         {Node(7), Node(8)}});

    aas::PathHandler handler(archi);

    aas::MatrixXu correct_distance_matrix(9, 9);
    correct_distance_matrix << 0, 1, 2, 1, 2, 3, 2, 3, 4,  // 0
        1, 0, 1, 2, 1, 2, 3, 2, 3,                         // 1
        2, 1, 0, 3, 2, 1, 4, 3, 2,                         // 2
        1, 2, 3, 0, 1, 2, 1, 2, 3,                         // 3
        2, 1, 2, 1, 0, 1, 2, 1, 2,                         // 4
        3, 2, 1, 2, 1, 0, 3, 2, 1,                         // 5
        2, 3, 4, 1, 2, 3, 0, 1, 2,                         // 6
        3, 2, 3, 2, 1, 2, 1, 0, 1,                         // 7
        4, 3, 2, 3, 2, 1, 2, 1, 0;                         // 8
    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);

    aas::MatrixXu correct_path_matrix(9, 9);
    correct_path_matrix << 0, 1, 1, 3, 1, 1, 3, 1, 1,  // 0
        0, 1, 2, 0, 4, 2, 0, 4, 2,                     // 1
        1, 1, 2, 1, 1, 5, 1, 1, 5,                     // 2
        0, 0, 0, 3, 4, 4, 6, 4, 4,                     // 3
        1, 1, 1, 3, 4, 5, 3, 7, 5,                     // 4
        2, 2, 2, 4, 4, 5, 4, 4, 8,                     // 5
        3, 3, 3, 3, 3, 3, 6, 7, 7,                     // 6
        4, 4, 4, 4, 4, 4, 6, 7, 8,                     // 7
        5, 5, 5, 5, 5, 5, 7, 7, 8;                     // 8
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);

    aas::PathHandler handi = handler.construct_acyclic_handler();

    aas::MatrixXu correct_distance_matrix_2(9, 9);
    correct_distance_matrix_2 << 0, 1, 2, 3, 2, 3, 4, 3, 4,  // 0
        1, 0, 1, 2, 1, 2, 3, 2, 3,                           // 1
        2, 1, 0, 3, 2, 3, 4, 3, 4,                           // 2
        3, 2, 3, 0, 1, 2, 1, 2, 3,                           // 3
        2, 1, 2, 1, 0, 1, 2, 1, 2,                           // 4
        3, 2, 3, 2, 1, 0, 3, 2, 1,                           // 5
        4, 3, 4, 1, 2, 3, 0, 3, 4,                           // 6
        3, 2, 3, 2, 1, 2, 3, 0, 3,                           // 7
        4, 3, 4, 3, 2, 1, 4, 3, 0;                           // 8
    REQUIRE(handi.get_distance_matrix() == correct_distance_matrix_2);

    aas::MatrixXu correct_path_matrix_2(9, 9);
    correct_path_matrix_2 << 0, 1, 1, 1, 1, 1, 1, 1, 1,  // 0
        0, 1, 2, 4, 4, 4, 4, 4, 4,                       // 1
        1, 1, 2, 1, 1, 1, 1, 1, 1,                       // 2
        4, 4, 4, 3, 4, 4, 6, 4, 4,                       // 3
        1, 1, 1, 3, 4, 5, 3, 7, 5,                       // 4
        4, 4, 4, 4, 4, 5, 4, 4, 8,                       // 5
        3, 3, 3, 3, 3, 3, 6, 3, 3,                       // 6
        4, 4, 4, 4, 4, 4, 4, 7, 4,                       // 7
        5, 5, 5, 5, 5, 5, 5, 5, 8;                       // 8
    REQUIRE(handi.get_path_matrix() == correct_path_matrix_2);
  }
  GIVEN("acyclic path generation IX") {
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

    aas::MatrixXu correct_distance_matrix(15, 15);
    correct_distance_matrix << 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
        3,                                            // 0
        1, 0, 2, 1, 1, 3, 3, 2, 2, 2, 2, 4, 4, 4, 4,  // 1
        1, 2, 0, 3, 3, 1, 1, 4, 4, 4, 4, 2, 2, 2, 2,  // 2
        2, 1, 3, 0, 2, 4, 4, 1, 1, 3, 3, 5, 5, 5, 5,  // 3
        2, 1, 3, 2, 0, 4, 4, 3, 3, 1, 1, 5, 5, 5, 5,  // 4
        2, 3, 1, 4, 4, 0, 2, 5, 5, 5, 5, 1, 1, 3, 3,  // 5
        2, 3, 1, 4, 4, 2, 0, 5, 5, 5, 5, 3, 3, 1, 1,  // 6
        3, 2, 4, 1, 3, 5, 5, 0, 2, 4, 4, 6, 6, 6, 6,  // 7
        3, 2, 4, 1, 3, 5, 5, 2, 0, 4, 4, 6, 6, 6, 6,  // 8
        3, 2, 4, 3, 1, 5, 5, 4, 4, 0, 2, 6, 6, 6, 6,  // 9
        3, 2, 4, 3, 1, 5, 5, 4, 4, 2, 0, 6, 6, 6, 6,  // 10
        3, 4, 2, 5, 5, 1, 3, 6, 6, 6, 6, 0, 2, 4, 4,  // 11
        3, 4, 2, 5, 5, 1, 3, 6, 6, 6, 6, 2, 0, 4, 4,  // 12
        3, 4, 2, 5, 5, 3, 1, 6, 6, 6, 6, 4, 4, 0, 2,  // 13
        3, 4, 2, 5, 5, 3, 1, 6, 6, 6, 6, 4, 4, 2, 0;  // 14

    REQUIRE(handler.get_distance_matrix() == correct_distance_matrix);

    aas::MatrixXu correct_path_matrix(15, 15);

    correct_path_matrix << 0, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2,
        2,                                              // 0
        0, 1, 0, 3, 4, 0, 0, 3, 3, 4, 4, 0, 0, 0, 0,    // 1
        0, 0, 2, 0, 0, 5, 6, 0, 0, 0, 0, 5, 5, 6, 6,    // 2
        1, 1, 1, 3, 1, 1, 1, 7, 8, 1, 1, 1, 1, 1, 1,    // 3
        1, 1, 1, 1, 4, 1, 1, 1, 1, 9, 10, 1, 1, 1, 1,   // 4
        2, 2, 2, 2, 2, 5, 2, 2, 2, 2, 2, 11, 12, 2, 2,  // 5
        2, 2, 2, 2, 2, 2, 6, 2, 2, 2, 2, 2, 2, 13, 14,  // 6
        3, 3, 3, 3, 3, 3, 3, 7, 3, 3, 3, 3, 3, 3, 3,    // 7
        3, 3, 3, 3, 3, 3, 3, 3, 8, 3, 3, 3, 3, 3, 3,    // 8
        4, 4, 4, 4, 4, 4, 4, 4, 4, 9, 4, 4, 4, 4, 4,    // 9
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4, 4,   // 10
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 11, 5, 5, 5,   // 11
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 12, 5, 5,   // 12
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 13, 6,   // 13
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 14;   // 14
    REQUIRE(handler.get_path_matrix() == correct_path_matrix);

    aas::PathHandler handi = handler.construct_acyclic_handler();

    aas::MatrixXu correct_distance_matrix_2(15, 15);
    correct_distance_matrix_2 << 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
        3,                                            // 0
        1, 0, 2, 1, 1, 3, 3, 2, 2, 2, 2, 4, 4, 4, 4,  // 1
        1, 2, 0, 3, 3, 1, 1, 4, 4, 4, 4, 2, 2, 2, 2,  // 2
        2, 1, 3, 0, 2, 4, 4, 1, 1, 3, 3, 5, 5, 5, 5,  // 3
        2, 1, 3, 2, 0, 4, 4, 3, 3, 1, 1, 5, 5, 5, 5,  // 4
        2, 3, 1, 4, 4, 0, 2, 5, 5, 5, 5, 1, 1, 3, 3,  // 5
        2, 3, 1, 4, 4, 2, 0, 5, 5, 5, 5, 3, 3, 1, 1,  // 6
        3, 2, 4, 1, 3, 5, 5, 0, 2, 4, 4, 6, 6, 6, 6,  // 7
        3, 2, 4, 1, 3, 5, 5, 2, 0, 4, 4, 6, 6, 6, 6,  // 8
        3, 2, 4, 3, 1, 5, 5, 4, 4, 0, 2, 6, 6, 6, 6,  // 9
        3, 2, 4, 3, 1, 5, 5, 4, 4, 2, 0, 6, 6, 6, 6,  // 10
        3, 4, 2, 5, 5, 1, 3, 6, 6, 6, 6, 0, 2, 4, 4,  // 11
        3, 4, 2, 5, 5, 1, 3, 6, 6, 6, 6, 2, 0, 4, 4,  // 12
        3, 4, 2, 5, 5, 3, 1, 6, 6, 6, 6, 4, 4, 0, 2,  // 13
        3, 4, 2, 5, 5, 3, 1, 6, 6, 6, 6, 4, 4, 2, 0;  // 14
    REQUIRE(handi.get_distance_matrix() == correct_distance_matrix_2);

    aas::MatrixXu correct_path_matrix_2(15, 15);
    correct_path_matrix_2 << 0, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2,
        2,                                              // 0
        0, 1, 0, 3, 4, 0, 0, 3, 3, 4, 4, 0, 0, 0, 0,    // 1
        0, 0, 2, 0, 0, 5, 6, 0, 0, 0, 0, 5, 5, 6, 6,    // 2
        1, 1, 1, 3, 1, 1, 1, 7, 8, 1, 1, 1, 1, 1, 1,    // 3
        1, 1, 1, 1, 4, 1, 1, 1, 1, 9, 10, 1, 1, 1, 1,   // 4
        2, 2, 2, 2, 2, 5, 2, 2, 2, 2, 2, 11, 12, 2, 2,  // 5
        2, 2, 2, 2, 2, 2, 6, 2, 2, 2, 2, 2, 2, 13, 14,  // 6
        3, 3, 3, 3, 3, 3, 3, 7, 3, 3, 3, 3, 3, 3, 3,    // 7
        3, 3, 3, 3, 3, 3, 3, 3, 8, 3, 3, 3, 3, 3, 3,    // 8
        4, 4, 4, 4, 4, 4, 4, 4, 4, 9, 4, 4, 4, 4, 4,    // 9
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 10, 4, 4, 4, 4,   // 10
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 11, 5, 5, 5,   // 11
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 12, 5, 5,   // 12
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 13, 6,   // 13
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 14;   // 14
    REQUIRE(handi.get_path_matrix() == correct_path_matrix_2);
  }
}
SCENARIO(
    "Check Hamiltonian path construction is correct - verylong",
    "[.verylong]") {
  GIVEN("acyclic path generation X - omp") {
    const Architecture archi(
        {{Node(0), Node(1)},      {Node(0), Node(2)},
         {Node(1), Node(2)},      {Node(1), Node(3)},
         {Node(2), Node(3)},      {Node(2), Node(4)},
         {Node(3), Node(4)},      {Node(3), Node(5)},
         {Node(4), Node(5)},      {Node(4), Node(6)},
         {Node(5), Node(6)},      {Node(5), Node(7)},
         {Node(6), Node(7)},      {Node(6), Node(8)},
         {Node(7), Node(8)},      {Node(7), Node(9)},
         {Node(8), Node(9)},      {Node(8), Node(10)},
         {Node(9), Node(10)},     {Node(9), Node(11)},
         {Node(10), Node(11)},    {Node(10), Node(12)},
         {Node(11), Node(12)},    {Node(11), Node(13)},
         {Node(12), Node(13)},    {Node(12), Node(14)},
         {Node(13), Node(14)},    {Node(13), Node(15)},
         {Node(14), Node(15)},    {Node(14), Node(16)},
         {Node(15), Node(16)},    {Node(15), Node(17)},
         {Node(16), Node(17)},    {Node(16), Node(18)},
         {Node(17), Node(18)},    {Node(17), Node(19)},
         {Node(18), Node(19)},    {Node(18), Node(20)},
         {Node(19), Node(20)},    {Node(19), Node(21)},
         {Node(20), Node(21)},    {Node(20), Node(22)},
         {Node(21), Node(22)},    {Node(21), Node(23)},
         {Node(22), Node(23)},    {Node(22), Node(24)},
         {Node(23), Node(24)},    {Node(23), Node(25)},
         {Node(24), Node(25)},    {Node(24), Node(26)},
         {Node(25), Node(26)},    {Node(25), Node(27)},
         {Node(26), Node(27)},    {Node(26), Node(28)},
         {Node(27), Node(28)},    {Node(27), Node(29)},
         {Node(28), Node(29)},    {Node(28), Node(30)},
         {Node(29), Node(30)},    {Node(29), Node(31)},
         {Node(30), Node(31)},    {Node(30), Node(32)},
         {Node(31), Node(32)},    {Node(31), Node(33)},
         {Node(32), Node(33)},    {Node(32), Node(34)},
         {Node(33), Node(34)},    {Node(33), Node(35)},
         {Node(34), Node(35)},    {Node(34), Node(36)},
         {Node(35), Node(36)},    {Node(35), Node(37)},
         {Node(36), Node(37)},    {Node(36), Node(38)},
         {Node(37), Node(38)},    {Node(37), Node(39)},
         {Node(38), Node(39)},    {Node(38), Node(40)},
         {Node(39), Node(40)},    {Node(39), Node(41)},
         {Node(40), Node(41)},    {Node(40), Node(42)},
         {Node(41), Node(42)},    {Node(41), Node(43)},
         {Node(42), Node(43)},    {Node(42), Node(44)},
         {Node(43), Node(44)},    {Node(43), Node(45)},
         {Node(44), Node(45)},    {Node(44), Node(46)},
         {Node(45), Node(46)},    {Node(45), Node(47)},
         {Node(46), Node(47)},    {Node(46), Node(48)},
         {Node(47), Node(48)},    {Node(47), Node(49)},
         {Node(48), Node(49)},    {Node(48), Node(50)},
         {Node(49), Node(50)},    {Node(49), Node(51)},
         {Node(50), Node(51)},    {Node(50), Node(52)},
         {Node(51), Node(52)},    {Node(51), Node(53)},
         {Node(52), Node(53)},    {Node(52), Node(54)},
         {Node(53), Node(54)},    {Node(53), Node(55)},
         {Node(54), Node(55)},    {Node(54), Node(56)},
         {Node(55), Node(56)},    {Node(55), Node(57)},
         {Node(56), Node(57)},    {Node(56), Node(58)},
         {Node(57), Node(58)},    {Node(57), Node(59)},
         {Node(58), Node(59)},    {Node(58), Node(60)},
         {Node(59), Node(60)},    {Node(59), Node(61)},
         {Node(60), Node(61)},    {Node(60), Node(62)},
         {Node(61), Node(62)},    {Node(61), Node(63)},
         {Node(62), Node(63)},    {Node(62), Node(64)},
         {Node(63), Node(64)},    {Node(63), Node(65)},
         {Node(64), Node(65)},    {Node(64), Node(66)},
         {Node(65), Node(66)},    {Node(65), Node(67)},
         {Node(66), Node(67)},    {Node(66), Node(68)},
         {Node(67), Node(68)},    {Node(67), Node(69)},
         {Node(68), Node(69)},    {Node(68), Node(70)},
         {Node(69), Node(70)},    {Node(69), Node(71)},
         {Node(70), Node(71)},    {Node(70), Node(72)},
         {Node(71), Node(72)},    {Node(71), Node(73)},
         {Node(72), Node(73)},    {Node(72), Node(74)},
         {Node(73), Node(74)},    {Node(73), Node(75)},
         {Node(74), Node(75)},    {Node(74), Node(76)},
         {Node(75), Node(76)},    {Node(75), Node(77)},
         {Node(76), Node(77)},    {Node(76), Node(78)},
         {Node(77), Node(78)},    {Node(77), Node(79)},
         {Node(78), Node(79)},    {Node(78), Node(80)},
         {Node(79), Node(80)},    {Node(79), Node(81)},
         {Node(80), Node(81)},    {Node(80), Node(82)},
         {Node(81), Node(82)},    {Node(81), Node(83)},
         {Node(82), Node(83)},    {Node(82), Node(84)},
         {Node(83), Node(84)},    {Node(83), Node(85)},
         {Node(84), Node(85)},    {Node(84), Node(86)},
         {Node(85), Node(86)},    {Node(85), Node(87)},
         {Node(86), Node(87)},    {Node(86), Node(88)},
         {Node(87), Node(88)},    {Node(87), Node(89)},
         {Node(88), Node(89)},    {Node(88), Node(90)},
         {Node(89), Node(90)},    {Node(89), Node(91)},
         {Node(90), Node(91)},    {Node(90), Node(92)},
         {Node(91), Node(92)},    {Node(91), Node(93)},
         {Node(92), Node(93)},    {Node(92), Node(94)},
         {Node(93), Node(94)},    {Node(93), Node(95)},
         {Node(94), Node(95)},    {Node(94), Node(96)},
         {Node(95), Node(96)},    {Node(95), Node(97)},
         {Node(96), Node(97)},    {Node(96), Node(98)},
         {Node(97), Node(98)},    {Node(97), Node(99)},
         {Node(98), Node(99)},    {Node(98), Node(100)},
         {Node(99), Node(100)},   {Node(99), Node(101)},
         {Node(100), Node(101)},  {Node(100), Node(102)},
         {Node(101), Node(102)},  {Node(101), Node(103)},
         {Node(102), Node(103)},  {Node(102), Node(104)},
         {Node(103), Node(104)},  {Node(103), Node(105)},
         {Node(104), Node(105)},  {Node(104), Node(106)},
         {Node(105), Node(106)},  {Node(105), Node(107)},
         {Node(106), Node(107)},  {Node(106), Node(108)},
         {Node(107), Node(108)},  {Node(107), Node(109)},
         {Node(108), Node(109)},  {Node(108), Node(110)},
         {Node(109), Node(110)},  {Node(109), Node(111)},
         {Node(110), Node(111)},  {Node(110), Node(112)},
         {Node(111), Node(112)},  {Node(111), Node(113)},
         {Node(112), Node(113)},  {Node(112), Node(114)},
         {Node(113), Node(114)},  {Node(113), Node(115)},
         {Node(114), Node(115)},  {Node(114), Node(116)},
         {Node(115), Node(116)},  {Node(115), Node(117)},
         {Node(116), Node(117)},  {Node(116), Node(118)},
         {Node(117), Node(118)},  {Node(117), Node(119)},
         {Node(118), Node(119)},  {Node(118), Node(120)},
         {Node(119), Node(120)},  {Node(119), Node(121)},
         {Node(120), Node(121)},  {Node(120), Node(122)},
         {Node(121), Node(122)},  {Node(121), Node(123)},
         {Node(122), Node(123)},  {Node(122), Node(124)},
         {Node(123), Node(124)},  {Node(123), Node(125)},
         {Node(124), Node(125)},  {Node(124), Node(126)},
         {Node(125), Node(126)},  {Node(125), Node(127)},
         {Node(126), Node(127)},  {Node(126), Node(128)},
         {Node(127), Node(128)},  {Node(127), Node(129)},
         {Node(128), Node(129)},  {Node(128), Node(130)},
         {Node(129), Node(130)},  {Node(129), Node(131)},
         {Node(130), Node(131)},  {Node(130), Node(132)},
         {Node(131), Node(132)},  {Node(131), Node(133)},
         {Node(132), Node(133)},  {Node(132), Node(134)},
         {Node(133), Node(134)},  {Node(133), Node(135)},
         {Node(134), Node(135)},  {Node(134), Node(136)},
         {Node(135), Node(136)},  {Node(135), Node(137)},
         {Node(136), Node(137)},  {Node(136), Node(138)},
         {Node(137), Node(138)},  {Node(137), Node(139)},
         {Node(138), Node(139)},  {Node(138), Node(140)},
         {Node(139), Node(140)},  {Node(139), Node(141)},
         {Node(140), Node(141)},  {Node(140), Node(142)},
         {Node(141), Node(142)},  {Node(141), Node(143)},
         {Node(142), Node(143)},  {Node(142), Node(144)},
         {Node(143), Node(144)},  {Node(143), Node(145)},
         {Node(144), Node(145)},  {Node(144), Node(146)},
         {Node(145), Node(146)},  {Node(145), Node(147)},
         {Node(146), Node(147)},  {Node(146), Node(148)},
         {Node(147), Node(148)},  {Node(147), Node(149)},
         {Node(148), Node(149)},  {Node(148), Node(150)},
         {Node(149), Node(150)},  {Node(149), Node(151)},
         {Node(150), Node(151)},  {Node(150), Node(152)},
         {Node(151), Node(152)},  {Node(151), Node(153)},
         {Node(152), Node(153)},  {Node(152), Node(154)},
         {Node(153), Node(154)},  {Node(153), Node(155)},
         {Node(154), Node(155)},  {Node(154), Node(156)},
         {Node(155), Node(156)},  {Node(155), Node(157)},
         {Node(156), Node(157)},  {Node(156), Node(158)},
         {Node(157), Node(158)},  {Node(157), Node(159)},
         {Node(158), Node(159)},  {Node(158), Node(160)},
         {Node(159), Node(160)},  {Node(159), Node(161)},
         {Node(160), Node(161)},  {Node(160), Node(162)},
         {Node(161), Node(162)},  {Node(161), Node(163)},
         {Node(162), Node(163)},  {Node(162), Node(164)},
         {Node(163), Node(164)},  {Node(163), Node(165)},
         {Node(164), Node(165)},  {Node(164), Node(166)},
         {Node(165), Node(166)},  {Node(165), Node(167)},
         {Node(166), Node(167)},  {Node(166), Node(168)},
         {Node(167), Node(168)},  {Node(167), Node(169)},
         {Node(168), Node(169)},  {Node(168), Node(170)},
         {Node(169), Node(170)},  {Node(169), Node(171)},
         {Node(170), Node(171)},  {Node(170), Node(172)},
         {Node(171), Node(172)},  {Node(171), Node(173)},
         {Node(172), Node(173)},  {Node(172), Node(174)},
         {Node(173), Node(174)},  {Node(173), Node(175)},
         {Node(174), Node(175)},  {Node(174), Node(176)},
         {Node(175), Node(176)},  {Node(175), Node(177)},
         {Node(176), Node(177)},  {Node(176), Node(178)},
         {Node(177), Node(178)},  {Node(177), Node(179)},
         {Node(178), Node(179)},  {Node(178), Node(180)},
         {Node(179), Node(180)},  {Node(179), Node(181)},
         {Node(180), Node(181)},  {Node(180), Node(182)},
         {Node(181), Node(182)},  {Node(181), Node(183)},
         {Node(182), Node(183)},  {Node(182), Node(184)},
         {Node(183), Node(184)},  {Node(183), Node(185)},
         {Node(184), Node(185)},  {Node(184), Node(186)},
         {Node(185), Node(186)},  {Node(185), Node(187)},
         {Node(186), Node(187)},  {Node(186), Node(188)},
         {Node(187), Node(188)},  {Node(187), Node(189)},
         {Node(188), Node(189)},  {Node(188), Node(190)},
         {Node(189), Node(190)},  {Node(189), Node(191)},
         {Node(190), Node(191)},  {Node(190), Node(192)},
         {Node(191), Node(192)},  {Node(191), Node(193)},
         {Node(192), Node(193)},  {Node(192), Node(194)},
         {Node(193), Node(194)},  {Node(193), Node(195)},
         {Node(194), Node(195)},  {Node(194), Node(196)},
         {Node(195), Node(196)},  {Node(195), Node(197)},
         {Node(196), Node(197)},  {Node(196), Node(198)},
         {Node(197), Node(198)},  {Node(197), Node(199)},
         {Node(198), Node(199)},  {Node(198), Node(200)},
         {Node(199), Node(200)},  {Node(199), Node(201)},
         {Node(200), Node(201)},  {Node(200), Node(202)},
         {Node(201), Node(202)},  {Node(201), Node(203)},
         {Node(202), Node(203)},  {Node(202), Node(204)},
         {Node(203), Node(204)},  {Node(203), Node(205)},
         {Node(204), Node(205)},  {Node(204), Node(206)},
         {Node(205), Node(206)},  {Node(205), Node(207)},
         {Node(206), Node(207)},  {Node(206), Node(208)},
         {Node(207), Node(208)},  {Node(207), Node(209)},
         {Node(208), Node(209)},  {Node(208), Node(210)},
         {Node(209), Node(210)},  {Node(209), Node(211)},
         {Node(210), Node(211)},  {Node(210), Node(212)},
         {Node(211), Node(212)},  {Node(211), Node(213)},
         {Node(212), Node(213)},  {Node(212), Node(214)},
         {Node(213), Node(214)},  {Node(213), Node(215)},
         {Node(214), Node(215)},  {Node(214), Node(216)},
         {Node(215), Node(216)},  {Node(215), Node(217)},
         {Node(216), Node(217)},  {Node(216), Node(218)},
         {Node(217), Node(218)},  {Node(217), Node(219)},
         {Node(218), Node(219)},  {Node(218), Node(220)},
         {Node(219), Node(220)},  {Node(219), Node(221)},
         {Node(220), Node(221)},  {Node(220), Node(222)},
         {Node(221), Node(222)},  {Node(221), Node(223)},
         {Node(222), Node(223)},  {Node(222), Node(224)},
         {Node(223), Node(224)},  {Node(223), Node(225)},
         {Node(224), Node(225)},  {Node(224), Node(226)},
         {Node(225), Node(226)},  {Node(225), Node(227)},
         {Node(226), Node(227)},  {Node(226), Node(228)},
         {Node(227), Node(228)},  {Node(227), Node(229)},
         {Node(228), Node(229)},  {Node(228), Node(230)},
         {Node(229), Node(230)},  {Node(229), Node(231)},
         {Node(230), Node(231)},  {Node(230), Node(232)},
         {Node(231), Node(232)},  {Node(231), Node(233)},
         {Node(232), Node(233)},  {Node(232), Node(234)},
         {Node(233), Node(234)},  {Node(233), Node(235)},
         {Node(234), Node(235)},  {Node(234), Node(236)},
         {Node(235), Node(236)},  {Node(235), Node(237)},
         {Node(236), Node(237)},  {Node(236), Node(238)},
         {Node(237), Node(238)},  {Node(237), Node(239)},
         {Node(238), Node(239)},  {Node(238), Node(240)},
         {Node(239), Node(240)},  {Node(239), Node(241)},
         {Node(240), Node(241)},  {Node(240), Node(242)},
         {Node(241), Node(242)},  {Node(241), Node(243)},
         {Node(242), Node(243)},  {Node(242), Node(244)},
         {Node(243), Node(244)},  {Node(243), Node(245)},
         {Node(244), Node(245)},  {Node(244), Node(246)},
         {Node(245), Node(246)},  {Node(245), Node(247)},
         {Node(246), Node(247)},  {Node(246), Node(248)},
         {Node(247), Node(248)},  {Node(247), Node(249)},
         {Node(248), Node(249)},  {Node(248), Node(250)},
         {Node(249), Node(250)},  {Node(249), Node(251)},
         {Node(250), Node(251)},  {Node(250), Node(252)},
         {Node(251), Node(252)},  {Node(251), Node(253)},
         {Node(252), Node(253)},  {Node(252), Node(254)},
         {Node(253), Node(254)},  {Node(253), Node(255)},
         {Node(254), Node(255)},  {Node(254), Node(256)},
         {Node(255), Node(256)},  {Node(255), Node(257)},
         {Node(256), Node(257)},  {Node(256), Node(258)},
         {Node(257), Node(258)},  {Node(257), Node(259)},
         {Node(258), Node(259)},  {Node(258), Node(260)},
         {Node(259), Node(260)},  {Node(259), Node(261)},
         {Node(260), Node(261)},  {Node(260), Node(262)},
         {Node(261), Node(262)},  {Node(261), Node(263)},
         {Node(262), Node(263)},  {Node(262), Node(264)},
         {Node(263), Node(264)},  {Node(263), Node(265)},
         {Node(264), Node(265)},  {Node(264), Node(266)},
         {Node(265), Node(266)},  {Node(265), Node(267)},
         {Node(266), Node(267)},  {Node(266), Node(268)},
         {Node(267), Node(268)},  {Node(267), Node(269)},
         {Node(268), Node(269)},  {Node(268), Node(270)},
         {Node(269), Node(270)},  {Node(269), Node(271)},
         {Node(270), Node(271)},  {Node(270), Node(272)},
         {Node(271), Node(272)},  {Node(271), Node(273)},
         {Node(272), Node(273)},  {Node(272), Node(274)},
         {Node(273), Node(274)},  {Node(273), Node(275)},
         {Node(274), Node(275)},  {Node(274), Node(276)},
         {Node(275), Node(276)},  {Node(275), Node(277)},
         {Node(276), Node(277)},  {Node(276), Node(278)},
         {Node(277), Node(278)},  {Node(277), Node(279)},
         {Node(278), Node(279)},  {Node(278), Node(280)},
         {Node(279), Node(280)},  {Node(279), Node(281)},
         {Node(280), Node(281)},  {Node(280), Node(282)},
         {Node(281), Node(282)},  {Node(281), Node(283)},
         {Node(282), Node(283)},  {Node(282), Node(284)},
         {Node(283), Node(284)},  {Node(283), Node(285)},
         {Node(284), Node(285)},  {Node(284), Node(286)},
         {Node(285), Node(286)},  {Node(285), Node(287)},
         {Node(286), Node(287)},  {Node(286), Node(288)},
         {Node(287), Node(288)},  {Node(287), Node(289)},
         {Node(288), Node(289)},  {Node(288), Node(290)},
         {Node(289), Node(290)},  {Node(289), Node(291)},
         {Node(290), Node(291)},  {Node(290), Node(292)},
         {Node(291), Node(292)},  {Node(291), Node(293)},
         {Node(292), Node(293)},  {Node(292), Node(294)},
         {Node(293), Node(294)},  {Node(293), Node(295)},
         {Node(294), Node(295)},  {Node(294), Node(296)},
         {Node(295), Node(296)},  {Node(295), Node(297)},
         {Node(296), Node(297)},  {Node(296), Node(298)},
         {Node(297), Node(298)},  {Node(297), Node(299)},
         {Node(298), Node(299)},  {Node(298), Node(300)},
         {Node(299), Node(300)},  {Node(299), Node(301)},
         {Node(300), Node(301)},  {Node(300), Node(302)},
         {Node(301), Node(302)},  {Node(301), Node(303)},
         {Node(302), Node(303)},  {Node(302), Node(304)},
         {Node(303), Node(304)},  {Node(303), Node(305)},
         {Node(304), Node(305)},  {Node(304), Node(306)},
         {Node(305), Node(306)},  {Node(305), Node(307)},
         {Node(306), Node(307)},  {Node(306), Node(308)},
         {Node(307), Node(308)},  {Node(307), Node(309)},
         {Node(308), Node(309)},  {Node(308), Node(310)},
         {Node(309), Node(310)},  {Node(309), Node(311)},
         {Node(310), Node(311)},  {Node(310), Node(312)},
         {Node(311), Node(312)},  {Node(311), Node(313)},
         {Node(312), Node(313)},  {Node(312), Node(314)},
         {Node(313), Node(314)},  {Node(313), Node(315)},
         {Node(314), Node(315)},  {Node(314), Node(316)},
         {Node(315), Node(316)},  {Node(315), Node(317)},
         {Node(316), Node(317)},  {Node(316), Node(318)},
         {Node(317), Node(318)},  {Node(317), Node(319)},
         {Node(318), Node(319)},  {Node(318), Node(320)},
         {Node(319), Node(320)},  {Node(319), Node(321)},
         {Node(320), Node(321)},  {Node(320), Node(322)},
         {Node(321), Node(322)},  {Node(321), Node(323)},
         {Node(322), Node(323)},  {Node(322), Node(324)},
         {Node(323), Node(324)},  {Node(323), Node(325)},
         {Node(324), Node(325)},  {Node(324), Node(326)},
         {Node(325), Node(326)},  {Node(325), Node(327)},
         {Node(326), Node(327)},  {Node(326), Node(328)},
         {Node(327), Node(328)},  {Node(327), Node(329)},
         {Node(328), Node(329)},  {Node(328), Node(330)},
         {Node(329), Node(330)},  {Node(329), Node(331)},
         {Node(330), Node(331)},  {Node(330), Node(332)},
         {Node(331), Node(332)},  {Node(331), Node(333)},
         {Node(332), Node(333)},  {Node(332), Node(334)},
         {Node(333), Node(334)},  {Node(333), Node(335)},
         {Node(334), Node(335)},  {Node(334), Node(336)},
         {Node(335), Node(336)},  {Node(335), Node(337)},
         {Node(336), Node(337)},  {Node(336), Node(338)},
         {Node(337), Node(338)},  {Node(337), Node(339)},
         {Node(338), Node(339)},  {Node(338), Node(340)},
         {Node(339), Node(340)},  {Node(339), Node(341)},
         {Node(340), Node(341)},  {Node(340), Node(342)},
         {Node(341), Node(342)},  {Node(341), Node(343)},
         {Node(342), Node(343)},  {Node(342), Node(344)},
         {Node(343), Node(344)},  {Node(343), Node(345)},
         {Node(344), Node(345)},  {Node(344), Node(346)},
         {Node(345), Node(346)},  {Node(345), Node(347)},
         {Node(346), Node(347)},  {Node(346), Node(348)},
         {Node(347), Node(348)},  {Node(347), Node(349)},
         {Node(348), Node(349)},  {Node(348), Node(350)},
         {Node(349), Node(350)},  {Node(349), Node(351)},
         {Node(350), Node(351)},  {Node(350), Node(352)},
         {Node(351), Node(352)},  {Node(351), Node(353)},
         {Node(352), Node(353)},  {Node(352), Node(354)},
         {Node(353), Node(354)},  {Node(353), Node(355)},
         {Node(354), Node(355)},  {Node(354), Node(356)},
         {Node(355), Node(356)},  {Node(355), Node(357)},
         {Node(356), Node(357)},  {Node(356), Node(358)},
         {Node(357), Node(358)},  {Node(357), Node(359)},
         {Node(358), Node(359)},  {Node(358), Node(360)},
         {Node(359), Node(360)},  {Node(359), Node(361)},
         {Node(360), Node(361)},  {Node(360), Node(362)},
         {Node(361), Node(362)},  {Node(361), Node(363)},
         {Node(362), Node(363)},  {Node(362), Node(364)},
         {Node(363), Node(364)},  {Node(363), Node(365)},
         {Node(364), Node(365)},  {Node(364), Node(366)},
         {Node(365), Node(366)},  {Node(365), Node(367)},
         {Node(366), Node(367)},  {Node(366), Node(368)},
         {Node(367), Node(368)},  {Node(367), Node(369)},
         {Node(368), Node(369)},  {Node(368), Node(370)},
         {Node(369), Node(370)},  {Node(369), Node(371)},
         {Node(370), Node(371)},  {Node(370), Node(372)},
         {Node(371), Node(372)},  {Node(371), Node(373)},
         {Node(372), Node(373)},  {Node(372), Node(374)},
         {Node(373), Node(374)},  {Node(373), Node(375)},
         {Node(374), Node(375)},  {Node(374), Node(376)},
         {Node(375), Node(376)},  {Node(375), Node(377)},
         {Node(376), Node(377)},  {Node(376), Node(378)},
         {Node(377), Node(378)},  {Node(377), Node(379)},
         {Node(378), Node(379)},  {Node(378), Node(380)},
         {Node(379), Node(380)},  {Node(379), Node(381)},
         {Node(380), Node(381)},  {Node(380), Node(382)},
         {Node(381), Node(382)},  {Node(381), Node(383)},
         {Node(382), Node(383)},  {Node(382), Node(384)},
         {Node(383), Node(384)},  {Node(383), Node(385)},
         {Node(384), Node(385)},  {Node(384), Node(386)},
         {Node(385), Node(386)},  {Node(385), Node(387)},
         {Node(386), Node(387)},  {Node(386), Node(388)},
         {Node(387), Node(388)},  {Node(387), Node(389)},
         {Node(388), Node(389)},  {Node(388), Node(390)},
         {Node(389), Node(390)},  {Node(389), Node(391)},
         {Node(390), Node(391)},  {Node(390), Node(392)},
         {Node(391), Node(392)},  {Node(391), Node(393)},
         {Node(392), Node(393)},  {Node(392), Node(394)},
         {Node(393), Node(394)},  {Node(393), Node(395)},
         {Node(394), Node(395)},  {Node(394), Node(396)},
         {Node(395), Node(396)},  {Node(395), Node(397)},
         {Node(396), Node(397)},  {Node(396), Node(398)},
         {Node(397), Node(398)},  {Node(397), Node(399)},
         {Node(398), Node(399)},  {Node(398), Node(400)},
         {Node(399), Node(400)},  {Node(399), Node(401)},
         {Node(400), Node(401)},  {Node(400), Node(402)},
         {Node(401), Node(402)},  {Node(401), Node(403)},
         {Node(402), Node(403)},  {Node(402), Node(404)},
         {Node(403), Node(404)},  {Node(403), Node(405)},
         {Node(404), Node(405)},  {Node(404), Node(406)},
         {Node(405), Node(406)},  {Node(405), Node(407)},
         {Node(406), Node(407)},  {Node(406), Node(408)},
         {Node(407), Node(408)},  {Node(407), Node(409)},
         {Node(408), Node(409)},  {Node(408), Node(410)},
         {Node(409), Node(410)},  {Node(409), Node(411)},
         {Node(410), Node(411)},  {Node(410), Node(412)},
         {Node(411), Node(412)},  {Node(411), Node(413)},
         {Node(412), Node(413)},  {Node(412), Node(414)},
         {Node(413), Node(414)},  {Node(413), Node(415)},
         {Node(414), Node(415)},  {Node(414), Node(416)},
         {Node(415), Node(416)},  {Node(415), Node(417)},
         {Node(416), Node(417)},  {Node(416), Node(418)},
         {Node(417), Node(418)},  {Node(417), Node(419)},
         {Node(418), Node(419)},  {Node(418), Node(420)},
         {Node(419), Node(420)},  {Node(419), Node(421)},
         {Node(420), Node(421)},  {Node(420), Node(422)},
         {Node(421), Node(422)},  {Node(421), Node(423)},
         {Node(422), Node(423)},  {Node(422), Node(424)},
         {Node(423), Node(424)},  {Node(423), Node(425)},
         {Node(424), Node(425)},  {Node(424), Node(426)},
         {Node(425), Node(426)},  {Node(425), Node(427)},
         {Node(426), Node(427)},  {Node(426), Node(428)},
         {Node(427), Node(428)},  {Node(427), Node(429)},
         {Node(428), Node(429)},  {Node(428), Node(430)},
         {Node(429), Node(430)},  {Node(429), Node(431)},
         {Node(430), Node(431)},  {Node(430), Node(432)},
         {Node(431), Node(432)},  {Node(431), Node(433)},
         {Node(432), Node(433)},  {Node(432), Node(434)},
         {Node(433), Node(434)},  {Node(433), Node(435)},
         {Node(434), Node(435)},  {Node(434), Node(436)},
         {Node(435), Node(436)},  {Node(435), Node(437)},
         {Node(436), Node(437)},  {Node(436), Node(438)},
         {Node(437), Node(438)},  {Node(437), Node(439)},
         {Node(438), Node(439)},  {Node(438), Node(440)},
         {Node(439), Node(440)},  {Node(439), Node(441)},
         {Node(440), Node(441)},  {Node(440), Node(442)},
         {Node(441), Node(442)},  {Node(441), Node(443)},
         {Node(442), Node(443)},  {Node(442), Node(444)},
         {Node(443), Node(444)},  {Node(443), Node(445)},
         {Node(444), Node(445)},  {Node(444), Node(446)},
         {Node(445), Node(446)},  {Node(445), Node(447)},
         {Node(446), Node(447)},  {Node(446), Node(448)},
         {Node(447), Node(448)},  {Node(447), Node(449)},
         {Node(448), Node(449)},  {Node(448), Node(450)},
         {Node(449), Node(450)},  {Node(449), Node(451)},
         {Node(450), Node(451)},  {Node(450), Node(452)},
         {Node(451), Node(452)},  {Node(451), Node(453)},
         {Node(452), Node(453)},  {Node(452), Node(454)},
         {Node(453), Node(454)},  {Node(453), Node(455)},
         {Node(454), Node(455)},  {Node(454), Node(456)},
         {Node(455), Node(456)},  {Node(455), Node(457)},
         {Node(456), Node(457)},  {Node(456), Node(458)},
         {Node(457), Node(458)},  {Node(457), Node(459)},
         {Node(458), Node(459)},  {Node(458), Node(460)},
         {Node(459), Node(460)},  {Node(459), Node(461)},
         {Node(460), Node(461)},  {Node(460), Node(462)},
         {Node(461), Node(462)},  {Node(461), Node(463)},
         {Node(462), Node(463)},  {Node(462), Node(464)},
         {Node(463), Node(464)},  {Node(463), Node(465)},
         {Node(464), Node(465)},  {Node(464), Node(466)},
         {Node(465), Node(466)},  {Node(465), Node(467)},
         {Node(466), Node(467)},  {Node(466), Node(468)},
         {Node(467), Node(468)},  {Node(467), Node(469)},
         {Node(468), Node(469)},  {Node(468), Node(470)},
         {Node(469), Node(470)},  {Node(469), Node(471)},
         {Node(470), Node(471)},  {Node(470), Node(472)},
         {Node(471), Node(472)},  {Node(471), Node(473)},
         {Node(472), Node(473)},  {Node(472), Node(474)},
         {Node(473), Node(474)},  {Node(473), Node(475)},
         {Node(474), Node(475)},  {Node(474), Node(476)},
         {Node(475), Node(476)},  {Node(475), Node(477)},
         {Node(476), Node(477)},  {Node(476), Node(478)},
         {Node(477), Node(478)},  {Node(477), Node(479)},
         {Node(478), Node(479)},  {Node(478), Node(480)},
         {Node(479), Node(480)},  {Node(479), Node(481)},
         {Node(480), Node(481)},  {Node(480), Node(482)},
         {Node(481), Node(482)},  {Node(481), Node(483)},
         {Node(482), Node(483)},  {Node(482), Node(484)},
         {Node(483), Node(484)},  {Node(483), Node(485)},
         {Node(484), Node(485)},  {Node(484), Node(486)},
         {Node(485), Node(486)},  {Node(485), Node(487)},
         {Node(486), Node(487)},  {Node(486), Node(488)},
         {Node(487), Node(488)},  {Node(487), Node(489)},
         {Node(488), Node(489)},  {Node(488), Node(490)},
         {Node(489), Node(490)},  {Node(489), Node(491)},
         {Node(490), Node(491)},  {Node(490), Node(492)},
         {Node(491), Node(492)},  {Node(491), Node(493)},
         {Node(492), Node(493)},  {Node(492), Node(494)},
         {Node(493), Node(494)},  {Node(493), Node(495)},
         {Node(494), Node(495)},  {Node(494), Node(496)},
         {Node(495), Node(496)},  {Node(495), Node(497)},
         {Node(496), Node(497)},  {Node(496), Node(498)},
         {Node(497), Node(498)},  {Node(497), Node(499)},
         {Node(498), Node(499)},  {Node(498), Node(500)},
         {Node(499), Node(500)},  {Node(499), Node(501)},
         {Node(500), Node(501)},  {Node(500), Node(502)},
         {Node(501), Node(502)},  {Node(501), Node(503)},
         {Node(502), Node(503)},  {Node(502), Node(504)},
         {Node(503), Node(504)},  {Node(503), Node(505)},
         {Node(504), Node(505)},  {Node(504), Node(506)},
         {Node(505), Node(506)},  {Node(505), Node(507)},
         {Node(506), Node(507)},  {Node(506), Node(508)},
         {Node(507), Node(508)},  {Node(507), Node(509)},
         {Node(508), Node(509)},  {Node(508), Node(510)},
         {Node(509), Node(510)},  {Node(509), Node(511)},
         {Node(510), Node(511)},  {Node(510), Node(512)},
         {Node(511), Node(512)},  {Node(511), Node(513)},
         {Node(512), Node(513)},  {Node(512), Node(514)},
         {Node(513), Node(514)},  {Node(513), Node(515)},
         {Node(514), Node(515)},  {Node(514), Node(516)},
         {Node(515), Node(516)},  {Node(515), Node(517)},
         {Node(516), Node(517)},  {Node(516), Node(518)},
         {Node(517), Node(518)},  {Node(517), Node(519)},
         {Node(518), Node(519)},  {Node(518), Node(520)},
         {Node(519), Node(520)},  {Node(519), Node(521)},
         {Node(520), Node(521)},  {Node(520), Node(522)},
         {Node(521), Node(522)},  {Node(521), Node(523)},
         {Node(522), Node(523)},  {Node(522), Node(524)},
         {Node(523), Node(524)},  {Node(523), Node(525)},
         {Node(524), Node(525)},  {Node(524), Node(526)},
         {Node(525), Node(526)},  {Node(525), Node(527)},
         {Node(526), Node(527)},  {Node(526), Node(528)},
         {Node(527), Node(528)},  {Node(527), Node(529)},
         {Node(528), Node(529)},  {Node(528), Node(530)},
         {Node(529), Node(530)},  {Node(529), Node(531)},
         {Node(530), Node(531)},  {Node(530), Node(532)},
         {Node(531), Node(532)},  {Node(531), Node(533)},
         {Node(532), Node(533)},  {Node(532), Node(534)},
         {Node(533), Node(534)},  {Node(533), Node(535)},
         {Node(534), Node(535)},  {Node(534), Node(536)},
         {Node(535), Node(536)},  {Node(535), Node(537)},
         {Node(536), Node(537)},  {Node(536), Node(538)},
         {Node(537), Node(538)},  {Node(537), Node(539)},
         {Node(538), Node(539)},  {Node(538), Node(540)},
         {Node(539), Node(540)},  {Node(539), Node(541)},
         {Node(540), Node(541)},  {Node(540), Node(542)},
         {Node(541), Node(542)},  {Node(541), Node(543)},
         {Node(542), Node(543)},  {Node(542), Node(544)},
         {Node(543), Node(544)},  {Node(543), Node(545)},
         {Node(544), Node(545)},  {Node(544), Node(546)},
         {Node(545), Node(546)},  {Node(545), Node(547)},
         {Node(546), Node(547)},  {Node(546), Node(548)},
         {Node(547), Node(548)},  {Node(547), Node(549)},
         {Node(548), Node(549)},  {Node(548), Node(550)},
         {Node(549), Node(550)},  {Node(549), Node(551)},
         {Node(550), Node(551)},  {Node(550), Node(552)},
         {Node(551), Node(552)},  {Node(551), Node(553)},
         {Node(552), Node(553)},  {Node(552), Node(554)},
         {Node(553), Node(554)},  {Node(553), Node(555)},
         {Node(554), Node(555)},  {Node(554), Node(556)},
         {Node(555), Node(556)},  {Node(555), Node(557)},
         {Node(556), Node(557)},  {Node(556), Node(558)},
         {Node(557), Node(558)},  {Node(557), Node(559)},
         {Node(558), Node(559)},  {Node(558), Node(560)},
         {Node(559), Node(560)},  {Node(559), Node(561)},
         {Node(560), Node(561)},  {Node(560), Node(562)},
         {Node(561), Node(562)},  {Node(561), Node(563)},
         {Node(562), Node(563)},  {Node(562), Node(564)},
         {Node(563), Node(564)},  {Node(563), Node(565)},
         {Node(564), Node(565)},  {Node(564), Node(566)},
         {Node(565), Node(566)},  {Node(565), Node(567)},
         {Node(566), Node(567)},  {Node(566), Node(568)},
         {Node(567), Node(568)},  {Node(567), Node(569)},
         {Node(568), Node(569)},  {Node(568), Node(570)},
         {Node(569), Node(570)},  {Node(569), Node(571)},
         {Node(570), Node(571)},  {Node(570), Node(572)},
         {Node(571), Node(572)},  {Node(571), Node(573)},
         {Node(572), Node(573)},  {Node(572), Node(574)},
         {Node(573), Node(574)},  {Node(573), Node(575)},
         {Node(574), Node(575)},  {Node(574), Node(576)},
         {Node(575), Node(576)},  {Node(575), Node(577)},
         {Node(576), Node(577)},  {Node(576), Node(578)},
         {Node(577), Node(578)},  {Node(577), Node(579)},
         {Node(578), Node(579)},  {Node(578), Node(580)},
         {Node(579), Node(580)},  {Node(579), Node(581)},
         {Node(580), Node(581)},  {Node(580), Node(582)},
         {Node(581), Node(582)},  {Node(581), Node(583)},
         {Node(582), Node(583)},  {Node(582), Node(584)},
         {Node(583), Node(584)},  {Node(583), Node(585)},
         {Node(584), Node(585)},  {Node(584), Node(586)},
         {Node(585), Node(586)},  {Node(585), Node(587)},
         {Node(586), Node(587)},  {Node(586), Node(588)},
         {Node(587), Node(588)},  {Node(587), Node(589)},
         {Node(588), Node(589)},  {Node(588), Node(590)},
         {Node(589), Node(590)},  {Node(589), Node(591)},
         {Node(590), Node(591)},  {Node(590), Node(592)},
         {Node(591), Node(592)},  {Node(591), Node(593)},
         {Node(592), Node(593)},  {Node(592), Node(594)},
         {Node(593), Node(594)},  {Node(593), Node(595)},
         {Node(594), Node(595)},  {Node(594), Node(596)},
         {Node(595), Node(596)},  {Node(595), Node(597)},
         {Node(596), Node(597)},  {Node(596), Node(598)},
         {Node(597), Node(598)},  {Node(597), Node(599)},
         {Node(598), Node(599)},  {Node(598), Node(600)},
         {Node(599), Node(600)},  {Node(599), Node(601)},
         {Node(600), Node(601)},  {Node(600), Node(602)},
         {Node(601), Node(602)},  {Node(601), Node(603)},
         {Node(602), Node(603)},  {Node(602), Node(604)},
         {Node(603), Node(604)},  {Node(603), Node(605)},
         {Node(604), Node(605)},  {Node(604), Node(606)},
         {Node(605), Node(606)},  {Node(605), Node(607)},
         {Node(606), Node(607)},  {Node(606), Node(608)},
         {Node(607), Node(608)},  {Node(607), Node(609)},
         {Node(608), Node(609)},  {Node(608), Node(610)},
         {Node(609), Node(610)},  {Node(609), Node(611)},
         {Node(610), Node(611)},  {Node(610), Node(612)},
         {Node(611), Node(612)},  {Node(611), Node(613)},
         {Node(612), Node(613)},  {Node(612), Node(614)},
         {Node(613), Node(614)},  {Node(613), Node(615)},
         {Node(614), Node(615)},  {Node(614), Node(616)},
         {Node(615), Node(616)},  {Node(615), Node(617)},
         {Node(616), Node(617)},  {Node(616), Node(618)},
         {Node(617), Node(618)},  {Node(617), Node(619)},
         {Node(618), Node(619)},  {Node(618), Node(620)},
         {Node(619), Node(620)},  {Node(619), Node(621)},
         {Node(620), Node(621)},  {Node(620), Node(622)},
         {Node(621), Node(622)},  {Node(621), Node(623)},
         {Node(622), Node(623)},  {Node(622), Node(624)},
         {Node(623), Node(624)},  {Node(623), Node(625)},
         {Node(624), Node(625)},  {Node(624), Node(626)},
         {Node(625), Node(626)},  {Node(625), Node(627)},
         {Node(626), Node(627)},  {Node(626), Node(628)},
         {Node(627), Node(628)},  {Node(627), Node(629)},
         {Node(628), Node(629)},  {Node(628), Node(630)},
         {Node(629), Node(630)},  {Node(629), Node(631)},
         {Node(630), Node(631)},  {Node(630), Node(632)},
         {Node(631), Node(632)},  {Node(631), Node(633)},
         {Node(632), Node(633)},  {Node(632), Node(634)},
         {Node(633), Node(634)},  {Node(633), Node(635)},
         {Node(634), Node(635)},  {Node(634), Node(636)},
         {Node(635), Node(636)},  {Node(635), Node(637)},
         {Node(636), Node(637)},  {Node(636), Node(638)},
         {Node(637), Node(638)},  {Node(637), Node(639)},
         {Node(638), Node(639)},  {Node(638), Node(640)},
         {Node(639), Node(640)},  {Node(639), Node(641)},
         {Node(640), Node(641)},  {Node(640), Node(642)},
         {Node(641), Node(642)},  {Node(641), Node(643)},
         {Node(642), Node(643)},  {Node(642), Node(644)},
         {Node(643), Node(644)},  {Node(643), Node(645)},
         {Node(644), Node(645)},  {Node(644), Node(646)},
         {Node(645), Node(646)},  {Node(645), Node(647)},
         {Node(646), Node(647)},  {Node(646), Node(648)},
         {Node(647), Node(648)},  {Node(647), Node(649)},
         {Node(648), Node(649)},  {Node(648), Node(650)},
         {Node(649), Node(650)},  {Node(649), Node(651)},
         {Node(650), Node(651)},  {Node(650), Node(652)},
         {Node(651), Node(652)},  {Node(651), Node(653)},
         {Node(652), Node(653)},  {Node(652), Node(654)},
         {Node(653), Node(654)},  {Node(653), Node(655)},
         {Node(654), Node(655)},  {Node(654), Node(656)},
         {Node(655), Node(656)},  {Node(655), Node(657)},
         {Node(656), Node(657)},  {Node(656), Node(658)},
         {Node(657), Node(658)},  {Node(657), Node(659)},
         {Node(658), Node(659)},  {Node(658), Node(660)},
         {Node(659), Node(660)},  {Node(659), Node(661)},
         {Node(660), Node(661)},  {Node(660), Node(662)},
         {Node(661), Node(662)},  {Node(661), Node(663)},
         {Node(662), Node(663)},  {Node(662), Node(664)},
         {Node(663), Node(664)},  {Node(663), Node(665)},
         {Node(664), Node(665)},  {Node(664), Node(666)},
         {Node(665), Node(666)},  {Node(665), Node(667)},
         {Node(666), Node(667)},  {Node(666), Node(668)},
         {Node(667), Node(668)},  {Node(667), Node(669)},
         {Node(668), Node(669)},  {Node(668), Node(670)},
         {Node(669), Node(670)},  {Node(669), Node(671)},
         {Node(670), Node(671)},  {Node(670), Node(672)},
         {Node(671), Node(672)},  {Node(671), Node(673)},
         {Node(672), Node(673)},  {Node(672), Node(674)},
         {Node(673), Node(674)},  {Node(673), Node(675)},
         {Node(674), Node(675)},  {Node(674), Node(676)},
         {Node(675), Node(676)},  {Node(675), Node(677)},
         {Node(676), Node(677)},  {Node(676), Node(678)},
         {Node(677), Node(678)},  {Node(677), Node(679)},
         {Node(678), Node(679)},  {Node(678), Node(680)},
         {Node(679), Node(680)},  {Node(679), Node(681)},
         {Node(680), Node(681)},  {Node(680), Node(682)},
         {Node(681), Node(682)},  {Node(681), Node(683)},
         {Node(682), Node(683)},  {Node(682), Node(684)},
         {Node(683), Node(684)},  {Node(683), Node(685)},
         {Node(684), Node(685)},  {Node(684), Node(686)},
         {Node(685), Node(686)},  {Node(685), Node(687)},
         {Node(686), Node(687)},  {Node(686), Node(688)},
         {Node(687), Node(688)},  {Node(687), Node(689)},
         {Node(688), Node(689)},  {Node(688), Node(690)},
         {Node(689), Node(690)},  {Node(689), Node(691)},
         {Node(690), Node(691)},  {Node(690), Node(692)},
         {Node(691), Node(692)},  {Node(691), Node(693)},
         {Node(692), Node(693)},  {Node(692), Node(694)},
         {Node(693), Node(694)},  {Node(693), Node(695)},
         {Node(694), Node(695)},  {Node(694), Node(696)},
         {Node(695), Node(696)},  {Node(695), Node(697)},
         {Node(696), Node(697)},  {Node(696), Node(698)},
         {Node(697), Node(698)},  {Node(697), Node(699)},
         {Node(698), Node(699)},  {Node(698), Node(700)},
         {Node(699), Node(700)},  {Node(699), Node(701)},
         {Node(700), Node(701)},  {Node(700), Node(702)},
         {Node(701), Node(702)},  {Node(701), Node(703)},
         {Node(702), Node(703)},  {Node(702), Node(704)},
         {Node(703), Node(704)},  {Node(703), Node(705)},
         {Node(704), Node(705)},  {Node(704), Node(706)},
         {Node(705), Node(706)},  {Node(705), Node(707)},
         {Node(706), Node(707)},  {Node(706), Node(708)},
         {Node(707), Node(708)},  {Node(707), Node(709)},
         {Node(708), Node(709)},  {Node(708), Node(710)},
         {Node(709), Node(710)},  {Node(709), Node(711)},
         {Node(710), Node(711)},  {Node(710), Node(712)},
         {Node(711), Node(712)},  {Node(711), Node(713)},
         {Node(712), Node(713)},  {Node(712), Node(714)},
         {Node(713), Node(714)},  {Node(713), Node(715)},
         {Node(714), Node(715)},  {Node(714), Node(716)},
         {Node(715), Node(716)},  {Node(715), Node(717)},
         {Node(716), Node(717)},  {Node(716), Node(718)},
         {Node(717), Node(718)},  {Node(717), Node(719)},
         {Node(718), Node(719)},  {Node(718), Node(720)},
         {Node(719), Node(720)},  {Node(719), Node(721)},
         {Node(720), Node(721)},  {Node(720), Node(722)},
         {Node(721), Node(722)},  {Node(721), Node(723)},
         {Node(722), Node(723)},  {Node(722), Node(724)},
         {Node(723), Node(724)},  {Node(723), Node(725)},
         {Node(724), Node(725)},  {Node(724), Node(726)},
         {Node(725), Node(726)},  {Node(725), Node(727)},
         {Node(726), Node(727)},  {Node(726), Node(728)},
         {Node(727), Node(728)},  {Node(727), Node(729)},
         {Node(728), Node(729)},  {Node(728), Node(730)},
         {Node(729), Node(730)},  {Node(729), Node(731)},
         {Node(730), Node(731)},  {Node(730), Node(732)},
         {Node(731), Node(732)},  {Node(731), Node(733)},
         {Node(732), Node(733)},  {Node(732), Node(734)},
         {Node(733), Node(734)},  {Node(733), Node(735)},
         {Node(734), Node(735)},  {Node(734), Node(736)},
         {Node(735), Node(736)},  {Node(735), Node(737)},
         {Node(736), Node(737)},  {Node(736), Node(738)},
         {Node(737), Node(738)},  {Node(737), Node(739)},
         {Node(738), Node(739)},  {Node(738), Node(740)},
         {Node(739), Node(740)},  {Node(739), Node(741)},
         {Node(740), Node(741)},  {Node(740), Node(742)},
         {Node(741), Node(742)},  {Node(741), Node(743)},
         {Node(742), Node(743)},  {Node(742), Node(744)},
         {Node(743), Node(744)},  {Node(743), Node(745)},
         {Node(744), Node(745)},  {Node(744), Node(746)},
         {Node(745), Node(746)},  {Node(745), Node(747)},
         {Node(746), Node(747)},  {Node(746), Node(748)},
         {Node(747), Node(748)},  {Node(747), Node(749)},
         {Node(748), Node(749)},  {Node(748), Node(750)},
         {Node(749), Node(750)},  {Node(749), Node(751)},
         {Node(750), Node(751)},  {Node(750), Node(752)},
         {Node(751), Node(752)},  {Node(751), Node(753)},
         {Node(752), Node(753)},  {Node(752), Node(754)},
         {Node(753), Node(754)},  {Node(753), Node(755)},
         {Node(754), Node(755)},  {Node(754), Node(756)},
         {Node(755), Node(756)},  {Node(755), Node(757)},
         {Node(756), Node(757)},  {Node(756), Node(758)},
         {Node(757), Node(758)},  {Node(757), Node(759)},
         {Node(758), Node(759)},  {Node(758), Node(760)},
         {Node(759), Node(760)},  {Node(759), Node(761)},
         {Node(760), Node(761)},  {Node(760), Node(762)},
         {Node(761), Node(762)},  {Node(761), Node(763)},
         {Node(762), Node(763)},  {Node(762), Node(764)},
         {Node(763), Node(764)},  {Node(763), Node(765)},
         {Node(764), Node(765)},  {Node(764), Node(766)},
         {Node(765), Node(766)},  {Node(765), Node(767)},
         {Node(766), Node(767)},  {Node(766), Node(768)},
         {Node(767), Node(768)},  {Node(767), Node(769)},
         {Node(768), Node(769)},  {Node(768), Node(770)},
         {Node(769), Node(770)},  {Node(769), Node(771)},
         {Node(770), Node(771)},  {Node(770), Node(772)},
         {Node(771), Node(772)},  {Node(771), Node(773)},
         {Node(772), Node(773)},  {Node(772), Node(774)},
         {Node(773), Node(774)},  {Node(773), Node(775)},
         {Node(774), Node(775)},  {Node(774), Node(776)},
         {Node(775), Node(776)},  {Node(775), Node(777)},
         {Node(776), Node(777)},  {Node(776), Node(778)},
         {Node(777), Node(778)},  {Node(777), Node(779)},
         {Node(778), Node(779)},  {Node(778), Node(780)},
         {Node(779), Node(780)},  {Node(779), Node(781)},
         {Node(780), Node(781)},  {Node(780), Node(782)},
         {Node(781), Node(782)},  {Node(781), Node(783)},
         {Node(782), Node(783)},  {Node(782), Node(784)},
         {Node(783), Node(784)},  {Node(783), Node(785)},
         {Node(784), Node(785)},  {Node(784), Node(786)},
         {Node(785), Node(786)},  {Node(785), Node(787)},
         {Node(786), Node(787)},  {Node(786), Node(788)},
         {Node(787), Node(788)},  {Node(787), Node(789)},
         {Node(788), Node(789)},  {Node(788), Node(790)},
         {Node(789), Node(790)},  {Node(789), Node(791)},
         {Node(790), Node(791)},  {Node(790), Node(792)},
         {Node(791), Node(792)},  {Node(791), Node(793)},
         {Node(792), Node(793)},  {Node(792), Node(794)},
         {Node(793), Node(794)},  {Node(793), Node(795)},
         {Node(794), Node(795)},  {Node(794), Node(796)},
         {Node(795), Node(796)},  {Node(795), Node(797)},
         {Node(796), Node(797)},  {Node(796), Node(798)},
         {Node(797), Node(798)},  {Node(797), Node(799)},
         {Node(798), Node(799)},  {Node(798), Node(800)},
         {Node(799), Node(800)},  {Node(799), Node(801)},
         {Node(800), Node(801)},  {Node(800), Node(802)},
         {Node(801), Node(802)},  {Node(801), Node(803)},
         {Node(802), Node(803)},  {Node(802), Node(804)},
         {Node(803), Node(804)},  {Node(803), Node(805)},
         {Node(804), Node(805)},  {Node(804), Node(806)},
         {Node(805), Node(806)},  {Node(805), Node(807)},
         {Node(806), Node(807)},  {Node(806), Node(808)},
         {Node(807), Node(808)},  {Node(807), Node(809)},
         {Node(808), Node(809)},  {Node(808), Node(810)},
         {Node(809), Node(810)},  {Node(809), Node(811)},
         {Node(810), Node(811)},  {Node(810), Node(812)},
         {Node(811), Node(812)},  {Node(811), Node(813)},
         {Node(812), Node(813)},  {Node(812), Node(814)},
         {Node(813), Node(814)},  {Node(813), Node(815)},
         {Node(814), Node(815)},  {Node(814), Node(816)},
         {Node(815), Node(816)},  {Node(815), Node(817)},
         {Node(816), Node(817)},  {Node(816), Node(818)},
         {Node(817), Node(818)},  {Node(817), Node(819)},
         {Node(818), Node(819)},  {Node(818), Node(820)},
         {Node(819), Node(820)},  {Node(819), Node(821)},
         {Node(820), Node(821)},  {Node(820), Node(822)},
         {Node(821), Node(822)},  {Node(821), Node(823)},
         {Node(822), Node(823)},  {Node(822), Node(824)},
         {Node(823), Node(824)},  {Node(823), Node(825)},
         {Node(824), Node(825)},  {Node(824), Node(826)},
         {Node(825), Node(826)},  {Node(825), Node(827)},
         {Node(826), Node(827)},  {Node(826), Node(828)},
         {Node(827), Node(828)},  {Node(827), Node(829)},
         {Node(828), Node(829)},  {Node(828), Node(830)},
         {Node(829), Node(830)},  {Node(829), Node(831)},
         {Node(830), Node(831)},  {Node(830), Node(832)},
         {Node(831), Node(832)},  {Node(831), Node(833)},
         {Node(832), Node(833)},  {Node(832), Node(834)},
         {Node(833), Node(834)},  {Node(833), Node(835)},
         {Node(834), Node(835)},  {Node(834), Node(836)},
         {Node(835), Node(836)},  {Node(835), Node(837)},
         {Node(836), Node(837)},  {Node(836), Node(838)},
         {Node(837), Node(838)},  {Node(837), Node(839)},
         {Node(838), Node(839)},  {Node(838), Node(840)},
         {Node(839), Node(840)},  {Node(839), Node(841)},
         {Node(840), Node(841)},  {Node(840), Node(842)},
         {Node(841), Node(842)},  {Node(841), Node(843)},
         {Node(842), Node(843)},  {Node(842), Node(844)},
         {Node(843), Node(844)},  {Node(843), Node(845)},
         {Node(844), Node(845)},  {Node(844), Node(846)},
         {Node(845), Node(846)},  {Node(845), Node(847)},
         {Node(846), Node(847)},  {Node(846), Node(848)},
         {Node(847), Node(848)},  {Node(847), Node(849)},
         {Node(848), Node(849)},  {Node(848), Node(850)},
         {Node(849), Node(850)},  {Node(849), Node(851)},
         {Node(850), Node(851)},  {Node(850), Node(852)},
         {Node(851), Node(852)},  {Node(851), Node(853)},
         {Node(852), Node(853)},  {Node(852), Node(854)},
         {Node(853), Node(854)},  {Node(853), Node(855)},
         {Node(854), Node(855)},  {Node(854), Node(856)},
         {Node(855), Node(856)},  {Node(855), Node(857)},
         {Node(856), Node(857)},  {Node(856), Node(858)},
         {Node(857), Node(858)},  {Node(857), Node(859)},
         {Node(858), Node(859)},  {Node(858), Node(860)},
         {Node(859), Node(860)},  {Node(859), Node(861)},
         {Node(860), Node(861)},  {Node(860), Node(862)},
         {Node(861), Node(862)},  {Node(861), Node(863)},
         {Node(862), Node(863)},  {Node(862), Node(864)},
         {Node(863), Node(864)},  {Node(863), Node(865)},
         {Node(864), Node(865)},  {Node(864), Node(866)},
         {Node(865), Node(866)},  {Node(865), Node(867)},
         {Node(866), Node(867)},  {Node(866), Node(868)},
         {Node(867), Node(868)},  {Node(867), Node(869)},
         {Node(868), Node(869)},  {Node(868), Node(870)},
         {Node(869), Node(870)},  {Node(869), Node(871)},
         {Node(870), Node(871)},  {Node(870), Node(872)},
         {Node(871), Node(872)},  {Node(871), Node(873)},
         {Node(872), Node(873)},  {Node(872), Node(874)},
         {Node(873), Node(874)},  {Node(873), Node(875)},
         {Node(874), Node(875)},  {Node(874), Node(876)},
         {Node(875), Node(876)},  {Node(875), Node(877)},
         {Node(876), Node(877)},  {Node(876), Node(878)},
         {Node(877), Node(878)},  {Node(877), Node(879)},
         {Node(878), Node(879)},  {Node(878), Node(880)},
         {Node(879), Node(880)},  {Node(879), Node(881)},
         {Node(880), Node(881)},  {Node(880), Node(882)},
         {Node(881), Node(882)},  {Node(881), Node(883)},
         {Node(882), Node(883)},  {Node(882), Node(884)},
         {Node(883), Node(884)},  {Node(883), Node(885)},
         {Node(884), Node(885)},  {Node(884), Node(886)},
         {Node(885), Node(886)},  {Node(885), Node(887)},
         {Node(886), Node(887)},  {Node(886), Node(888)},
         {Node(887), Node(888)},  {Node(887), Node(889)},
         {Node(888), Node(889)},  {Node(888), Node(890)},
         {Node(889), Node(890)},  {Node(889), Node(891)},
         {Node(890), Node(891)},  {Node(890), Node(892)},
         {Node(891), Node(892)},  {Node(891), Node(893)},
         {Node(892), Node(893)},  {Node(892), Node(894)},
         {Node(893), Node(894)},  {Node(893), Node(895)},
         {Node(894), Node(895)},  {Node(894), Node(896)},
         {Node(895), Node(896)},  {Node(895), Node(897)},
         {Node(896), Node(897)},  {Node(896), Node(898)},
         {Node(897), Node(898)},  {Node(897), Node(899)},
         {Node(898), Node(899)},  {Node(898), Node(900)},
         {Node(899), Node(900)},  {Node(899), Node(901)},
         {Node(900), Node(901)},  {Node(900), Node(902)},
         {Node(901), Node(902)},  {Node(901), Node(903)},
         {Node(902), Node(903)},  {Node(902), Node(904)},
         {Node(903), Node(904)},  {Node(903), Node(905)},
         {Node(904), Node(905)},  {Node(904), Node(906)},
         {Node(905), Node(906)},  {Node(905), Node(907)},
         {Node(906), Node(907)},  {Node(906), Node(908)},
         {Node(907), Node(908)},  {Node(907), Node(909)},
         {Node(908), Node(909)},  {Node(908), Node(910)},
         {Node(909), Node(910)},  {Node(909), Node(911)},
         {Node(910), Node(911)},  {Node(910), Node(912)},
         {Node(911), Node(912)},  {Node(911), Node(913)},
         {Node(912), Node(913)},  {Node(912), Node(914)},
         {Node(913), Node(914)},  {Node(913), Node(915)},
         {Node(914), Node(915)},  {Node(914), Node(916)},
         {Node(915), Node(916)},  {Node(915), Node(917)},
         {Node(916), Node(917)},  {Node(916), Node(918)},
         {Node(917), Node(918)},  {Node(917), Node(919)},
         {Node(918), Node(919)},  {Node(918), Node(920)},
         {Node(919), Node(920)},  {Node(919), Node(921)},
         {Node(920), Node(921)},  {Node(920), Node(922)},
         {Node(921), Node(922)},  {Node(921), Node(923)},
         {Node(922), Node(923)},  {Node(922), Node(924)},
         {Node(923), Node(924)},  {Node(923), Node(925)},
         {Node(924), Node(925)},  {Node(924), Node(926)},
         {Node(925), Node(926)},  {Node(925), Node(927)},
         {Node(926), Node(927)},  {Node(926), Node(928)},
         {Node(927), Node(928)},  {Node(927), Node(929)},
         {Node(928), Node(929)},  {Node(928), Node(930)},
         {Node(929), Node(930)},  {Node(929), Node(931)},
         {Node(930), Node(931)},  {Node(930), Node(932)},
         {Node(931), Node(932)},  {Node(931), Node(933)},
         {Node(932), Node(933)},  {Node(932), Node(934)},
         {Node(933), Node(934)},  {Node(933), Node(935)},
         {Node(934), Node(935)},  {Node(934), Node(936)},
         {Node(935), Node(936)},  {Node(935), Node(937)},
         {Node(936), Node(937)},  {Node(936), Node(938)},
         {Node(937), Node(938)},  {Node(937), Node(939)},
         {Node(938), Node(939)},  {Node(938), Node(940)},
         {Node(939), Node(940)},  {Node(939), Node(941)},
         {Node(940), Node(941)},  {Node(940), Node(942)},
         {Node(941), Node(942)},  {Node(941), Node(943)},
         {Node(942), Node(943)},  {Node(942), Node(944)},
         {Node(943), Node(944)},  {Node(943), Node(945)},
         {Node(944), Node(945)},  {Node(944), Node(946)},
         {Node(945), Node(946)},  {Node(945), Node(947)},
         {Node(946), Node(947)},  {Node(946), Node(948)},
         {Node(947), Node(948)},  {Node(947), Node(949)},
         {Node(948), Node(949)},  {Node(948), Node(950)},
         {Node(949), Node(950)},  {Node(949), Node(951)},
         {Node(950), Node(951)},  {Node(950), Node(952)},
         {Node(951), Node(952)},  {Node(951), Node(953)},
         {Node(952), Node(953)},  {Node(952), Node(954)},
         {Node(953), Node(954)},  {Node(953), Node(955)},
         {Node(954), Node(955)},  {Node(954), Node(956)},
         {Node(955), Node(956)},  {Node(955), Node(957)},
         {Node(956), Node(957)},  {Node(956), Node(958)},
         {Node(957), Node(958)},  {Node(957), Node(959)},
         {Node(958), Node(959)},  {Node(958), Node(960)},
         {Node(959), Node(960)},  {Node(959), Node(961)},
         {Node(960), Node(961)},  {Node(960), Node(962)},
         {Node(961), Node(962)},  {Node(961), Node(963)},
         {Node(962), Node(963)},  {Node(962), Node(964)},
         {Node(963), Node(964)},  {Node(963), Node(965)},
         {Node(964), Node(965)},  {Node(964), Node(966)},
         {Node(965), Node(966)},  {Node(965), Node(967)},
         {Node(966), Node(967)},  {Node(966), Node(968)},
         {Node(967), Node(968)},  {Node(967), Node(969)},
         {Node(968), Node(969)},  {Node(968), Node(970)},
         {Node(969), Node(970)},  {Node(969), Node(971)},
         {Node(970), Node(971)},  {Node(970), Node(972)},
         {Node(971), Node(972)},  {Node(971), Node(973)},
         {Node(972), Node(973)},  {Node(972), Node(974)},
         {Node(973), Node(974)},  {Node(973), Node(975)},
         {Node(974), Node(975)},  {Node(974), Node(976)},
         {Node(975), Node(976)},  {Node(975), Node(977)},
         {Node(976), Node(977)},  {Node(976), Node(978)},
         {Node(977), Node(978)},  {Node(977), Node(979)},
         {Node(978), Node(979)},  {Node(978), Node(980)},
         {Node(979), Node(980)},  {Node(979), Node(981)},
         {Node(980), Node(981)},  {Node(980), Node(982)},
         {Node(981), Node(982)},  {Node(981), Node(983)},
         {Node(982), Node(983)},  {Node(982), Node(984)},
         {Node(983), Node(984)},  {Node(983), Node(985)},
         {Node(984), Node(985)},  {Node(984), Node(986)},
         {Node(985), Node(986)},  {Node(985), Node(987)},
         {Node(986), Node(987)},  {Node(986), Node(988)},
         {Node(987), Node(988)},  {Node(987), Node(989)},
         {Node(988), Node(989)},  {Node(988), Node(990)},
         {Node(989), Node(990)},  {Node(989), Node(991)},
         {Node(990), Node(991)},  {Node(990), Node(992)},
         {Node(991), Node(992)},  {Node(991), Node(993)},
         {Node(992), Node(993)},  {Node(992), Node(994)},
         {Node(993), Node(994)},  {Node(993), Node(995)},
         {Node(994), Node(995)},  {Node(994), Node(996)},
         {Node(995), Node(996)},  {Node(995), Node(997)},
         {Node(996), Node(997)},  {Node(996), Node(998)},
         {Node(997), Node(998)},  {Node(997), Node(999)},
         {Node(998), Node(999)},  {Node(998), Node(1000)},
         {Node(999), Node(1000)}, {Node(999), Node(1001)}});

    aas::PathHandler handler(archi);

    aas::PathHandler handi = handler.construct_acyclic_handler();

    handi.get_distance_matrix();

    handi.get_path_matrix();

    // REQUIRE(handi.get_distance_matrix() == correct_distance_matrix_2);

    // REQUIRE(handi.get_path_matrix() == correct_path_matrix_2);
  }
  GIVEN("acyclic path generation X - omp") {
    const Architecture archi(
        {{Node(0), Node(10)},     {Node(0), Node(5)},
         {Node(0), Node(17)},     {Node(1), Node(11)},
         {Node(1), Node(6)},      {Node(1), Node(18)},
         {Node(2), Node(12)},     {Node(2), Node(7)},
         {Node(2), Node(19)},     {Node(3), Node(13)},
         {Node(3), Node(8)},      {Node(3), Node(20)},
         {Node(4), Node(14)},     {Node(4), Node(9)},
         {Node(4), Node(21)},     {Node(5), Node(15)},
         {Node(5), Node(10)},     {Node(5), Node(22)},
         {Node(6), Node(16)},     {Node(6), Node(11)},
         {Node(6), Node(23)},     {Node(7), Node(17)},
         {Node(7), Node(12)},     {Node(7), Node(24)},
         {Node(8), Node(18)},     {Node(8), Node(13)},
         {Node(8), Node(25)},     {Node(9), Node(19)},
         {Node(9), Node(14)},     {Node(9), Node(26)},
         {Node(10), Node(20)},    {Node(10), Node(15)},
         {Node(10), Node(27)},    {Node(11), Node(21)},
         {Node(11), Node(16)},    {Node(11), Node(28)},
         {Node(12), Node(22)},    {Node(12), Node(17)},
         {Node(12), Node(29)},    {Node(13), Node(23)},
         {Node(13), Node(18)},    {Node(13), Node(30)},
         {Node(14), Node(24)},    {Node(14), Node(19)},
         {Node(14), Node(31)},    {Node(15), Node(25)},
         {Node(15), Node(20)},    {Node(15), Node(32)},
         {Node(16), Node(26)},    {Node(16), Node(21)},
         {Node(16), Node(33)},    {Node(17), Node(27)},
         {Node(17), Node(22)},    {Node(17), Node(34)},
         {Node(18), Node(28)},    {Node(18), Node(23)},
         {Node(18), Node(35)},    {Node(19), Node(29)},
         {Node(19), Node(24)},    {Node(19), Node(36)},
         {Node(20), Node(30)},    {Node(20), Node(25)},
         {Node(20), Node(37)},    {Node(21), Node(31)},
         {Node(21), Node(26)},    {Node(21), Node(38)},
         {Node(22), Node(32)},    {Node(22), Node(27)},
         {Node(22), Node(39)},    {Node(23), Node(33)},
         {Node(23), Node(28)},    {Node(23), Node(40)},
         {Node(24), Node(34)},    {Node(24), Node(29)},
         {Node(24), Node(41)},    {Node(25), Node(35)},
         {Node(25), Node(30)},    {Node(25), Node(42)},
         {Node(26), Node(36)},    {Node(26), Node(31)},
         {Node(26), Node(43)},    {Node(27), Node(37)},
         {Node(27), Node(32)},    {Node(27), Node(44)},
         {Node(28), Node(38)},    {Node(28), Node(33)},
         {Node(28), Node(45)},    {Node(29), Node(39)},
         {Node(29), Node(34)},    {Node(29), Node(46)},
         {Node(30), Node(40)},    {Node(30), Node(35)},
         {Node(30), Node(47)},    {Node(31), Node(41)},
         {Node(31), Node(36)},    {Node(31), Node(48)},
         {Node(32), Node(42)},    {Node(32), Node(37)},
         {Node(32), Node(49)},    {Node(33), Node(43)},
         {Node(33), Node(38)},    {Node(33), Node(50)},
         {Node(34), Node(44)},    {Node(34), Node(39)},
         {Node(34), Node(51)},    {Node(35), Node(45)},
         {Node(35), Node(40)},    {Node(35), Node(52)},
         {Node(36), Node(46)},    {Node(36), Node(41)},
         {Node(36), Node(53)},    {Node(37), Node(47)},
         {Node(37), Node(42)},    {Node(37), Node(54)},
         {Node(38), Node(48)},    {Node(38), Node(43)},
         {Node(38), Node(55)},    {Node(39), Node(49)},
         {Node(39), Node(44)},    {Node(39), Node(56)},
         {Node(40), Node(50)},    {Node(40), Node(45)},
         {Node(40), Node(57)},    {Node(41), Node(51)},
         {Node(41), Node(46)},    {Node(41), Node(58)},
         {Node(42), Node(52)},    {Node(42), Node(47)},
         {Node(42), Node(59)},    {Node(43), Node(53)},
         {Node(43), Node(48)},    {Node(43), Node(60)},
         {Node(44), Node(54)},    {Node(44), Node(49)},
         {Node(44), Node(61)},    {Node(45), Node(55)},
         {Node(45), Node(50)},    {Node(45), Node(62)},
         {Node(46), Node(56)},    {Node(46), Node(51)},
         {Node(46), Node(63)},    {Node(47), Node(57)},
         {Node(47), Node(52)},    {Node(47), Node(64)},
         {Node(48), Node(58)},    {Node(48), Node(53)},
         {Node(48), Node(65)},    {Node(49), Node(59)},
         {Node(49), Node(54)},    {Node(49), Node(66)},
         {Node(50), Node(60)},    {Node(50), Node(55)},
         {Node(50), Node(67)},    {Node(51), Node(61)},
         {Node(51), Node(56)},    {Node(51), Node(68)},
         {Node(52), Node(62)},    {Node(52), Node(57)},
         {Node(52), Node(69)},    {Node(53), Node(63)},
         {Node(53), Node(58)},    {Node(53), Node(70)},
         {Node(54), Node(64)},    {Node(54), Node(59)},
         {Node(54), Node(71)},    {Node(55), Node(65)},
         {Node(55), Node(60)},    {Node(55), Node(72)},
         {Node(56), Node(66)},    {Node(56), Node(61)},
         {Node(56), Node(73)},    {Node(57), Node(67)},
         {Node(57), Node(62)},    {Node(57), Node(74)},
         {Node(58), Node(68)},    {Node(58), Node(63)},
         {Node(58), Node(75)},    {Node(59), Node(69)},
         {Node(59), Node(64)},    {Node(59), Node(76)},
         {Node(60), Node(70)},    {Node(60), Node(65)},
         {Node(60), Node(77)},    {Node(61), Node(71)},
         {Node(61), Node(66)},    {Node(61), Node(78)},
         {Node(62), Node(72)},    {Node(62), Node(67)},
         {Node(62), Node(79)},    {Node(63), Node(73)},
         {Node(63), Node(68)},    {Node(63), Node(80)},
         {Node(64), Node(74)},    {Node(64), Node(69)},
         {Node(64), Node(81)},    {Node(65), Node(75)},
         {Node(65), Node(70)},    {Node(65), Node(82)},
         {Node(66), Node(76)},    {Node(66), Node(71)},
         {Node(66), Node(83)},    {Node(67), Node(77)},
         {Node(67), Node(72)},    {Node(67), Node(84)},
         {Node(68), Node(78)},    {Node(68), Node(73)},
         {Node(68), Node(85)},    {Node(69), Node(79)},
         {Node(69), Node(74)},    {Node(69), Node(86)},
         {Node(70), Node(80)},    {Node(70), Node(75)},
         {Node(70), Node(87)},    {Node(71), Node(81)},
         {Node(71), Node(76)},    {Node(71), Node(88)},
         {Node(72), Node(82)},    {Node(72), Node(77)},
         {Node(72), Node(89)},    {Node(73), Node(83)},
         {Node(73), Node(78)},    {Node(73), Node(90)},
         {Node(74), Node(84)},    {Node(74), Node(79)},
         {Node(74), Node(91)},    {Node(75), Node(85)},
         {Node(75), Node(80)},    {Node(75), Node(92)},
         {Node(76), Node(86)},    {Node(76), Node(81)},
         {Node(76), Node(93)},    {Node(77), Node(87)},
         {Node(77), Node(82)},    {Node(77), Node(94)},
         {Node(78), Node(88)},    {Node(78), Node(83)},
         {Node(78), Node(95)},    {Node(79), Node(89)},
         {Node(79), Node(84)},    {Node(79), Node(96)},
         {Node(80), Node(90)},    {Node(80), Node(85)},
         {Node(80), Node(97)},    {Node(81), Node(91)},
         {Node(81), Node(86)},    {Node(81), Node(98)},
         {Node(82), Node(92)},    {Node(82), Node(87)},
         {Node(82), Node(99)},    {Node(83), Node(93)},
         {Node(83), Node(88)},    {Node(83), Node(100)},
         {Node(84), Node(94)},    {Node(84), Node(89)},
         {Node(84), Node(101)},   {Node(85), Node(95)},
         {Node(85), Node(90)},    {Node(85), Node(102)},
         {Node(86), Node(96)},    {Node(86), Node(91)},
         {Node(86), Node(103)},   {Node(87), Node(97)},
         {Node(87), Node(92)},    {Node(87), Node(104)},
         {Node(88), Node(98)},    {Node(88), Node(93)},
         {Node(88), Node(105)},   {Node(89), Node(99)},
         {Node(89), Node(94)},    {Node(89), Node(106)},
         {Node(90), Node(100)},   {Node(90), Node(95)},
         {Node(90), Node(107)},   {Node(91), Node(101)},
         {Node(91), Node(96)},    {Node(91), Node(108)},
         {Node(92), Node(102)},   {Node(92), Node(97)},
         {Node(92), Node(109)},   {Node(93), Node(103)},
         {Node(93), Node(98)},    {Node(93), Node(110)},
         {Node(94), Node(104)},   {Node(94), Node(99)},
         {Node(94), Node(111)},   {Node(95), Node(105)},
         {Node(95), Node(100)},   {Node(95), Node(112)},
         {Node(96), Node(106)},   {Node(96), Node(101)},
         {Node(96), Node(113)},   {Node(97), Node(107)},
         {Node(97), Node(102)},   {Node(97), Node(114)},
         {Node(98), Node(108)},   {Node(98), Node(103)},
         {Node(98), Node(115)},   {Node(99), Node(109)},
         {Node(99), Node(104)},   {Node(99), Node(116)},
         {Node(100), Node(110)},  {Node(100), Node(105)},
         {Node(100), Node(117)},  {Node(101), Node(111)},
         {Node(101), Node(106)},  {Node(101), Node(118)},
         {Node(102), Node(112)},  {Node(102), Node(107)},
         {Node(102), Node(119)},  {Node(103), Node(113)},
         {Node(103), Node(108)},  {Node(103), Node(120)},
         {Node(104), Node(114)},  {Node(104), Node(109)},
         {Node(104), Node(121)},  {Node(105), Node(115)},
         {Node(105), Node(110)},  {Node(105), Node(122)},
         {Node(106), Node(116)},  {Node(106), Node(111)},
         {Node(106), Node(123)},  {Node(107), Node(117)},
         {Node(107), Node(112)},  {Node(107), Node(124)},
         {Node(108), Node(118)},  {Node(108), Node(113)},
         {Node(108), Node(125)},  {Node(109), Node(119)},
         {Node(109), Node(114)},  {Node(109), Node(126)},
         {Node(110), Node(120)},  {Node(110), Node(115)},
         {Node(110), Node(127)},  {Node(111), Node(121)},
         {Node(111), Node(116)},  {Node(111), Node(128)},
         {Node(112), Node(122)},  {Node(112), Node(117)},
         {Node(112), Node(129)},  {Node(113), Node(123)},
         {Node(113), Node(118)},  {Node(113), Node(130)},
         {Node(114), Node(124)},  {Node(114), Node(119)},
         {Node(114), Node(131)},  {Node(115), Node(125)},
         {Node(115), Node(120)},  {Node(115), Node(132)},
         {Node(116), Node(126)},  {Node(116), Node(121)},
         {Node(116), Node(133)},  {Node(117), Node(127)},
         {Node(117), Node(122)},  {Node(117), Node(134)},
         {Node(118), Node(128)},  {Node(118), Node(123)},
         {Node(118), Node(135)},  {Node(119), Node(129)},
         {Node(119), Node(124)},  {Node(119), Node(136)},
         {Node(120), Node(130)},  {Node(120), Node(125)},
         {Node(120), Node(137)},  {Node(121), Node(131)},
         {Node(121), Node(126)},  {Node(121), Node(138)},
         {Node(122), Node(132)},  {Node(122), Node(127)},
         {Node(122), Node(139)},  {Node(123), Node(133)},
         {Node(123), Node(128)},  {Node(123), Node(140)},
         {Node(124), Node(134)},  {Node(124), Node(129)},
         {Node(124), Node(141)},  {Node(125), Node(135)},
         {Node(125), Node(130)},  {Node(125), Node(142)},
         {Node(126), Node(136)},  {Node(126), Node(131)},
         {Node(126), Node(143)},  {Node(127), Node(137)},
         {Node(127), Node(132)},  {Node(127), Node(144)},
         {Node(128), Node(138)},  {Node(128), Node(133)},
         {Node(128), Node(145)},  {Node(129), Node(139)},
         {Node(129), Node(134)},  {Node(129), Node(146)},
         {Node(130), Node(140)},  {Node(130), Node(135)},
         {Node(130), Node(147)},  {Node(131), Node(141)},
         {Node(131), Node(136)},  {Node(131), Node(148)},
         {Node(132), Node(142)},  {Node(132), Node(137)},
         {Node(132), Node(149)},  {Node(133), Node(143)},
         {Node(133), Node(138)},  {Node(133), Node(150)},
         {Node(134), Node(144)},  {Node(134), Node(139)},
         {Node(134), Node(151)},  {Node(135), Node(145)},
         {Node(135), Node(140)},  {Node(135), Node(152)},
         {Node(136), Node(146)},  {Node(136), Node(141)},
         {Node(136), Node(153)},  {Node(137), Node(147)},
         {Node(137), Node(142)},  {Node(137), Node(154)},
         {Node(138), Node(148)},  {Node(138), Node(143)},
         {Node(138), Node(155)},  {Node(139), Node(149)},
         {Node(139), Node(144)},  {Node(139), Node(156)},
         {Node(140), Node(150)},  {Node(140), Node(145)},
         {Node(140), Node(157)},  {Node(141), Node(151)},
         {Node(141), Node(146)},  {Node(141), Node(158)},
         {Node(142), Node(152)},  {Node(142), Node(147)},
         {Node(142), Node(159)},  {Node(143), Node(153)},
         {Node(143), Node(148)},  {Node(143), Node(160)},
         {Node(144), Node(154)},  {Node(144), Node(149)},
         {Node(144), Node(161)},  {Node(145), Node(155)},
         {Node(145), Node(150)},  {Node(145), Node(162)},
         {Node(146), Node(156)},  {Node(146), Node(151)},
         {Node(146), Node(163)},  {Node(147), Node(157)},
         {Node(147), Node(152)},  {Node(147), Node(164)},
         {Node(148), Node(158)},  {Node(148), Node(153)},
         {Node(148), Node(165)},  {Node(149), Node(159)},
         {Node(149), Node(154)},  {Node(149), Node(166)},
         {Node(150), Node(160)},  {Node(150), Node(155)},
         {Node(150), Node(167)},  {Node(151), Node(161)},
         {Node(151), Node(156)},  {Node(151), Node(168)},
         {Node(152), Node(162)},  {Node(152), Node(157)},
         {Node(152), Node(169)},  {Node(153), Node(163)},
         {Node(153), Node(158)},  {Node(153), Node(170)},
         {Node(154), Node(164)},  {Node(154), Node(159)},
         {Node(154), Node(171)},  {Node(155), Node(165)},
         {Node(155), Node(160)},  {Node(155), Node(172)},
         {Node(156), Node(166)},  {Node(156), Node(161)},
         {Node(156), Node(173)},  {Node(157), Node(167)},
         {Node(157), Node(162)},  {Node(157), Node(174)},
         {Node(158), Node(168)},  {Node(158), Node(163)},
         {Node(158), Node(175)},  {Node(159), Node(169)},
         {Node(159), Node(164)},  {Node(159), Node(176)},
         {Node(160), Node(170)},  {Node(160), Node(165)},
         {Node(160), Node(177)},  {Node(161), Node(171)},
         {Node(161), Node(166)},  {Node(161), Node(178)},
         {Node(162), Node(172)},  {Node(162), Node(167)},
         {Node(162), Node(179)},  {Node(163), Node(173)},
         {Node(163), Node(168)},  {Node(163), Node(180)},
         {Node(164), Node(174)},  {Node(164), Node(169)},
         {Node(164), Node(181)},  {Node(165), Node(175)},
         {Node(165), Node(170)},  {Node(165), Node(182)},
         {Node(166), Node(176)},  {Node(166), Node(171)},
         {Node(166), Node(183)},  {Node(167), Node(177)},
         {Node(167), Node(172)},  {Node(167), Node(184)},
         {Node(168), Node(178)},  {Node(168), Node(173)},
         {Node(168), Node(185)},  {Node(169), Node(179)},
         {Node(169), Node(174)},  {Node(169), Node(186)},
         {Node(170), Node(180)},  {Node(170), Node(175)},
         {Node(170), Node(187)},  {Node(171), Node(181)},
         {Node(171), Node(176)},  {Node(171), Node(188)},
         {Node(172), Node(182)},  {Node(172), Node(177)},
         {Node(172), Node(189)},  {Node(173), Node(183)},
         {Node(173), Node(178)},  {Node(173), Node(190)},
         {Node(174), Node(184)},  {Node(174), Node(179)},
         {Node(174), Node(191)},  {Node(175), Node(185)},
         {Node(175), Node(180)},  {Node(175), Node(192)},
         {Node(176), Node(186)},  {Node(176), Node(181)},
         {Node(176), Node(193)},  {Node(177), Node(187)},
         {Node(177), Node(182)},  {Node(177), Node(194)},
         {Node(178), Node(188)},  {Node(178), Node(183)},
         {Node(178), Node(195)},  {Node(179), Node(189)},
         {Node(179), Node(184)},  {Node(179), Node(196)},
         {Node(180), Node(190)},  {Node(180), Node(185)},
         {Node(180), Node(197)},  {Node(181), Node(191)},
         {Node(181), Node(186)},  {Node(181), Node(198)},
         {Node(182), Node(192)},  {Node(182), Node(187)},
         {Node(182), Node(199)},  {Node(183), Node(193)},
         {Node(183), Node(188)},  {Node(183), Node(200)},
         {Node(184), Node(194)},  {Node(184), Node(189)},
         {Node(184), Node(201)},  {Node(185), Node(195)},
         {Node(185), Node(190)},  {Node(185), Node(202)},
         {Node(186), Node(196)},  {Node(186), Node(191)},
         {Node(186), Node(203)},  {Node(187), Node(197)},
         {Node(187), Node(192)},  {Node(187), Node(204)},
         {Node(188), Node(198)},  {Node(188), Node(193)},
         {Node(188), Node(205)},  {Node(189), Node(199)},
         {Node(189), Node(194)},  {Node(189), Node(206)},
         {Node(190), Node(200)},  {Node(190), Node(195)},
         {Node(190), Node(207)},  {Node(191), Node(201)},
         {Node(191), Node(196)},  {Node(191), Node(208)},
         {Node(192), Node(202)},  {Node(192), Node(197)},
         {Node(192), Node(209)},  {Node(193), Node(203)},
         {Node(193), Node(198)},  {Node(193), Node(210)},
         {Node(194), Node(204)},  {Node(194), Node(199)},
         {Node(194), Node(211)},  {Node(195), Node(205)},
         {Node(195), Node(200)},  {Node(195), Node(212)},
         {Node(196), Node(206)},  {Node(196), Node(201)},
         {Node(196), Node(213)},  {Node(197), Node(207)},
         {Node(197), Node(202)},  {Node(197), Node(214)},
         {Node(198), Node(208)},  {Node(198), Node(203)},
         {Node(198), Node(215)},  {Node(199), Node(209)},
         {Node(199), Node(204)},  {Node(199), Node(216)},
         {Node(200), Node(210)},  {Node(200), Node(205)},
         {Node(200), Node(217)},  {Node(201), Node(211)},
         {Node(201), Node(206)},  {Node(201), Node(218)},
         {Node(202), Node(212)},  {Node(202), Node(207)},
         {Node(202), Node(219)},  {Node(203), Node(213)},
         {Node(203), Node(208)},  {Node(203), Node(220)},
         {Node(204), Node(214)},  {Node(204), Node(209)},
         {Node(204), Node(221)},  {Node(205), Node(215)},
         {Node(205), Node(210)},  {Node(205), Node(222)},
         {Node(206), Node(216)},  {Node(206), Node(211)},
         {Node(206), Node(223)},  {Node(207), Node(217)},
         {Node(207), Node(212)},  {Node(207), Node(224)},
         {Node(208), Node(218)},  {Node(208), Node(213)},
         {Node(208), Node(225)},  {Node(209), Node(219)},
         {Node(209), Node(214)},  {Node(209), Node(226)},
         {Node(210), Node(220)},  {Node(210), Node(215)},
         {Node(210), Node(227)},  {Node(211), Node(221)},
         {Node(211), Node(216)},  {Node(211), Node(228)},
         {Node(212), Node(222)},  {Node(212), Node(217)},
         {Node(212), Node(229)},  {Node(213), Node(223)},
         {Node(213), Node(218)},  {Node(213), Node(230)},
         {Node(214), Node(224)},  {Node(214), Node(219)},
         {Node(214), Node(231)},  {Node(215), Node(225)},
         {Node(215), Node(220)},  {Node(215), Node(232)},
         {Node(216), Node(226)},  {Node(216), Node(221)},
         {Node(216), Node(233)},  {Node(217), Node(227)},
         {Node(217), Node(222)},  {Node(217), Node(234)},
         {Node(218), Node(228)},  {Node(218), Node(223)},
         {Node(218), Node(235)},  {Node(219), Node(229)},
         {Node(219), Node(224)},  {Node(219), Node(236)},
         {Node(220), Node(230)},  {Node(220), Node(225)},
         {Node(220), Node(237)},  {Node(221), Node(231)},
         {Node(221), Node(226)},  {Node(221), Node(238)},
         {Node(222), Node(232)},  {Node(222), Node(227)},
         {Node(222), Node(239)},  {Node(223), Node(233)},
         {Node(223), Node(228)},  {Node(223), Node(240)},
         {Node(224), Node(234)},  {Node(224), Node(229)},
         {Node(224), Node(241)},  {Node(225), Node(235)},
         {Node(225), Node(230)},  {Node(225), Node(242)},
         {Node(226), Node(236)},  {Node(226), Node(231)},
         {Node(226), Node(243)},  {Node(227), Node(237)},
         {Node(227), Node(232)},  {Node(227), Node(244)},
         {Node(228), Node(238)},  {Node(228), Node(233)},
         {Node(228), Node(245)},  {Node(229), Node(239)},
         {Node(229), Node(234)},  {Node(229), Node(246)},
         {Node(230), Node(240)},  {Node(230), Node(235)},
         {Node(230), Node(247)},  {Node(231), Node(241)},
         {Node(231), Node(236)},  {Node(231), Node(248)},
         {Node(232), Node(242)},  {Node(232), Node(237)},
         {Node(232), Node(249)},  {Node(233), Node(243)},
         {Node(233), Node(238)},  {Node(233), Node(250)},
         {Node(234), Node(244)},  {Node(234), Node(239)},
         {Node(234), Node(251)},  {Node(235), Node(245)},
         {Node(235), Node(240)},  {Node(235), Node(252)},
         {Node(236), Node(246)},  {Node(236), Node(241)},
         {Node(236), Node(253)},  {Node(237), Node(247)},
         {Node(237), Node(242)},  {Node(237), Node(254)},
         {Node(238), Node(248)},  {Node(238), Node(243)},
         {Node(238), Node(255)},  {Node(239), Node(249)},
         {Node(239), Node(244)},  {Node(239), Node(256)},
         {Node(240), Node(250)},  {Node(240), Node(245)},
         {Node(240), Node(257)},  {Node(241), Node(251)},
         {Node(241), Node(246)},  {Node(241), Node(258)},
         {Node(242), Node(252)},  {Node(242), Node(247)},
         {Node(242), Node(259)},  {Node(243), Node(253)},
         {Node(243), Node(248)},  {Node(243), Node(260)},
         {Node(244), Node(254)},  {Node(244), Node(249)},
         {Node(244), Node(261)},  {Node(245), Node(255)},
         {Node(245), Node(250)},  {Node(245), Node(262)},
         {Node(246), Node(256)},  {Node(246), Node(251)},
         {Node(246), Node(263)},  {Node(247), Node(257)},
         {Node(247), Node(252)},  {Node(247), Node(264)},
         {Node(248), Node(258)},  {Node(248), Node(253)},
         {Node(248), Node(265)},  {Node(249), Node(259)},
         {Node(249), Node(254)},  {Node(249), Node(266)},
         {Node(250), Node(260)},  {Node(250), Node(255)},
         {Node(250), Node(267)},  {Node(251), Node(261)},
         {Node(251), Node(256)},  {Node(251), Node(268)},
         {Node(252), Node(262)},  {Node(252), Node(257)},
         {Node(252), Node(269)},  {Node(253), Node(263)},
         {Node(253), Node(258)},  {Node(253), Node(270)},
         {Node(254), Node(264)},  {Node(254), Node(259)},
         {Node(254), Node(271)},  {Node(255), Node(265)},
         {Node(255), Node(260)},  {Node(255), Node(272)},
         {Node(256), Node(266)},  {Node(256), Node(261)},
         {Node(256), Node(273)},  {Node(257), Node(267)},
         {Node(257), Node(262)},  {Node(257), Node(274)},
         {Node(258), Node(268)},  {Node(258), Node(263)},
         {Node(258), Node(275)},  {Node(259), Node(269)},
         {Node(259), Node(264)},  {Node(259), Node(276)},
         {Node(260), Node(270)},  {Node(260), Node(265)},
         {Node(260), Node(277)},  {Node(261), Node(271)},
         {Node(261), Node(266)},  {Node(261), Node(278)},
         {Node(262), Node(272)},  {Node(262), Node(267)},
         {Node(262), Node(279)},  {Node(263), Node(273)},
         {Node(263), Node(268)},  {Node(263), Node(280)},
         {Node(264), Node(274)},  {Node(264), Node(269)},
         {Node(264), Node(281)},  {Node(265), Node(275)},
         {Node(265), Node(270)},  {Node(265), Node(282)},
         {Node(266), Node(276)},  {Node(266), Node(271)},
         {Node(266), Node(283)},  {Node(267), Node(277)},
         {Node(267), Node(272)},  {Node(267), Node(284)},
         {Node(268), Node(278)},  {Node(268), Node(273)},
         {Node(268), Node(285)},  {Node(269), Node(279)},
         {Node(269), Node(274)},  {Node(269), Node(286)},
         {Node(270), Node(280)},  {Node(270), Node(275)},
         {Node(270), Node(287)},  {Node(271), Node(281)},
         {Node(271), Node(276)},  {Node(271), Node(288)},
         {Node(272), Node(282)},  {Node(272), Node(277)},
         {Node(272), Node(289)},  {Node(273), Node(283)},
         {Node(273), Node(278)},  {Node(273), Node(290)},
         {Node(274), Node(284)},  {Node(274), Node(279)},
         {Node(274), Node(291)},  {Node(275), Node(285)},
         {Node(275), Node(280)},  {Node(275), Node(292)},
         {Node(276), Node(286)},  {Node(276), Node(281)},
         {Node(276), Node(293)},  {Node(277), Node(287)},
         {Node(277), Node(282)},  {Node(277), Node(294)},
         {Node(278), Node(288)},  {Node(278), Node(283)},
         {Node(278), Node(295)},  {Node(279), Node(289)},
         {Node(279), Node(284)},  {Node(279), Node(296)},
         {Node(280), Node(290)},  {Node(280), Node(285)},
         {Node(280), Node(297)},  {Node(281), Node(291)},
         {Node(281), Node(286)},  {Node(281), Node(298)},
         {Node(282), Node(292)},  {Node(282), Node(287)},
         {Node(282), Node(299)},  {Node(283), Node(293)},
         {Node(283), Node(288)},  {Node(283), Node(300)},
         {Node(284), Node(294)},  {Node(284), Node(289)},
         {Node(284), Node(301)},  {Node(285), Node(295)},
         {Node(285), Node(290)},  {Node(285), Node(302)},
         {Node(286), Node(296)},  {Node(286), Node(291)},
         {Node(286), Node(303)},  {Node(287), Node(297)},
         {Node(287), Node(292)},  {Node(287), Node(304)},
         {Node(288), Node(298)},  {Node(288), Node(293)},
         {Node(288), Node(305)},  {Node(289), Node(299)},
         {Node(289), Node(294)},  {Node(289), Node(306)},
         {Node(290), Node(300)},  {Node(290), Node(295)},
         {Node(290), Node(307)},  {Node(291), Node(301)},
         {Node(291), Node(296)},  {Node(291), Node(308)},
         {Node(292), Node(302)},  {Node(292), Node(297)},
         {Node(292), Node(309)},  {Node(293), Node(303)},
         {Node(293), Node(298)},  {Node(293), Node(310)},
         {Node(294), Node(304)},  {Node(294), Node(299)},
         {Node(294), Node(311)},  {Node(295), Node(305)},
         {Node(295), Node(300)},  {Node(295), Node(312)},
         {Node(296), Node(306)},  {Node(296), Node(301)},
         {Node(296), Node(313)},  {Node(297), Node(307)},
         {Node(297), Node(302)},  {Node(297), Node(314)},
         {Node(298), Node(308)},  {Node(298), Node(303)},
         {Node(298), Node(315)},  {Node(299), Node(309)},
         {Node(299), Node(304)},  {Node(299), Node(316)},
         {Node(300), Node(310)},  {Node(300), Node(305)},
         {Node(300), Node(317)},  {Node(301), Node(311)},
         {Node(301), Node(306)},  {Node(301), Node(318)},
         {Node(302), Node(312)},  {Node(302), Node(307)},
         {Node(302), Node(319)},  {Node(303), Node(313)},
         {Node(303), Node(308)},  {Node(303), Node(320)},
         {Node(304), Node(314)},  {Node(304), Node(309)},
         {Node(304), Node(321)},  {Node(305), Node(315)},
         {Node(305), Node(310)},  {Node(305), Node(322)},
         {Node(306), Node(316)},  {Node(306), Node(311)},
         {Node(306), Node(323)},  {Node(307), Node(317)},
         {Node(307), Node(312)},  {Node(307), Node(324)},
         {Node(308), Node(318)},  {Node(308), Node(313)},
         {Node(308), Node(325)},  {Node(309), Node(319)},
         {Node(309), Node(314)},  {Node(309), Node(326)},
         {Node(310), Node(320)},  {Node(310), Node(315)},
         {Node(310), Node(327)},  {Node(311), Node(321)},
         {Node(311), Node(316)},  {Node(311), Node(328)},
         {Node(312), Node(322)},  {Node(312), Node(317)},
         {Node(312), Node(329)},  {Node(313), Node(323)},
         {Node(313), Node(318)},  {Node(313), Node(330)},
         {Node(314), Node(324)},  {Node(314), Node(319)},
         {Node(314), Node(331)},  {Node(315), Node(325)},
         {Node(315), Node(320)},  {Node(315), Node(332)},
         {Node(316), Node(326)},  {Node(316), Node(321)},
         {Node(316), Node(333)},  {Node(317), Node(327)},
         {Node(317), Node(322)},  {Node(317), Node(334)},
         {Node(318), Node(328)},  {Node(318), Node(323)},
         {Node(318), Node(335)},  {Node(319), Node(329)},
         {Node(319), Node(324)},  {Node(319), Node(336)},
         {Node(320), Node(330)},  {Node(320), Node(325)},
         {Node(320), Node(337)},  {Node(321), Node(331)},
         {Node(321), Node(326)},  {Node(321), Node(338)},
         {Node(322), Node(332)},  {Node(322), Node(327)},
         {Node(322), Node(339)},  {Node(323), Node(333)},
         {Node(323), Node(328)},  {Node(323), Node(340)},
         {Node(324), Node(334)},  {Node(324), Node(329)},
         {Node(324), Node(341)},  {Node(325), Node(335)},
         {Node(325), Node(330)},  {Node(325), Node(342)},
         {Node(326), Node(336)},  {Node(326), Node(331)},
         {Node(326), Node(343)},  {Node(327), Node(337)},
         {Node(327), Node(332)},  {Node(327), Node(344)},
         {Node(328), Node(338)},  {Node(328), Node(333)},
         {Node(328), Node(345)},  {Node(329), Node(339)},
         {Node(329), Node(334)},  {Node(329), Node(346)},
         {Node(330), Node(340)},  {Node(330), Node(335)},
         {Node(330), Node(347)},  {Node(331), Node(341)},
         {Node(331), Node(336)},  {Node(331), Node(348)},
         {Node(332), Node(342)},  {Node(332), Node(337)},
         {Node(332), Node(349)},  {Node(333), Node(343)},
         {Node(333), Node(338)},  {Node(333), Node(350)},
         {Node(334), Node(344)},  {Node(334), Node(339)},
         {Node(334), Node(351)},  {Node(335), Node(345)},
         {Node(335), Node(340)},  {Node(335), Node(352)},
         {Node(336), Node(346)},  {Node(336), Node(341)},
         {Node(336), Node(353)},  {Node(337), Node(347)},
         {Node(337), Node(342)},  {Node(337), Node(354)},
         {Node(338), Node(348)},  {Node(338), Node(343)},
         {Node(338), Node(355)},  {Node(339), Node(349)},
         {Node(339), Node(344)},  {Node(339), Node(356)},
         {Node(340), Node(350)},  {Node(340), Node(345)},
         {Node(340), Node(357)},  {Node(341), Node(351)},
         {Node(341), Node(346)},  {Node(341), Node(358)},
         {Node(342), Node(352)},  {Node(342), Node(347)},
         {Node(342), Node(359)},  {Node(343), Node(353)},
         {Node(343), Node(348)},  {Node(343), Node(360)},
         {Node(344), Node(354)},  {Node(344), Node(349)},
         {Node(344), Node(361)},  {Node(345), Node(355)},
         {Node(345), Node(350)},  {Node(345), Node(362)},
         {Node(346), Node(356)},  {Node(346), Node(351)},
         {Node(346), Node(363)},  {Node(347), Node(357)},
         {Node(347), Node(352)},  {Node(347), Node(364)},
         {Node(348), Node(358)},  {Node(348), Node(353)},
         {Node(348), Node(365)},  {Node(349), Node(359)},
         {Node(349), Node(354)},  {Node(349), Node(366)},
         {Node(350), Node(360)},  {Node(350), Node(355)},
         {Node(350), Node(367)},  {Node(351), Node(361)},
         {Node(351), Node(356)},  {Node(351), Node(368)},
         {Node(352), Node(362)},  {Node(352), Node(357)},
         {Node(352), Node(369)},  {Node(353), Node(363)},
         {Node(353), Node(358)},  {Node(353), Node(370)},
         {Node(354), Node(364)},  {Node(354), Node(359)},
         {Node(354), Node(371)},  {Node(355), Node(365)},
         {Node(355), Node(360)},  {Node(355), Node(372)},
         {Node(356), Node(366)},  {Node(356), Node(361)},
         {Node(356), Node(373)},  {Node(357), Node(367)},
         {Node(357), Node(362)},  {Node(357), Node(374)},
         {Node(358), Node(368)},  {Node(358), Node(363)},
         {Node(358), Node(375)},  {Node(359), Node(369)},
         {Node(359), Node(364)},  {Node(359), Node(376)},
         {Node(360), Node(370)},  {Node(360), Node(365)},
         {Node(360), Node(377)},  {Node(361), Node(371)},
         {Node(361), Node(366)},  {Node(361), Node(378)},
         {Node(362), Node(372)},  {Node(362), Node(367)},
         {Node(362), Node(379)},  {Node(363), Node(373)},
         {Node(363), Node(368)},  {Node(363), Node(380)},
         {Node(364), Node(374)},  {Node(364), Node(369)},
         {Node(364), Node(381)},  {Node(365), Node(375)},
         {Node(365), Node(370)},  {Node(365), Node(382)},
         {Node(366), Node(376)},  {Node(366), Node(371)},
         {Node(366), Node(383)},  {Node(367), Node(377)},
         {Node(367), Node(372)},  {Node(367), Node(384)},
         {Node(368), Node(378)},  {Node(368), Node(373)},
         {Node(368), Node(385)},  {Node(369), Node(379)},
         {Node(369), Node(374)},  {Node(369), Node(386)},
         {Node(370), Node(380)},  {Node(370), Node(375)},
         {Node(370), Node(387)},  {Node(371), Node(381)},
         {Node(371), Node(376)},  {Node(371), Node(388)},
         {Node(372), Node(382)},  {Node(372), Node(377)},
         {Node(372), Node(389)},  {Node(373), Node(383)},
         {Node(373), Node(378)},  {Node(373), Node(390)},
         {Node(374), Node(384)},  {Node(374), Node(379)},
         {Node(374), Node(391)},  {Node(375), Node(385)},
         {Node(375), Node(380)},  {Node(375), Node(392)},
         {Node(376), Node(386)},  {Node(376), Node(381)},
         {Node(376), Node(393)},  {Node(377), Node(387)},
         {Node(377), Node(382)},  {Node(377), Node(394)},
         {Node(378), Node(388)},  {Node(378), Node(383)},
         {Node(378), Node(395)},  {Node(379), Node(389)},
         {Node(379), Node(384)},  {Node(379), Node(396)},
         {Node(380), Node(390)},  {Node(380), Node(385)},
         {Node(380), Node(397)},  {Node(381), Node(391)},
         {Node(381), Node(386)},  {Node(381), Node(398)},
         {Node(382), Node(392)},  {Node(382), Node(387)},
         {Node(382), Node(399)},  {Node(383), Node(393)},
         {Node(383), Node(388)},  {Node(383), Node(400)},
         {Node(384), Node(394)},  {Node(384), Node(389)},
         {Node(384), Node(401)},  {Node(385), Node(395)},
         {Node(385), Node(390)},  {Node(385), Node(402)},
         {Node(386), Node(396)},  {Node(386), Node(391)},
         {Node(386), Node(403)},  {Node(387), Node(397)},
         {Node(387), Node(392)},  {Node(387), Node(404)},
         {Node(388), Node(398)},  {Node(388), Node(393)},
         {Node(388), Node(405)},  {Node(389), Node(399)},
         {Node(389), Node(394)},  {Node(389), Node(406)},
         {Node(390), Node(400)},  {Node(390), Node(395)},
         {Node(390), Node(407)},  {Node(391), Node(401)},
         {Node(391), Node(396)},  {Node(391), Node(408)},
         {Node(392), Node(402)},  {Node(392), Node(397)},
         {Node(392), Node(409)},  {Node(393), Node(403)},
         {Node(393), Node(398)},  {Node(393), Node(410)},
         {Node(394), Node(404)},  {Node(394), Node(399)},
         {Node(394), Node(411)},  {Node(395), Node(405)},
         {Node(395), Node(400)},  {Node(395), Node(412)},
         {Node(396), Node(406)},  {Node(396), Node(401)},
         {Node(396), Node(413)},  {Node(397), Node(407)},
         {Node(397), Node(402)},  {Node(397), Node(414)},
         {Node(398), Node(408)},  {Node(398), Node(403)},
         {Node(398), Node(415)},  {Node(399), Node(409)},
         {Node(399), Node(404)},  {Node(399), Node(416)},
         {Node(400), Node(410)},  {Node(400), Node(405)},
         {Node(400), Node(417)},  {Node(401), Node(411)},
         {Node(401), Node(406)},  {Node(401), Node(418)},
         {Node(402), Node(412)},  {Node(402), Node(407)},
         {Node(402), Node(419)},  {Node(403), Node(413)},
         {Node(403), Node(408)},  {Node(403), Node(420)},
         {Node(404), Node(414)},  {Node(404), Node(409)},
         {Node(404), Node(421)},  {Node(405), Node(415)},
         {Node(405), Node(410)},  {Node(405), Node(422)},
         {Node(406), Node(416)},  {Node(406), Node(411)},
         {Node(406), Node(423)},  {Node(407), Node(417)},
         {Node(407), Node(412)},  {Node(407), Node(424)},
         {Node(408), Node(418)},  {Node(408), Node(413)},
         {Node(408), Node(425)},  {Node(409), Node(419)},
         {Node(409), Node(414)},  {Node(409), Node(426)},
         {Node(410), Node(420)},  {Node(410), Node(415)},
         {Node(410), Node(427)},  {Node(411), Node(421)},
         {Node(411), Node(416)},  {Node(411), Node(428)},
         {Node(412), Node(422)},  {Node(412), Node(417)},
         {Node(412), Node(429)},  {Node(413), Node(423)},
         {Node(413), Node(418)},  {Node(413), Node(430)},
         {Node(414), Node(424)},  {Node(414), Node(419)},
         {Node(414), Node(431)},  {Node(415), Node(425)},
         {Node(415), Node(420)},  {Node(415), Node(432)},
         {Node(416), Node(426)},  {Node(416), Node(421)},
         {Node(416), Node(433)},  {Node(417), Node(427)},
         {Node(417), Node(422)},  {Node(417), Node(434)},
         {Node(418), Node(428)},  {Node(418), Node(423)},
         {Node(418), Node(435)},  {Node(419), Node(429)},
         {Node(419), Node(424)},  {Node(419), Node(436)},
         {Node(420), Node(430)},  {Node(420), Node(425)},
         {Node(420), Node(437)},  {Node(421), Node(431)},
         {Node(421), Node(426)},  {Node(421), Node(438)},
         {Node(422), Node(432)},  {Node(422), Node(427)},
         {Node(422), Node(439)},  {Node(423), Node(433)},
         {Node(423), Node(428)},  {Node(423), Node(440)},
         {Node(424), Node(434)},  {Node(424), Node(429)},
         {Node(424), Node(441)},  {Node(425), Node(435)},
         {Node(425), Node(430)},  {Node(425), Node(442)},
         {Node(426), Node(436)},  {Node(426), Node(431)},
         {Node(426), Node(443)},  {Node(427), Node(437)},
         {Node(427), Node(432)},  {Node(427), Node(444)},
         {Node(428), Node(438)},  {Node(428), Node(433)},
         {Node(428), Node(445)},  {Node(429), Node(439)},
         {Node(429), Node(434)},  {Node(429), Node(446)},
         {Node(430), Node(440)},  {Node(430), Node(435)},
         {Node(430), Node(447)},  {Node(431), Node(441)},
         {Node(431), Node(436)},  {Node(431), Node(448)},
         {Node(432), Node(442)},  {Node(432), Node(437)},
         {Node(432), Node(449)},  {Node(433), Node(443)},
         {Node(433), Node(438)},  {Node(433), Node(450)},
         {Node(434), Node(444)},  {Node(434), Node(439)},
         {Node(434), Node(451)},  {Node(435), Node(445)},
         {Node(435), Node(440)},  {Node(435), Node(452)},
         {Node(436), Node(446)},  {Node(436), Node(441)},
         {Node(436), Node(453)},  {Node(437), Node(447)},
         {Node(437), Node(442)},  {Node(437), Node(454)},
         {Node(438), Node(448)},  {Node(438), Node(443)},
         {Node(438), Node(455)},  {Node(439), Node(449)},
         {Node(439), Node(444)},  {Node(439), Node(456)},
         {Node(440), Node(450)},  {Node(440), Node(445)},
         {Node(440), Node(457)},  {Node(441), Node(451)},
         {Node(441), Node(446)},  {Node(441), Node(458)},
         {Node(442), Node(452)},  {Node(442), Node(447)},
         {Node(442), Node(459)},  {Node(443), Node(453)},
         {Node(443), Node(448)},  {Node(443), Node(460)},
         {Node(444), Node(454)},  {Node(444), Node(449)},
         {Node(444), Node(461)},  {Node(445), Node(455)},
         {Node(445), Node(450)},  {Node(445), Node(462)},
         {Node(446), Node(456)},  {Node(446), Node(451)},
         {Node(446), Node(463)},  {Node(447), Node(457)},
         {Node(447), Node(452)},  {Node(447), Node(464)},
         {Node(448), Node(458)},  {Node(448), Node(453)},
         {Node(448), Node(465)},  {Node(449), Node(459)},
         {Node(449), Node(454)},  {Node(449), Node(466)},
         {Node(450), Node(460)},  {Node(450), Node(455)},
         {Node(450), Node(467)},  {Node(451), Node(461)},
         {Node(451), Node(456)},  {Node(451), Node(468)},
         {Node(452), Node(462)},  {Node(452), Node(457)},
         {Node(452), Node(469)},  {Node(453), Node(463)},
         {Node(453), Node(458)},  {Node(453), Node(470)},
         {Node(454), Node(464)},  {Node(454), Node(459)},
         {Node(454), Node(471)},  {Node(455), Node(465)},
         {Node(455), Node(460)},  {Node(455), Node(472)},
         {Node(456), Node(466)},  {Node(456), Node(461)},
         {Node(456), Node(473)},  {Node(457), Node(467)},
         {Node(457), Node(462)},  {Node(457), Node(474)},
         {Node(458), Node(468)},  {Node(458), Node(463)},
         {Node(458), Node(475)},  {Node(459), Node(469)},
         {Node(459), Node(464)},  {Node(459), Node(476)},
         {Node(460), Node(470)},  {Node(460), Node(465)},
         {Node(460), Node(477)},  {Node(461), Node(471)},
         {Node(461), Node(466)},  {Node(461), Node(478)},
         {Node(462), Node(472)},  {Node(462), Node(467)},
         {Node(462), Node(479)},  {Node(463), Node(473)},
         {Node(463), Node(468)},  {Node(463), Node(480)},
         {Node(464), Node(474)},  {Node(464), Node(469)},
         {Node(464), Node(481)},  {Node(465), Node(475)},
         {Node(465), Node(470)},  {Node(465), Node(482)},
         {Node(466), Node(476)},  {Node(466), Node(471)},
         {Node(466), Node(483)},  {Node(467), Node(477)},
         {Node(467), Node(472)},  {Node(467), Node(484)},
         {Node(468), Node(478)},  {Node(468), Node(473)},
         {Node(468), Node(485)},  {Node(469), Node(479)},
         {Node(469), Node(474)},  {Node(469), Node(486)},
         {Node(470), Node(480)},  {Node(470), Node(475)},
         {Node(470), Node(487)},  {Node(471), Node(481)},
         {Node(471), Node(476)},  {Node(471), Node(488)},
         {Node(472), Node(482)},  {Node(472), Node(477)},
         {Node(472), Node(489)},  {Node(473), Node(483)},
         {Node(473), Node(478)},  {Node(473), Node(490)},
         {Node(474), Node(484)},  {Node(474), Node(479)},
         {Node(474), Node(491)},  {Node(475), Node(485)},
         {Node(475), Node(480)},  {Node(475), Node(492)},
         {Node(476), Node(486)},  {Node(476), Node(481)},
         {Node(476), Node(493)},  {Node(477), Node(487)},
         {Node(477), Node(482)},  {Node(477), Node(494)},
         {Node(478), Node(488)},  {Node(478), Node(483)},
         {Node(478), Node(495)},  {Node(479), Node(489)},
         {Node(479), Node(484)},  {Node(479), Node(496)},
         {Node(480), Node(490)},  {Node(480), Node(485)},
         {Node(480), Node(497)},  {Node(481), Node(491)},
         {Node(481), Node(486)},  {Node(481), Node(498)},
         {Node(482), Node(492)},  {Node(482), Node(487)},
         {Node(482), Node(499)},  {Node(483), Node(493)},
         {Node(483), Node(488)},  {Node(483), Node(500)},
         {Node(484), Node(494)},  {Node(484), Node(489)},
         {Node(484), Node(501)},  {Node(485), Node(495)},
         {Node(485), Node(490)},  {Node(485), Node(502)},
         {Node(486), Node(496)},  {Node(486), Node(491)},
         {Node(486), Node(503)},  {Node(487), Node(497)},
         {Node(487), Node(492)},  {Node(487), Node(504)},
         {Node(488), Node(498)},  {Node(488), Node(493)},
         {Node(488), Node(505)},  {Node(489), Node(499)},
         {Node(489), Node(494)},  {Node(489), Node(506)},
         {Node(490), Node(500)},  {Node(490), Node(495)},
         {Node(490), Node(507)},  {Node(491), Node(501)},
         {Node(491), Node(496)},  {Node(491), Node(508)},
         {Node(492), Node(502)},  {Node(492), Node(497)},
         {Node(492), Node(509)},  {Node(493), Node(503)},
         {Node(493), Node(498)},  {Node(493), Node(510)},
         {Node(494), Node(504)},  {Node(494), Node(499)},
         {Node(494), Node(511)},  {Node(495), Node(505)},
         {Node(495), Node(500)},  {Node(495), Node(512)},
         {Node(496), Node(506)},  {Node(496), Node(501)},
         {Node(496), Node(513)},  {Node(497), Node(507)},
         {Node(497), Node(502)},  {Node(497), Node(514)},
         {Node(498), Node(508)},  {Node(498), Node(503)},
         {Node(498), Node(515)},  {Node(499), Node(509)},
         {Node(499), Node(504)},  {Node(499), Node(516)},
         {Node(500), Node(510)},  {Node(500), Node(505)},
         {Node(500), Node(517)},  {Node(501), Node(511)},
         {Node(501), Node(506)},  {Node(501), Node(518)},
         {Node(502), Node(512)},  {Node(502), Node(507)},
         {Node(502), Node(519)},  {Node(503), Node(513)},
         {Node(503), Node(508)},  {Node(503), Node(520)},
         {Node(504), Node(514)},  {Node(504), Node(509)},
         {Node(504), Node(521)},  {Node(505), Node(515)},
         {Node(505), Node(510)},  {Node(505), Node(522)},
         {Node(506), Node(516)},  {Node(506), Node(511)},
         {Node(506), Node(523)},  {Node(507), Node(517)},
         {Node(507), Node(512)},  {Node(507), Node(524)},
         {Node(508), Node(518)},  {Node(508), Node(513)},
         {Node(508), Node(525)},  {Node(509), Node(519)},
         {Node(509), Node(514)},  {Node(509), Node(526)},
         {Node(510), Node(520)},  {Node(510), Node(515)},
         {Node(510), Node(527)},  {Node(511), Node(521)},
         {Node(511), Node(516)},  {Node(511), Node(528)},
         {Node(512), Node(522)},  {Node(512), Node(517)},
         {Node(512), Node(529)},  {Node(513), Node(523)},
         {Node(513), Node(518)},  {Node(513), Node(530)},
         {Node(514), Node(524)},  {Node(514), Node(519)},
         {Node(514), Node(531)},  {Node(515), Node(525)},
         {Node(515), Node(520)},  {Node(515), Node(532)},
         {Node(516), Node(526)},  {Node(516), Node(521)},
         {Node(516), Node(533)},  {Node(517), Node(527)},
         {Node(517), Node(522)},  {Node(517), Node(534)},
         {Node(518), Node(528)},  {Node(518), Node(523)},
         {Node(518), Node(535)},  {Node(519), Node(529)},
         {Node(519), Node(524)},  {Node(519), Node(536)},
         {Node(520), Node(530)},  {Node(520), Node(525)},
         {Node(520), Node(537)},  {Node(521), Node(531)},
         {Node(521), Node(526)},  {Node(521), Node(538)},
         {Node(522), Node(532)},  {Node(522), Node(527)},
         {Node(522), Node(539)},  {Node(523), Node(533)},
         {Node(523), Node(528)},  {Node(523), Node(540)},
         {Node(524), Node(534)},  {Node(524), Node(529)},
         {Node(524), Node(541)},  {Node(525), Node(535)},
         {Node(525), Node(530)},  {Node(525), Node(542)},
         {Node(526), Node(536)},  {Node(526), Node(531)},
         {Node(526), Node(543)},  {Node(527), Node(537)},
         {Node(527), Node(532)},  {Node(527), Node(544)},
         {Node(528), Node(538)},  {Node(528), Node(533)},
         {Node(528), Node(545)},  {Node(529), Node(539)},
         {Node(529), Node(534)},  {Node(529), Node(546)},
         {Node(530), Node(540)},  {Node(530), Node(535)},
         {Node(530), Node(547)},  {Node(531), Node(541)},
         {Node(531), Node(536)},  {Node(531), Node(548)},
         {Node(532), Node(542)},  {Node(532), Node(537)},
         {Node(532), Node(549)},  {Node(533), Node(543)},
         {Node(533), Node(538)},  {Node(533), Node(550)},
         {Node(534), Node(544)},  {Node(534), Node(539)},
         {Node(534), Node(551)},  {Node(535), Node(545)},
         {Node(535), Node(540)},  {Node(535), Node(552)},
         {Node(536), Node(546)},  {Node(536), Node(541)},
         {Node(536), Node(553)},  {Node(537), Node(547)},
         {Node(537), Node(542)},  {Node(537), Node(554)},
         {Node(538), Node(548)},  {Node(538), Node(543)},
         {Node(538), Node(555)},  {Node(539), Node(549)},
         {Node(539), Node(544)},  {Node(539), Node(556)},
         {Node(540), Node(550)},  {Node(540), Node(545)},
         {Node(540), Node(557)},  {Node(541), Node(551)},
         {Node(541), Node(546)},  {Node(541), Node(558)},
         {Node(542), Node(552)},  {Node(542), Node(547)},
         {Node(542), Node(559)},  {Node(543), Node(553)},
         {Node(543), Node(548)},  {Node(543), Node(560)},
         {Node(544), Node(554)},  {Node(544), Node(549)},
         {Node(544), Node(561)},  {Node(545), Node(555)},
         {Node(545), Node(550)},  {Node(545), Node(562)},
         {Node(546), Node(556)},  {Node(546), Node(551)},
         {Node(546), Node(563)},  {Node(547), Node(557)},
         {Node(547), Node(552)},  {Node(547), Node(564)},
         {Node(548), Node(558)},  {Node(548), Node(553)},
         {Node(548), Node(565)},  {Node(549), Node(559)},
         {Node(549), Node(554)},  {Node(549), Node(566)},
         {Node(550), Node(560)},  {Node(550), Node(555)},
         {Node(550), Node(567)},  {Node(551), Node(561)},
         {Node(551), Node(556)},  {Node(551), Node(568)},
         {Node(552), Node(562)},  {Node(552), Node(557)},
         {Node(552), Node(569)},  {Node(553), Node(563)},
         {Node(553), Node(558)},  {Node(553), Node(570)},
         {Node(554), Node(564)},  {Node(554), Node(559)},
         {Node(554), Node(571)},  {Node(555), Node(565)},
         {Node(555), Node(560)},  {Node(555), Node(572)},
         {Node(556), Node(566)},  {Node(556), Node(561)},
         {Node(556), Node(573)},  {Node(557), Node(567)},
         {Node(557), Node(562)},  {Node(557), Node(574)},
         {Node(558), Node(568)},  {Node(558), Node(563)},
         {Node(558), Node(575)},  {Node(559), Node(569)},
         {Node(559), Node(564)},  {Node(559), Node(576)},
         {Node(560), Node(570)},  {Node(560), Node(565)},
         {Node(560), Node(577)},  {Node(561), Node(571)},
         {Node(561), Node(566)},  {Node(561), Node(578)},
         {Node(562), Node(572)},  {Node(562), Node(567)},
         {Node(562), Node(579)},  {Node(563), Node(573)},
         {Node(563), Node(568)},  {Node(563), Node(580)},
         {Node(564), Node(574)},  {Node(564), Node(569)},
         {Node(564), Node(581)},  {Node(565), Node(575)},
         {Node(565), Node(570)},  {Node(565), Node(582)},
         {Node(566), Node(576)},  {Node(566), Node(571)},
         {Node(566), Node(583)},  {Node(567), Node(577)},
         {Node(567), Node(572)},  {Node(567), Node(584)},
         {Node(568), Node(578)},  {Node(568), Node(573)},
         {Node(568), Node(585)},  {Node(569), Node(579)},
         {Node(569), Node(574)},  {Node(569), Node(586)},
         {Node(570), Node(580)},  {Node(570), Node(575)},
         {Node(570), Node(587)},  {Node(571), Node(581)},
         {Node(571), Node(576)},  {Node(571), Node(588)},
         {Node(572), Node(582)},  {Node(572), Node(577)},
         {Node(572), Node(589)},  {Node(573), Node(583)},
         {Node(573), Node(578)},  {Node(573), Node(590)},
         {Node(574), Node(584)},  {Node(574), Node(579)},
         {Node(574), Node(591)},  {Node(575), Node(585)},
         {Node(575), Node(580)},  {Node(575), Node(592)},
         {Node(576), Node(586)},  {Node(576), Node(581)},
         {Node(576), Node(593)},  {Node(577), Node(587)},
         {Node(577), Node(582)},  {Node(577), Node(594)},
         {Node(578), Node(588)},  {Node(578), Node(583)},
         {Node(578), Node(595)},  {Node(579), Node(589)},
         {Node(579), Node(584)},  {Node(579), Node(596)},
         {Node(580), Node(590)},  {Node(580), Node(585)},
         {Node(580), Node(597)},  {Node(581), Node(591)},
         {Node(581), Node(586)},  {Node(581), Node(598)},
         {Node(582), Node(592)},  {Node(582), Node(587)},
         {Node(582), Node(599)},  {Node(583), Node(593)},
         {Node(583), Node(588)},  {Node(583), Node(600)},
         {Node(584), Node(594)},  {Node(584), Node(589)},
         {Node(584), Node(601)},  {Node(585), Node(595)},
         {Node(585), Node(590)},  {Node(585), Node(602)},
         {Node(586), Node(596)},  {Node(586), Node(591)},
         {Node(586), Node(603)},  {Node(587), Node(597)},
         {Node(587), Node(592)},  {Node(587), Node(604)},
         {Node(588), Node(598)},  {Node(588), Node(593)},
         {Node(588), Node(605)},  {Node(589), Node(599)},
         {Node(589), Node(594)},  {Node(589), Node(606)},
         {Node(590), Node(600)},  {Node(590), Node(595)},
         {Node(590), Node(607)},  {Node(591), Node(601)},
         {Node(591), Node(596)},  {Node(591), Node(608)},
         {Node(592), Node(602)},  {Node(592), Node(597)},
         {Node(592), Node(609)},  {Node(593), Node(603)},
         {Node(593), Node(598)},  {Node(593), Node(610)},
         {Node(594), Node(604)},  {Node(594), Node(599)},
         {Node(594), Node(611)},  {Node(595), Node(605)},
         {Node(595), Node(600)},  {Node(595), Node(612)},
         {Node(596), Node(606)},  {Node(596), Node(601)},
         {Node(596), Node(613)},  {Node(597), Node(607)},
         {Node(597), Node(602)},  {Node(597), Node(614)},
         {Node(598), Node(608)},  {Node(598), Node(603)},
         {Node(598), Node(615)},  {Node(599), Node(609)},
         {Node(599), Node(604)},  {Node(599), Node(616)},
         {Node(600), Node(610)},  {Node(600), Node(605)},
         {Node(600), Node(617)},  {Node(601), Node(611)},
         {Node(601), Node(606)},  {Node(601), Node(618)},
         {Node(602), Node(612)},  {Node(602), Node(607)},
         {Node(602), Node(619)},  {Node(603), Node(613)},
         {Node(603), Node(608)},  {Node(603), Node(620)},
         {Node(604), Node(614)},  {Node(604), Node(609)},
         {Node(604), Node(621)},  {Node(605), Node(615)},
         {Node(605), Node(610)},  {Node(605), Node(622)},
         {Node(606), Node(616)},  {Node(606), Node(611)},
         {Node(606), Node(623)},  {Node(607), Node(617)},
         {Node(607), Node(612)},  {Node(607), Node(624)},
         {Node(608), Node(618)},  {Node(608), Node(613)},
         {Node(608), Node(625)},  {Node(609), Node(619)},
         {Node(609), Node(614)},  {Node(609), Node(626)},
         {Node(610), Node(620)},  {Node(610), Node(615)},
         {Node(610), Node(627)},  {Node(611), Node(621)},
         {Node(611), Node(616)},  {Node(611), Node(628)},
         {Node(612), Node(622)},  {Node(612), Node(617)},
         {Node(612), Node(629)},  {Node(613), Node(623)},
         {Node(613), Node(618)},  {Node(613), Node(630)},
         {Node(614), Node(624)},  {Node(614), Node(619)},
         {Node(614), Node(631)},  {Node(615), Node(625)},
         {Node(615), Node(620)},  {Node(615), Node(632)},
         {Node(616), Node(626)},  {Node(616), Node(621)},
         {Node(616), Node(633)},  {Node(617), Node(627)},
         {Node(617), Node(622)},  {Node(617), Node(634)},
         {Node(618), Node(628)},  {Node(618), Node(623)},
         {Node(618), Node(635)},  {Node(619), Node(629)},
         {Node(619), Node(624)},  {Node(619), Node(636)},
         {Node(620), Node(630)},  {Node(620), Node(625)},
         {Node(620), Node(637)},  {Node(621), Node(631)},
         {Node(621), Node(626)},  {Node(621), Node(638)},
         {Node(622), Node(632)},  {Node(622), Node(627)},
         {Node(622), Node(639)},  {Node(623), Node(633)},
         {Node(623), Node(628)},  {Node(623), Node(640)},
         {Node(624), Node(634)},  {Node(624), Node(629)},
         {Node(624), Node(641)},  {Node(625), Node(635)},
         {Node(625), Node(630)},  {Node(625), Node(642)},
         {Node(626), Node(636)},  {Node(626), Node(631)},
         {Node(626), Node(643)},  {Node(627), Node(637)},
         {Node(627), Node(632)},  {Node(627), Node(644)},
         {Node(628), Node(638)},  {Node(628), Node(633)},
         {Node(628), Node(645)},  {Node(629), Node(639)},
         {Node(629), Node(634)},  {Node(629), Node(646)},
         {Node(630), Node(640)},  {Node(630), Node(635)},
         {Node(630), Node(647)},  {Node(631), Node(641)},
         {Node(631), Node(636)},  {Node(631), Node(648)},
         {Node(632), Node(642)},  {Node(632), Node(637)},
         {Node(632), Node(649)},  {Node(633), Node(643)},
         {Node(633), Node(638)},  {Node(633), Node(650)},
         {Node(634), Node(644)},  {Node(634), Node(639)},
         {Node(634), Node(651)},  {Node(635), Node(645)},
         {Node(635), Node(640)},  {Node(635), Node(652)},
         {Node(636), Node(646)},  {Node(636), Node(641)},
         {Node(636), Node(653)},  {Node(637), Node(647)},
         {Node(637), Node(642)},  {Node(637), Node(654)},
         {Node(638), Node(648)},  {Node(638), Node(643)},
         {Node(638), Node(655)},  {Node(639), Node(649)},
         {Node(639), Node(644)},  {Node(639), Node(656)},
         {Node(640), Node(650)},  {Node(640), Node(645)},
         {Node(640), Node(657)},  {Node(641), Node(651)},
         {Node(641), Node(646)},  {Node(641), Node(658)},
         {Node(642), Node(652)},  {Node(642), Node(647)},
         {Node(642), Node(659)},  {Node(643), Node(653)},
         {Node(643), Node(648)},  {Node(643), Node(660)},
         {Node(644), Node(654)},  {Node(644), Node(649)},
         {Node(644), Node(661)},  {Node(645), Node(655)},
         {Node(645), Node(650)},  {Node(645), Node(662)},
         {Node(646), Node(656)},  {Node(646), Node(651)},
         {Node(646), Node(663)},  {Node(647), Node(657)},
         {Node(647), Node(652)},  {Node(647), Node(664)},
         {Node(648), Node(658)},  {Node(648), Node(653)},
         {Node(648), Node(665)},  {Node(649), Node(659)},
         {Node(649), Node(654)},  {Node(649), Node(666)},
         {Node(650), Node(660)},  {Node(650), Node(655)},
         {Node(650), Node(667)},  {Node(651), Node(661)},
         {Node(651), Node(656)},  {Node(651), Node(668)},
         {Node(652), Node(662)},  {Node(652), Node(657)},
         {Node(652), Node(669)},  {Node(653), Node(663)},
         {Node(653), Node(658)},  {Node(653), Node(670)},
         {Node(654), Node(664)},  {Node(654), Node(659)},
         {Node(654), Node(671)},  {Node(655), Node(665)},
         {Node(655), Node(660)},  {Node(655), Node(672)},
         {Node(656), Node(666)},  {Node(656), Node(661)},
         {Node(656), Node(673)},  {Node(657), Node(667)},
         {Node(657), Node(662)},  {Node(657), Node(674)},
         {Node(658), Node(668)},  {Node(658), Node(663)},
         {Node(658), Node(675)},  {Node(659), Node(669)},
         {Node(659), Node(664)},  {Node(659), Node(676)},
         {Node(660), Node(670)},  {Node(660), Node(665)},
         {Node(660), Node(677)},  {Node(661), Node(671)},
         {Node(661), Node(666)},  {Node(661), Node(678)},
         {Node(662), Node(672)},  {Node(662), Node(667)},
         {Node(662), Node(679)},  {Node(663), Node(673)},
         {Node(663), Node(668)},  {Node(663), Node(680)},
         {Node(664), Node(674)},  {Node(664), Node(669)},
         {Node(664), Node(681)},  {Node(665), Node(675)},
         {Node(665), Node(670)},  {Node(665), Node(682)},
         {Node(666), Node(676)},  {Node(666), Node(671)},
         {Node(666), Node(683)},  {Node(667), Node(677)},
         {Node(667), Node(672)},  {Node(667), Node(684)},
         {Node(668), Node(678)},  {Node(668), Node(673)},
         {Node(668), Node(685)},  {Node(669), Node(679)},
         {Node(669), Node(674)},  {Node(669), Node(686)},
         {Node(670), Node(680)},  {Node(670), Node(675)},
         {Node(670), Node(687)},  {Node(671), Node(681)},
         {Node(671), Node(676)},  {Node(671), Node(688)},
         {Node(672), Node(682)},  {Node(672), Node(677)},
         {Node(672), Node(689)},  {Node(673), Node(683)},
         {Node(673), Node(678)},  {Node(673), Node(690)},
         {Node(674), Node(684)},  {Node(674), Node(679)},
         {Node(674), Node(691)},  {Node(675), Node(685)},
         {Node(675), Node(680)},  {Node(675), Node(692)},
         {Node(676), Node(686)},  {Node(676), Node(681)},
         {Node(676), Node(693)},  {Node(677), Node(687)},
         {Node(677), Node(682)},  {Node(677), Node(694)},
         {Node(678), Node(688)},  {Node(678), Node(683)},
         {Node(678), Node(695)},  {Node(679), Node(689)},
         {Node(679), Node(684)},  {Node(679), Node(696)},
         {Node(680), Node(690)},  {Node(680), Node(685)},
         {Node(680), Node(697)},  {Node(681), Node(691)},
         {Node(681), Node(686)},  {Node(681), Node(698)},
         {Node(682), Node(692)},  {Node(682), Node(687)},
         {Node(682), Node(699)},  {Node(683), Node(693)},
         {Node(683), Node(688)},  {Node(683), Node(700)},
         {Node(684), Node(694)},  {Node(684), Node(689)},
         {Node(684), Node(701)},  {Node(685), Node(695)},
         {Node(685), Node(690)},  {Node(685), Node(702)},
         {Node(686), Node(696)},  {Node(686), Node(691)},
         {Node(686), Node(703)},  {Node(687), Node(697)},
         {Node(687), Node(692)},  {Node(687), Node(704)},
         {Node(688), Node(698)},  {Node(688), Node(693)},
         {Node(688), Node(705)},  {Node(689), Node(699)},
         {Node(689), Node(694)},  {Node(689), Node(706)},
         {Node(690), Node(700)},  {Node(690), Node(695)},
         {Node(690), Node(707)},  {Node(691), Node(701)},
         {Node(691), Node(696)},  {Node(691), Node(708)},
         {Node(692), Node(702)},  {Node(692), Node(697)},
         {Node(692), Node(709)},  {Node(693), Node(703)},
         {Node(693), Node(698)},  {Node(693), Node(710)},
         {Node(694), Node(704)},  {Node(694), Node(699)},
         {Node(694), Node(711)},  {Node(695), Node(705)},
         {Node(695), Node(700)},  {Node(695), Node(712)},
         {Node(696), Node(706)},  {Node(696), Node(701)},
         {Node(696), Node(713)},  {Node(697), Node(707)},
         {Node(697), Node(702)},  {Node(697), Node(714)},
         {Node(698), Node(708)},  {Node(698), Node(703)},
         {Node(698), Node(715)},  {Node(699), Node(709)},
         {Node(699), Node(704)},  {Node(699), Node(716)},
         {Node(700), Node(710)},  {Node(700), Node(705)},
         {Node(700), Node(717)},  {Node(701), Node(711)},
         {Node(701), Node(706)},  {Node(701), Node(718)},
         {Node(702), Node(712)},  {Node(702), Node(707)},
         {Node(702), Node(719)},  {Node(703), Node(713)},
         {Node(703), Node(708)},  {Node(703), Node(720)},
         {Node(704), Node(714)},  {Node(704), Node(709)},
         {Node(704), Node(721)},  {Node(705), Node(715)},
         {Node(705), Node(710)},  {Node(705), Node(722)},
         {Node(706), Node(716)},  {Node(706), Node(711)},
         {Node(706), Node(723)},  {Node(707), Node(717)},
         {Node(707), Node(712)},  {Node(707), Node(724)},
         {Node(708), Node(718)},  {Node(708), Node(713)},
         {Node(708), Node(725)},  {Node(709), Node(719)},
         {Node(709), Node(714)},  {Node(709), Node(726)},
         {Node(710), Node(720)},  {Node(710), Node(715)},
         {Node(710), Node(727)},  {Node(711), Node(721)},
         {Node(711), Node(716)},  {Node(711), Node(728)},
         {Node(712), Node(722)},  {Node(712), Node(717)},
         {Node(712), Node(729)},  {Node(713), Node(723)},
         {Node(713), Node(718)},  {Node(713), Node(730)},
         {Node(714), Node(724)},  {Node(714), Node(719)},
         {Node(714), Node(731)},  {Node(715), Node(725)},
         {Node(715), Node(720)},  {Node(715), Node(732)},
         {Node(716), Node(726)},  {Node(716), Node(721)},
         {Node(716), Node(733)},  {Node(717), Node(727)},
         {Node(717), Node(722)},  {Node(717), Node(734)},
         {Node(718), Node(728)},  {Node(718), Node(723)},
         {Node(718), Node(735)},  {Node(719), Node(729)},
         {Node(719), Node(724)},  {Node(719), Node(736)},
         {Node(720), Node(730)},  {Node(720), Node(725)},
         {Node(720), Node(737)},  {Node(721), Node(731)},
         {Node(721), Node(726)},  {Node(721), Node(738)},
         {Node(722), Node(732)},  {Node(722), Node(727)},
         {Node(722), Node(739)},  {Node(723), Node(733)},
         {Node(723), Node(728)},  {Node(723), Node(740)},
         {Node(724), Node(734)},  {Node(724), Node(729)},
         {Node(724), Node(741)},  {Node(725), Node(735)},
         {Node(725), Node(730)},  {Node(725), Node(742)},
         {Node(726), Node(736)},  {Node(726), Node(731)},
         {Node(726), Node(743)},  {Node(727), Node(737)},
         {Node(727), Node(732)},  {Node(727), Node(744)},
         {Node(728), Node(738)},  {Node(728), Node(733)},
         {Node(728), Node(745)},  {Node(729), Node(739)},
         {Node(729), Node(734)},  {Node(729), Node(746)},
         {Node(730), Node(740)},  {Node(730), Node(735)},
         {Node(730), Node(747)},  {Node(731), Node(741)},
         {Node(731), Node(736)},  {Node(731), Node(748)},
         {Node(732), Node(742)},  {Node(732), Node(737)},
         {Node(732), Node(749)},  {Node(733), Node(743)},
         {Node(733), Node(738)},  {Node(733), Node(750)},
         {Node(734), Node(744)},  {Node(734), Node(739)},
         {Node(734), Node(751)},  {Node(735), Node(745)},
         {Node(735), Node(740)},  {Node(735), Node(752)},
         {Node(736), Node(746)},  {Node(736), Node(741)},
         {Node(736), Node(753)},  {Node(737), Node(747)},
         {Node(737), Node(742)},  {Node(737), Node(754)},
         {Node(738), Node(748)},  {Node(738), Node(743)},
         {Node(738), Node(755)},  {Node(739), Node(749)},
         {Node(739), Node(744)},  {Node(739), Node(756)},
         {Node(740), Node(750)},  {Node(740), Node(745)},
         {Node(740), Node(757)},  {Node(741), Node(751)},
         {Node(741), Node(746)},  {Node(741), Node(758)},
         {Node(742), Node(752)},  {Node(742), Node(747)},
         {Node(742), Node(759)},  {Node(743), Node(753)},
         {Node(743), Node(748)},  {Node(743), Node(760)},
         {Node(744), Node(754)},  {Node(744), Node(749)},
         {Node(744), Node(761)},  {Node(745), Node(755)},
         {Node(745), Node(750)},  {Node(745), Node(762)},
         {Node(746), Node(756)},  {Node(746), Node(751)},
         {Node(746), Node(763)},  {Node(747), Node(757)},
         {Node(747), Node(752)},  {Node(747), Node(764)},
         {Node(748), Node(758)},  {Node(748), Node(753)},
         {Node(748), Node(765)},  {Node(749), Node(759)},
         {Node(749), Node(754)},  {Node(749), Node(766)},
         {Node(750), Node(760)},  {Node(750), Node(755)},
         {Node(750), Node(767)},  {Node(751), Node(761)},
         {Node(751), Node(756)},  {Node(751), Node(768)},
         {Node(752), Node(762)},  {Node(752), Node(757)},
         {Node(752), Node(769)},  {Node(753), Node(763)},
         {Node(753), Node(758)},  {Node(753), Node(770)},
         {Node(754), Node(764)},  {Node(754), Node(759)},
         {Node(754), Node(771)},  {Node(755), Node(765)},
         {Node(755), Node(760)},  {Node(755), Node(772)},
         {Node(756), Node(766)},  {Node(756), Node(761)},
         {Node(756), Node(773)},  {Node(757), Node(767)},
         {Node(757), Node(762)},  {Node(757), Node(774)},
         {Node(758), Node(768)},  {Node(758), Node(763)},
         {Node(758), Node(775)},  {Node(759), Node(769)},
         {Node(759), Node(764)},  {Node(759), Node(776)},
         {Node(760), Node(770)},  {Node(760), Node(765)},
         {Node(760), Node(777)},  {Node(761), Node(771)},
         {Node(761), Node(766)},  {Node(761), Node(778)},
         {Node(762), Node(772)},  {Node(762), Node(767)},
         {Node(762), Node(779)},  {Node(763), Node(773)},
         {Node(763), Node(768)},  {Node(763), Node(780)},
         {Node(764), Node(774)},  {Node(764), Node(769)},
         {Node(764), Node(781)},  {Node(765), Node(775)},
         {Node(765), Node(770)},  {Node(765), Node(782)},
         {Node(766), Node(776)},  {Node(766), Node(771)},
         {Node(766), Node(783)},  {Node(767), Node(777)},
         {Node(767), Node(772)},  {Node(767), Node(784)},
         {Node(768), Node(778)},  {Node(768), Node(773)},
         {Node(768), Node(785)},  {Node(769), Node(779)},
         {Node(769), Node(774)},  {Node(769), Node(786)},
         {Node(770), Node(780)},  {Node(770), Node(775)},
         {Node(770), Node(787)},  {Node(771), Node(781)},
         {Node(771), Node(776)},  {Node(771), Node(788)},
         {Node(772), Node(782)},  {Node(772), Node(777)},
         {Node(772), Node(789)},  {Node(773), Node(783)},
         {Node(773), Node(778)},  {Node(773), Node(790)},
         {Node(774), Node(784)},  {Node(774), Node(779)},
         {Node(774), Node(791)},  {Node(775), Node(785)},
         {Node(775), Node(780)},  {Node(775), Node(792)},
         {Node(776), Node(786)},  {Node(776), Node(781)},
         {Node(776), Node(793)},  {Node(777), Node(787)},
         {Node(777), Node(782)},  {Node(777), Node(794)},
         {Node(778), Node(788)},  {Node(778), Node(783)},
         {Node(778), Node(795)},  {Node(779), Node(789)},
         {Node(779), Node(784)},  {Node(779), Node(796)},
         {Node(780), Node(790)},  {Node(780), Node(785)},
         {Node(780), Node(797)},  {Node(781), Node(791)},
         {Node(781), Node(786)},  {Node(781), Node(798)},
         {Node(782), Node(792)},  {Node(782), Node(787)},
         {Node(782), Node(799)},  {Node(783), Node(793)},
         {Node(783), Node(788)},  {Node(783), Node(800)},
         {Node(784), Node(794)},  {Node(784), Node(789)},
         {Node(784), Node(801)},  {Node(785), Node(795)},
         {Node(785), Node(790)},  {Node(785), Node(802)},
         {Node(786), Node(796)},  {Node(786), Node(791)},
         {Node(786), Node(803)},  {Node(787), Node(797)},
         {Node(787), Node(792)},  {Node(787), Node(804)},
         {Node(788), Node(798)},  {Node(788), Node(793)},
         {Node(788), Node(805)},  {Node(789), Node(799)},
         {Node(789), Node(794)},  {Node(789), Node(806)},
         {Node(790), Node(800)},  {Node(790), Node(795)},
         {Node(790), Node(807)},  {Node(791), Node(801)},
         {Node(791), Node(796)},  {Node(791), Node(808)},
         {Node(792), Node(802)},  {Node(792), Node(797)},
         {Node(792), Node(809)},  {Node(793), Node(803)},
         {Node(793), Node(798)},  {Node(793), Node(810)},
         {Node(794), Node(804)},  {Node(794), Node(799)},
         {Node(794), Node(811)},  {Node(795), Node(805)},
         {Node(795), Node(800)},  {Node(795), Node(812)},
         {Node(796), Node(806)},  {Node(796), Node(801)},
         {Node(796), Node(813)},  {Node(797), Node(807)},
         {Node(797), Node(802)},  {Node(797), Node(814)},
         {Node(798), Node(808)},  {Node(798), Node(803)},
         {Node(798), Node(815)},  {Node(799), Node(809)},
         {Node(799), Node(804)},  {Node(799), Node(816)},
         {Node(800), Node(810)},  {Node(800), Node(805)},
         {Node(800), Node(817)},  {Node(801), Node(811)},
         {Node(801), Node(806)},  {Node(801), Node(818)},
         {Node(802), Node(812)},  {Node(802), Node(807)},
         {Node(802), Node(819)},  {Node(803), Node(813)},
         {Node(803), Node(808)},  {Node(803), Node(820)},
         {Node(804), Node(814)},  {Node(804), Node(809)},
         {Node(804), Node(821)},  {Node(805), Node(815)},
         {Node(805), Node(810)},  {Node(805), Node(822)},
         {Node(806), Node(816)},  {Node(806), Node(811)},
         {Node(806), Node(823)},  {Node(807), Node(817)},
         {Node(807), Node(812)},  {Node(807), Node(824)},
         {Node(808), Node(818)},  {Node(808), Node(813)},
         {Node(808), Node(825)},  {Node(809), Node(819)},
         {Node(809), Node(814)},  {Node(809), Node(826)},
         {Node(810), Node(820)},  {Node(810), Node(815)},
         {Node(810), Node(827)},  {Node(811), Node(821)},
         {Node(811), Node(816)},  {Node(811), Node(828)},
         {Node(812), Node(822)},  {Node(812), Node(817)},
         {Node(812), Node(829)},  {Node(813), Node(823)},
         {Node(813), Node(818)},  {Node(813), Node(830)},
         {Node(814), Node(824)},  {Node(814), Node(819)},
         {Node(814), Node(831)},  {Node(815), Node(825)},
         {Node(815), Node(820)},  {Node(815), Node(832)},
         {Node(816), Node(826)},  {Node(816), Node(821)},
         {Node(816), Node(833)},  {Node(817), Node(827)},
         {Node(817), Node(822)},  {Node(817), Node(834)},
         {Node(818), Node(828)},  {Node(818), Node(823)},
         {Node(818), Node(835)},  {Node(819), Node(829)},
         {Node(819), Node(824)},  {Node(819), Node(836)},
         {Node(820), Node(830)},  {Node(820), Node(825)},
         {Node(820), Node(837)},  {Node(821), Node(831)},
         {Node(821), Node(826)},  {Node(821), Node(838)},
         {Node(822), Node(832)},  {Node(822), Node(827)},
         {Node(822), Node(839)},  {Node(823), Node(833)},
         {Node(823), Node(828)},  {Node(823), Node(840)},
         {Node(824), Node(834)},  {Node(824), Node(829)},
         {Node(824), Node(841)},  {Node(825), Node(835)},
         {Node(825), Node(830)},  {Node(825), Node(842)},
         {Node(826), Node(836)},  {Node(826), Node(831)},
         {Node(826), Node(843)},  {Node(827), Node(837)},
         {Node(827), Node(832)},  {Node(827), Node(844)},
         {Node(828), Node(838)},  {Node(828), Node(833)},
         {Node(828), Node(845)},  {Node(829), Node(839)},
         {Node(829), Node(834)},  {Node(829), Node(846)},
         {Node(830), Node(840)},  {Node(830), Node(835)},
         {Node(830), Node(847)},  {Node(831), Node(841)},
         {Node(831), Node(836)},  {Node(831), Node(848)},
         {Node(832), Node(842)},  {Node(832), Node(837)},
         {Node(832), Node(849)},  {Node(833), Node(843)},
         {Node(833), Node(838)},  {Node(833), Node(850)},
         {Node(834), Node(844)},  {Node(834), Node(839)},
         {Node(834), Node(851)},  {Node(835), Node(845)},
         {Node(835), Node(840)},  {Node(835), Node(852)},
         {Node(836), Node(846)},  {Node(836), Node(841)},
         {Node(836), Node(853)},  {Node(837), Node(847)},
         {Node(837), Node(842)},  {Node(837), Node(854)},
         {Node(838), Node(848)},  {Node(838), Node(843)},
         {Node(838), Node(855)},  {Node(839), Node(849)},
         {Node(839), Node(844)},  {Node(839), Node(856)},
         {Node(840), Node(850)},  {Node(840), Node(845)},
         {Node(840), Node(857)},  {Node(841), Node(851)},
         {Node(841), Node(846)},  {Node(841), Node(858)},
         {Node(842), Node(852)},  {Node(842), Node(847)},
         {Node(842), Node(859)},  {Node(843), Node(853)},
         {Node(843), Node(848)},  {Node(843), Node(860)},
         {Node(844), Node(854)},  {Node(844), Node(849)},
         {Node(844), Node(861)},  {Node(845), Node(855)},
         {Node(845), Node(850)},  {Node(845), Node(862)},
         {Node(846), Node(856)},  {Node(846), Node(851)},
         {Node(846), Node(863)},  {Node(847), Node(857)},
         {Node(847), Node(852)},  {Node(847), Node(864)},
         {Node(848), Node(858)},  {Node(848), Node(853)},
         {Node(848), Node(865)},  {Node(849), Node(859)},
         {Node(849), Node(854)},  {Node(849), Node(866)},
         {Node(850), Node(860)},  {Node(850), Node(855)},
         {Node(850), Node(867)},  {Node(851), Node(861)},
         {Node(851), Node(856)},  {Node(851), Node(868)},
         {Node(852), Node(862)},  {Node(852), Node(857)},
         {Node(852), Node(869)},  {Node(853), Node(863)},
         {Node(853), Node(858)},  {Node(853), Node(870)},
         {Node(854), Node(864)},  {Node(854), Node(859)},
         {Node(854), Node(871)},  {Node(855), Node(865)},
         {Node(855), Node(860)},  {Node(855), Node(872)},
         {Node(856), Node(866)},  {Node(856), Node(861)},
         {Node(856), Node(873)},  {Node(857), Node(867)},
         {Node(857), Node(862)},  {Node(857), Node(874)},
         {Node(858), Node(868)},  {Node(858), Node(863)},
         {Node(858), Node(875)},  {Node(859), Node(869)},
         {Node(859), Node(864)},  {Node(859), Node(876)},
         {Node(860), Node(870)},  {Node(860), Node(865)},
         {Node(860), Node(877)},  {Node(861), Node(871)},
         {Node(861), Node(866)},  {Node(861), Node(878)},
         {Node(862), Node(872)},  {Node(862), Node(867)},
         {Node(862), Node(879)},  {Node(863), Node(873)},
         {Node(863), Node(868)},  {Node(863), Node(880)},
         {Node(864), Node(874)},  {Node(864), Node(869)},
         {Node(864), Node(881)},  {Node(865), Node(875)},
         {Node(865), Node(870)},  {Node(865), Node(882)},
         {Node(866), Node(876)},  {Node(866), Node(871)},
         {Node(866), Node(883)},  {Node(867), Node(877)},
         {Node(867), Node(872)},  {Node(867), Node(884)},
         {Node(868), Node(878)},  {Node(868), Node(873)},
         {Node(868), Node(885)},  {Node(869), Node(879)},
         {Node(869), Node(874)},  {Node(869), Node(886)},
         {Node(870), Node(880)},  {Node(870), Node(875)},
         {Node(870), Node(887)},  {Node(871), Node(881)},
         {Node(871), Node(876)},  {Node(871), Node(888)},
         {Node(872), Node(882)},  {Node(872), Node(877)},
         {Node(872), Node(889)},  {Node(873), Node(883)},
         {Node(873), Node(878)},  {Node(873), Node(890)},
         {Node(874), Node(884)},  {Node(874), Node(879)},
         {Node(874), Node(891)},  {Node(875), Node(885)},
         {Node(875), Node(880)},  {Node(875), Node(892)},
         {Node(876), Node(886)},  {Node(876), Node(881)},
         {Node(876), Node(893)},  {Node(877), Node(887)},
         {Node(877), Node(882)},  {Node(877), Node(894)},
         {Node(878), Node(888)},  {Node(878), Node(883)},
         {Node(878), Node(895)},  {Node(879), Node(889)},
         {Node(879), Node(884)},  {Node(879), Node(896)},
         {Node(880), Node(890)},  {Node(880), Node(885)},
         {Node(880), Node(897)},  {Node(881), Node(891)},
         {Node(881), Node(886)},  {Node(881), Node(898)},
         {Node(882), Node(892)},  {Node(882), Node(887)},
         {Node(882), Node(899)},  {Node(883), Node(893)},
         {Node(883), Node(888)},  {Node(883), Node(900)},
         {Node(884), Node(894)},  {Node(884), Node(889)},
         {Node(884), Node(901)},  {Node(885), Node(895)},
         {Node(885), Node(890)},  {Node(885), Node(902)},
         {Node(886), Node(896)},  {Node(886), Node(891)},
         {Node(886), Node(903)},  {Node(887), Node(897)},
         {Node(887), Node(892)},  {Node(887), Node(904)},
         {Node(888), Node(898)},  {Node(888), Node(893)},
         {Node(888), Node(905)},  {Node(889), Node(899)},
         {Node(889), Node(894)},  {Node(889), Node(906)},
         {Node(890), Node(900)},  {Node(890), Node(895)},
         {Node(890), Node(907)},  {Node(891), Node(901)},
         {Node(891), Node(896)},  {Node(891), Node(908)},
         {Node(892), Node(902)},  {Node(892), Node(897)},
         {Node(892), Node(909)},  {Node(893), Node(903)},
         {Node(893), Node(898)},  {Node(893), Node(910)},
         {Node(894), Node(904)},  {Node(894), Node(899)},
         {Node(894), Node(911)},  {Node(895), Node(905)},
         {Node(895), Node(900)},  {Node(895), Node(912)},
         {Node(896), Node(906)},  {Node(896), Node(901)},
         {Node(896), Node(913)},  {Node(897), Node(907)},
         {Node(897), Node(902)},  {Node(897), Node(914)},
         {Node(898), Node(908)},  {Node(898), Node(903)},
         {Node(898), Node(915)},  {Node(899), Node(909)},
         {Node(899), Node(904)},  {Node(899), Node(916)},
         {Node(900), Node(910)},  {Node(900), Node(905)},
         {Node(900), Node(917)},  {Node(901), Node(911)},
         {Node(901), Node(906)},  {Node(901), Node(918)},
         {Node(902), Node(912)},  {Node(902), Node(907)},
         {Node(902), Node(919)},  {Node(903), Node(913)},
         {Node(903), Node(908)},  {Node(903), Node(920)},
         {Node(904), Node(914)},  {Node(904), Node(909)},
         {Node(904), Node(921)},  {Node(905), Node(915)},
         {Node(905), Node(910)},  {Node(905), Node(922)},
         {Node(906), Node(916)},  {Node(906), Node(911)},
         {Node(906), Node(923)},  {Node(907), Node(917)},
         {Node(907), Node(912)},  {Node(907), Node(924)},
         {Node(908), Node(918)},  {Node(908), Node(913)},
         {Node(908), Node(925)},  {Node(909), Node(919)},
         {Node(909), Node(914)},  {Node(909), Node(926)},
         {Node(910), Node(920)},  {Node(910), Node(915)},
         {Node(910), Node(927)},  {Node(911), Node(921)},
         {Node(911), Node(916)},  {Node(911), Node(928)},
         {Node(912), Node(922)},  {Node(912), Node(917)},
         {Node(912), Node(929)},  {Node(913), Node(923)},
         {Node(913), Node(918)},  {Node(913), Node(930)},
         {Node(914), Node(924)},  {Node(914), Node(919)},
         {Node(914), Node(931)},  {Node(915), Node(925)},
         {Node(915), Node(920)},  {Node(915), Node(932)},
         {Node(916), Node(926)},  {Node(916), Node(921)},
         {Node(916), Node(933)},  {Node(917), Node(927)},
         {Node(917), Node(922)},  {Node(917), Node(934)},
         {Node(918), Node(928)},  {Node(918), Node(923)},
         {Node(918), Node(935)},  {Node(919), Node(929)},
         {Node(919), Node(924)},  {Node(919), Node(936)},
         {Node(920), Node(930)},  {Node(920), Node(925)},
         {Node(920), Node(937)},  {Node(921), Node(931)},
         {Node(921), Node(926)},  {Node(921), Node(938)},
         {Node(922), Node(932)},  {Node(922), Node(927)},
         {Node(922), Node(939)},  {Node(923), Node(933)},
         {Node(923), Node(928)},  {Node(923), Node(940)},
         {Node(924), Node(934)},  {Node(924), Node(929)},
         {Node(924), Node(941)},  {Node(925), Node(935)},
         {Node(925), Node(930)},  {Node(925), Node(942)},
         {Node(926), Node(936)},  {Node(926), Node(931)},
         {Node(926), Node(943)},  {Node(927), Node(937)},
         {Node(927), Node(932)},  {Node(927), Node(944)},
         {Node(928), Node(938)},  {Node(928), Node(933)},
         {Node(928), Node(945)},  {Node(929), Node(939)},
         {Node(929), Node(934)},  {Node(929), Node(946)},
         {Node(930), Node(940)},  {Node(930), Node(935)},
         {Node(930), Node(947)},  {Node(931), Node(941)},
         {Node(931), Node(936)},  {Node(931), Node(948)},
         {Node(932), Node(942)},  {Node(932), Node(937)},
         {Node(932), Node(949)},  {Node(933), Node(943)},
         {Node(933), Node(938)},  {Node(933), Node(950)},
         {Node(934), Node(944)},  {Node(934), Node(939)},
         {Node(934), Node(951)},  {Node(935), Node(945)},
         {Node(935), Node(940)},  {Node(935), Node(952)},
         {Node(936), Node(946)},  {Node(936), Node(941)},
         {Node(936), Node(953)},  {Node(937), Node(947)},
         {Node(937), Node(942)},  {Node(937), Node(954)},
         {Node(938), Node(948)},  {Node(938), Node(943)},
         {Node(938), Node(955)},  {Node(939), Node(949)},
         {Node(939), Node(944)},  {Node(939), Node(956)},
         {Node(940), Node(950)},  {Node(940), Node(945)},
         {Node(940), Node(957)},  {Node(941), Node(951)},
         {Node(941), Node(946)},  {Node(941), Node(958)},
         {Node(942), Node(952)},  {Node(942), Node(947)},
         {Node(942), Node(959)},  {Node(943), Node(953)},
         {Node(943), Node(948)},  {Node(943), Node(960)},
         {Node(944), Node(954)},  {Node(944), Node(949)},
         {Node(944), Node(961)},  {Node(945), Node(955)},
         {Node(945), Node(950)},  {Node(945), Node(962)},
         {Node(946), Node(956)},  {Node(946), Node(951)},
         {Node(946), Node(963)},  {Node(947), Node(957)},
         {Node(947), Node(952)},  {Node(947), Node(964)},
         {Node(948), Node(958)},  {Node(948), Node(953)},
         {Node(948), Node(965)},  {Node(949), Node(959)},
         {Node(949), Node(954)},  {Node(949), Node(966)},
         {Node(950), Node(960)},  {Node(950), Node(955)},
         {Node(950), Node(967)},  {Node(951), Node(961)},
         {Node(951), Node(956)},  {Node(951), Node(968)},
         {Node(952), Node(962)},  {Node(952), Node(957)},
         {Node(952), Node(969)},  {Node(953), Node(963)},
         {Node(953), Node(958)},  {Node(953), Node(970)},
         {Node(954), Node(964)},  {Node(954), Node(959)},
         {Node(954), Node(971)},  {Node(955), Node(965)},
         {Node(955), Node(960)},  {Node(955), Node(972)},
         {Node(956), Node(966)},  {Node(956), Node(961)},
         {Node(956), Node(973)},  {Node(957), Node(967)},
         {Node(957), Node(962)},  {Node(957), Node(974)},
         {Node(958), Node(968)},  {Node(958), Node(963)},
         {Node(958), Node(975)},  {Node(959), Node(969)},
         {Node(959), Node(964)},  {Node(959), Node(976)},
         {Node(960), Node(970)},  {Node(960), Node(965)},
         {Node(960), Node(977)},  {Node(961), Node(971)},
         {Node(961), Node(966)},  {Node(961), Node(978)},
         {Node(962), Node(972)},  {Node(962), Node(967)},
         {Node(962), Node(979)},  {Node(963), Node(973)},
         {Node(963), Node(968)},  {Node(963), Node(980)},
         {Node(964), Node(974)},  {Node(964), Node(969)},
         {Node(964), Node(981)},  {Node(965), Node(975)},
         {Node(965), Node(970)},  {Node(965), Node(982)},
         {Node(966), Node(976)},  {Node(966), Node(971)},
         {Node(966), Node(983)},  {Node(967), Node(977)},
         {Node(967), Node(972)},  {Node(967), Node(984)},
         {Node(968), Node(978)},  {Node(968), Node(973)},
         {Node(968), Node(985)},  {Node(969), Node(979)},
         {Node(969), Node(974)},  {Node(969), Node(986)},
         {Node(970), Node(980)},  {Node(970), Node(975)},
         {Node(970), Node(987)},  {Node(971), Node(981)},
         {Node(971), Node(976)},  {Node(971), Node(988)},
         {Node(972), Node(982)},  {Node(972), Node(977)},
         {Node(972), Node(989)},  {Node(973), Node(983)},
         {Node(973), Node(978)},  {Node(973), Node(990)},
         {Node(974), Node(984)},  {Node(974), Node(979)},
         {Node(974), Node(991)},  {Node(975), Node(985)},
         {Node(975), Node(980)},  {Node(975), Node(992)},
         {Node(976), Node(986)},  {Node(976), Node(981)},
         {Node(976), Node(993)},  {Node(977), Node(987)},
         {Node(977), Node(982)},  {Node(977), Node(994)},
         {Node(978), Node(988)},  {Node(978), Node(983)},
         {Node(978), Node(995)},  {Node(979), Node(989)},
         {Node(979), Node(984)},  {Node(979), Node(996)},
         {Node(980), Node(990)},  {Node(980), Node(985)},
         {Node(980), Node(997)},  {Node(981), Node(991)},
         {Node(981), Node(986)},  {Node(981), Node(998)},
         {Node(982), Node(992)},  {Node(982), Node(987)},
         {Node(982), Node(999)},  {Node(983), Node(993)},
         {Node(983), Node(988)},  {Node(983), Node(1000)},
         {Node(984), Node(994)},  {Node(984), Node(989)},
         {Node(984), Node(1001)}, {Node(985), Node(995)},
         {Node(985), Node(990)},  {Node(985), Node(1002)},
         {Node(986), Node(996)},  {Node(986), Node(991)},
         {Node(986), Node(1003)}, {Node(987), Node(997)},
         {Node(987), Node(992)},  {Node(987), Node(1004)},
         {Node(988), Node(998)},  {Node(988), Node(993)},
         {Node(988), Node(1005)}, {Node(989), Node(999)},
         {Node(989), Node(994)},  {Node(989), Node(1006)},
         {Node(990), Node(1000)}, {Node(990), Node(995)},
         {Node(990), Node(1007)}, {Node(991), Node(1001)},
         {Node(991), Node(996)},  {Node(991), Node(1008)},
         {Node(992), Node(1002)}, {Node(992), Node(997)},
         {Node(992), Node(1009)}, {Node(993), Node(1003)},
         {Node(993), Node(998)},  {Node(993), Node(1010)},
         {Node(994), Node(1004)}, {Node(994), Node(999)},
         {Node(994), Node(1011)}, {Node(995), Node(1005)},
         {Node(995), Node(1000)}, {Node(995), Node(1012)},
         {Node(996), Node(1006)}, {Node(996), Node(1001)},
         {Node(996), Node(1013)}, {Node(997), Node(1007)},
         {Node(997), Node(1002)}, {Node(997), Node(1014)},
         {Node(998), Node(1008)}, {Node(998), Node(1003)},
         {Node(998), Node(1015)}, {Node(999), Node(1009)},
         {Node(999), Node(1004)}, {Node(999), Node(1016)}});

    aas::PathHandler handler(archi);

    aas::PathHandler handi = handler.construct_acyclic_handler();

    handi.get_distance_matrix();

    handi.get_path_matrix();

    std::list<unsigned> nodes_to_add{
        0,   11,  17,  22,  34,  33,  51,  44,  68,  55,  85,  66,  102,
        77,  119, 88,  136, 99,  153, 110, 170, 121, 187, 132, 204, 143,
        221, 154, 238, 165, 255, 176, 272, 187, 289, 198, 306, 209, 323,
        220, 340, 231, 357, 242, 374, 253, 391, 264, 408, 275, 425, 286,
        442, 297, 459, 308, 476, 319, 493, 330, 510, 341, 527, 352, 544,
        363, 561, 374, 578, 385, 595, 396, 612, 407, 629, 418, 646, 429,
        663, 440, 680, 451, 697, 462, 714, 473, 731, 484, 748};

    aas::SteinerTree st(handi, nodes_to_add, 0);

    st.operations_available(handi);
    st.operations_available(handi);
    st.operations_available(handi);
  }
}
SCENARIO("Check Hamiltonian path construction is correct") {
  GIVEN("1 edge Architecture") {
    std::vector<std::pair<Node, Node>> arch_constructor{
        {Node(0), Node(1)}};  // verbose to avoid ambiguous type
    Architecture arch(arch_constructor);
    std::vector<Node> ham = aas::find_hampath(arch);
    std::vector<Node> correct_ham1{Node(0), Node(1)};
    std::vector<Node> correct_ham2{Node(1), Node(0)};
    REQUIRE((ham == correct_ham1) | (ham == correct_ham2));
  }
  GIVEN("2 edge Architecture") {
    Architecture arch({{Node(0), Node(1)}, {Node(1), Node(2)}});
    std::vector<Node> ham = aas::find_hampath(arch);
    REQUIRE(ham.size() == 3);
    std::vector<Node> correct_ham1{Node(0), Node(1), Node(2)};
    std::vector<Node> correct_ham2{Node(2), Node(1), Node(0)};
    REQUIRE((ham == correct_ham1) | (ham == correct_ham2));
  }
  GIVEN("3 edge line Architecture") {
    Architecture arch(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(2), Node(3)}});
    std::vector<Node> ham = aas::find_hampath(arch);
    REQUIRE(ham.size() == 4);
    std::vector<Node> correct_ham1{Node(0), Node(1), Node(2), Node(3)};
    std::vector<Node> correct_ham2{Node(3), Node(2), Node(1), Node(0)};
    REQUIRE((ham == correct_ham1) | (ham == correct_ham2));
  }
  GIVEN("3 edge star Architecture") {
    Architecture arch(
        {{Node(0), Node(1)}, {Node(1), Node(2)}, {Node(1), Node(3)}});
    REQUIRE_THROWS_AS(aas::find_hampath(arch), aas::NoHamiltonPath);
  }
  GIVEN("5 edge cycle Architecture") {
    Architecture arch(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)}});
    std::vector<Node> ham = aas::find_hampath(arch);
    REQUIRE(ham.size() == 5);
  }
  GIVEN("6 edge star Architecture") {
    Architecture arch(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(0), Node(4)},
         {Node(0), Node(5)}});
    REQUIRE_THROWS_AS(aas::find_hampath(arch), aas::NoHamiltonPath);
  }
  GIVEN("8 edge line Architecture") {
    Architecture arch(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(5), Node(6)},
         {Node(6), Node(7)}});
    std::vector<Node> ham = aas::find_hampath(arch);
    REQUIRE(ham.size() == 8);
  }
  GIVEN("8 edge line shuffled Architecture") {
    Architecture arch(
        {{Node(6), Node(4)},
         {Node(4), Node(2)},
         {Node(2), Node(5)},
         {Node(5), Node(3)},
         {Node(3), Node(1)},
         {Node(1), Node(7)},
         {Node(7), Node(0)}});
    std::vector<Node> ham = aas::find_hampath(arch);
    REQUIRE(ham.size() == 8);
    std::pair<Node, Node> edge;
    for (unsigned i = 1; i < 8; ++i) {
      edge = std::pair<Node, Node>(ham[i - 1], ham[i]);
      bool check = false;
      for (auto e : arch.get_all_edges_vec()) {
        if (edge.first == e.first && edge.second == e.second) {
          check = true;
        } else if (edge.first == e.second && edge.second == e.first) {
          check = true;
        }
      }
      REQUIRE(check);
    }
  }
  GIVEN("20 edge line shuffled Architecture") {
    Architecture arch(
        {{Node(18), Node(0)},  {Node(0), Node(12)},  {Node(12), Node(16)},
         {Node(16), Node(13)}, {Node(13), Node(19)}, {Node(19), Node(4)},
         {Node(4), Node(11)},  {Node(11), Node(7)},  {Node(7), Node(15)},
         {Node(15), Node(10)}, {Node(10), Node(5)},  {Node(5), Node(1)},
         {Node(1), Node(17)},  {Node(17), Node(6)},  {Node(6), Node(8)},
         {Node(8), Node(3)},   {Node(3), Node(9)},   {Node(9), Node(14)},
         {Node(14), Node(2)},  {Node(10), Node(9)},  {Node(7), Node(18)},
         {Node(13), Node(14)}, {Node(0), Node(19)},  {Node(11), Node(16)},
         {Node(3), Node(17)},  {Node(12), Node(6)},  {Node(1), Node(2)},
         {Node(5), Node(4)},   {Node(8), Node(15)},  {Node(8), Node(15)}});
    std::vector<Node> ham = aas::find_hampath(arch);
    REQUIRE(ham.size() == 20);
    std::pair<Node, Node> edge;
    for (unsigned i = 1; i < 20; ++i) {
      edge = std::pair<Node, Node>(ham[i - 1], ham[i]);
      bool check = false;
      for (auto e : arch.get_all_edges_vec()) {
        if (edge.first == e.first && edge.second == e.second) {
          check = true;
        } else if (edge.first == e.second && edge.second == e.first) {
          check = true;
        }
      }
      REQUIRE(check);
    }
  }
  GIVEN("20 edge line shuffled Architecture") {
    const Architecture arch(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(6), Node(7)},
         {Node(7), Node(8)},
         {Node(0), Node(3)},
         {Node(3), Node(6)},
         {Node(1), Node(4)},
         {Node(4), Node(7)},
         {Node(2), Node(5)},
         {Node(5), Node(8)}});
    std::vector<Node> ham = aas::find_hampath(arch);
    REQUIRE(ham.size() == 9);
    std::pair<Node, Node> edge;
    for (unsigned i = 1; i < 9; ++i) {
      edge = std::pair<Node, Node>(ham[i - 1], ham[i]);
      bool check = false;
      for (auto e : arch.get_all_edges_vec()) {
        if (edge.first == e.first && edge.second == e.second) {
          check = true;
        } else if (edge.first == e.second && edge.second == e.first) {
          check = true;
        }
      }
      REQUIRE(check);
    }
  }
  GIVEN("20 edge line shuffled Architecture") {
    const Architecture arch(
        {{Node(0), Node(1)},
         {Node(1), Node(2)},
         {Node(3), Node(4)},
         {Node(4), Node(5)},
         {Node(6), Node(7)},
         {Node(7), Node(8)},
         {Node(0), Node(3)},
         {Node(3), Node(6)},
         {Node(1), Node(4)},
         {Node(4), Node(7)},
         {Node(2), Node(5)},
         {Node(5), Node(8)}});
    std::vector<Node> ham = aas::find_hampath(arch);
    REQUIRE(ham.size() == 9);
    std::pair<Node, Node> edge;
    for (unsigned i = 1; i < 9; ++i) {
      edge = std::pair<Node, Node>(ham[i - 1], ham[i]);
      bool check = false;
      for (auto e : arch.get_all_edges_vec()) {
        if (edge.first == e.first && edge.second == e.second) {
          check = true;
        } else if (edge.first == e.second && edge.second == e.first) {
          check = true;
        }
      }
      REQUIRE(check);
    }
  }
}
SCENARIO("Check iteration order construciton") {
  GIVEN("iteration order - simple example") {
    Architecture arch(
        {{Node(0), Node(1)}, {Node(0), Node(2)}, {Node(0), Node(3)}});

    aas::IterationOrder iter_order(arch);

    std::vector<Node> node_order = iter_order.get_iterationorder();
    std::vector<std::pair<Node, Node>> edgelist = iter_order.get_edgelist();
    REQUIRE(node_order.size() == 4);
    REQUIRE(edgelist.size() == 3);
  }
  GIVEN("iteration order - complex example") {
    Architecture arch(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(1)}});

    aas::IterationOrder iter_order(arch);

    std::vector<Node> node_order = iter_order.get_iterationorder();
    std::vector<std::pair<Node, Node>> edgelist = iter_order.get_edgelist();
    REQUIRE(node_order.size() == 4);
    REQUIRE(edgelist.size() == 3);
  }
  GIVEN("iteration order - complex example II ") {
    Architecture arch(
        {{Node(0), Node(1)},
         {Node(0), Node(2)},
         {Node(0), Node(3)},
         {Node(1), Node(2)},
         {Node(2), Node(3)},
         {Node(3), Node(1)},
         {Node(1), Node(4)},
         {Node(2), Node(5)},
         {Node(3), Node(6)}});

    aas::IterationOrder iter_order(arch);

    std::vector<Node> node_order = iter_order.get_iterationorder();
    std::vector<std::pair<Node, Node>> edgelist = iter_order.get_edgelist();
    REQUIRE(node_order.size() == 7);
    REQUIRE(edgelist.size() == 6);
  }
  GIVEN("iteration order - complex example III ") {
    Architecture arch(
        {{Node(100), Node(10)},
         {Node(100), Node(20)},
         {Node(100), Node(30)},
         {Node(10), Node(20)},
         {Node(20), Node(30)},
         {Node(30), Node(10)},
         {Node(10), Node(40)},
         {Node(20), Node(50)},
         {Node(30), Node(60)}});

    aas::IterationOrder iter_order(arch);

    std::vector<Node> node_order = iter_order.get_iterationorder();
    std::vector<std::pair<Node, Node>> edgelist = iter_order.get_edgelist();
    REQUIRE(node_order.size() == 7);
    REQUIRE(edgelist.size() == 6);
  }
}
}  // namespace tket
