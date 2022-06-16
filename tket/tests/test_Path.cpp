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

#include "ArchAwareSynth/Path.hpp"

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
    std::vector<Node> ham = aas::find_hampath(arch);
    REQUIRE(ham.empty());
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
    std::vector<Node> ham = aas::find_hampath(arch);
    REQUIRE(ham.empty());
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
