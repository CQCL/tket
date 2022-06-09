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

#include "Circuit/Circuit.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Transformations/Decomposition.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_GlobalisePhasedX {

static bool is_global(Vertex v, const Circuit& c) {
  unsigned in = c.n_in_edges_of_type(v, EdgeType::Quantum);
  unsigned out = c.n_out_edges_of_type(v, EdgeType::Quantum);
  return in == c.n_qubits() && out == c.n_qubits();
}

SCENARIO("globalise_PhasedX transform") {
  GIVEN("A simple circuit") {
    Circuit c1(2);
    c1.add_op<unsigned>(OpType::NPhasedX, {0.4, 0.3}, {0, 1});
    c1.add_op<unsigned>(OpType::NPhasedX, {0.4, 0.3}, {0, 1});

    WHEN("Applying the transform") {
      Circuit c2 = c1;
      REQUIRE(Transforms::globalise_PhasedX().apply(c2));
      THEN("There are two global NPhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 1);
        BGL_FORALL_VERTICES(v, c2.dag, DAG) {
          OpType type = c2.get_OpType_from_Vertex(v);
          if (type == OpType::NPhasedX) {
            REQUIRE(is_global(v, c2));
          } else if (c2.detect_singleq_unitary_op(v)) {
            // single-qb gates are Rz
            REQUIRE(type == OpType::Rz);
          }
        }
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
    WHEN("Applying the transform (without squashing)") {
      Circuit c2 = c1;
      Transforms::globalise_PhasedX(false).apply(c2);
      THEN("There are two global NPhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
        BGL_FORALL_VERTICES(v, c2.dag, DAG) {
          OpType type = c2.get_OpType_from_Vertex(v);
          if (type == OpType::NPhasedX) {
            REQUIRE(is_global(v, c2));
          } else if (c2.detect_singleq_unitary_op(v)) {
            // single-qb gates are Rz
            REQUIRE(type == OpType::Rz);
          }
        }
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
  }
  GIVEN("Another simple circuit") {
    Circuit c1(2);
    c1.add_op<unsigned>(OpType::CZ, {0, 1});
    c1.add_op<unsigned>(OpType::PhasedX, {0.3, 0.5}, {1});
    c1.add_op<unsigned>(OpType::CZ, {0, 1});

    WHEN("Applying the transform") {
      Circuit c2 = c1;
      Transforms::globalise_PhasedX().apply(c2);
      THEN("There are two global NPhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
        BGL_FORALL_VERTICES(v, c2.dag, DAG) {
          OpType type = c2.get_OpType_from_Vertex(v);
          if (type == OpType::NPhasedX) {
            REQUIRE(is_global(v, c2));
          } else if (c2.detect_singleq_unitary_op(v)) {
            // single-qb gates are Rz
            REQUIRE(type == OpType::Rz);
          }
        }
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
    WHEN("Applying the transform (without squashing)") {
      Circuit c2 = c1;
      Transforms::globalise_PhasedX(false).apply(c2);
      THEN("There are two global NPhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
        BGL_FORALL_VERTICES(v, c2.dag, DAG) {
          OpType type = c2.get_OpType_from_Vertex(v);
          if (type == OpType::NPhasedX) {
            REQUIRE(is_global(v, c2));
          } else if (c2.detect_singleq_unitary_op(v)) {
            // single-qb gates are Rz
            REQUIRE(type == OpType::Rz);
          }
        }
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
  }
  GIVEN("Another simple circuit with parallel NPhasedX") {
    Circuit c1(4);
    c1.add_op<unsigned>(OpType::NPhasedX, {0.3, 0.5}, {0, 1});
    c1.add_op<unsigned>(OpType::NPhasedX, {0.7, 0.8}, {2, 3});
    c1.add_op<unsigned>(OpType::CZ, {0, 1});

    WHEN("Applying the transform") {
      Circuit c2 = c1;
      Transforms::globalise_PhasedX().apply(c2);
      THEN("There are two global NPhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
        BGL_FORALL_VERTICES(v, c2.dag, DAG) {
          OpType type = c2.get_OpType_from_Vertex(v);
          if (type == OpType::NPhasedX) {
            REQUIRE(is_global(v, c2));
          } else if (c2.detect_singleq_unitary_op(v)) {
            // single-qb gates are Rz
            REQUIRE(type == OpType::Rz);
          }
        }
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
    WHEN("Applying the transform (without squashing)") {
      Circuit c2 = c1;
      Transforms::globalise_PhasedX(false).apply(c2);
      THEN("There are two global NPhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
        BGL_FORALL_VERTICES(v, c2.dag, DAG) {
          OpType type = c2.get_OpType_from_Vertex(v);
          if (type == OpType::NPhasedX) {
            REQUIRE(is_global(v, c2));
          } else if (c2.detect_singleq_unitary_op(v)) {
            // single-qb gates are Rz
            REQUIRE(type == OpType::Rz);
          }
        }
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
  }
  GIVEN("A circuit with floating point error") {
    Circuit c1(3);
    c1.add_op<unsigned>(OpType::PhasedX, {0.3, 0.0}, {0});
    c1.add_op<unsigned>(OpType::PhasedX, {0.3, 0.0}, {1});
    c1.add_op<unsigned>(OpType::PhasedX, {0.4, 0.0}, {2});
    c1.add_op<unsigned>(OpType::CZ, {0, 1});
    c1.add_op<unsigned>(OpType::PhasedX, {0.1, 0.0}, {0});
    c1.add_op<unsigned>(OpType::PhasedX, {0.1, 0.0}, {1});

    WHEN("Applying the transform") {
      Circuit c2 = c1;
      Transforms::globalise_PhasedX().apply(c2);
      THEN("There are two global NPhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
        BGL_FORALL_VERTICES(v, c2.dag, DAG) {
          OpType type = c2.get_OpType_from_Vertex(v);
          if (type == OpType::NPhasedX) {
            REQUIRE(is_global(v, c2));
          } else if (c2.detect_singleq_unitary_op(v)) {
            // single-qb gates are Rz
            REQUIRE(type == OpType::Rz);
          }
        }
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
  }
  GIVEN("A circuit with just a global gate") {
    Circuit c1(2);
    c1.add_op<unsigned>(OpType::NPhasedX, {1.231, 4.21}, {0, 1});
    THEN("Nothing happens") {
      Circuit c2 = c1;
      REQUIRE(!Transforms::globalise_PhasedX(false).apply(c2));
      REQUIRE(c2 == c1);
    }
  }
  GIVEN("A circuit with lots of gates to just ignore (or squash)") {
    Circuit c1(3);
    c1.add_op<unsigned>(OpType::NPhasedX, {1.231, 4.21}, {0, 1});
    c1.add_op<unsigned>(OpType::H, {0});
    c1.add_op<unsigned>(OpType::Rz, 0.23, {0});
    c1.add_op<unsigned>(OpType::PhasedX, {0.1231, 5.123}, {0});
    WHEN("GlobalisePhasedX(squash=false)") {
      Circuit c2 = c1;
      REQUIRE(Transforms::globalise_PhasedX(false).apply(c2));
      THEN("There are 3x NPhasedX and other gates are preserved") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 4);
        REQUIRE(c2.count_gates(OpType::H) == 1);
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
    WHEN("GlobalisePhasedX(squash=true)") {
      Circuit c2 = c1;
      REQUIRE(Transforms::globalise_PhasedX(true).apply(c2));
      THEN("Everything got squashed to 2x NPhasedX") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 2);
        REQUIRE(c2.count_gates(OpType::H) == 0);
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
  }
  GIVEN("A slightly more elaborate circuit") {
    Circuit c1(4);
    c1.add_op<unsigned>(OpType::NPhasedX, {0.4, 0.3}, {0, 1, 2});
    c1.add_op<unsigned>(OpType::NPhasedX, {0.4, 0.3}, {0, 1, 2});
    c1.add_op<unsigned>(OpType::CX, {0, 1});
    c1.add_op<unsigned>(OpType::NPhasedX, {0.2, 0.5}, {0, 1, 3});

    WHEN("Applying the transform") {
      Circuit c2 = c1;
      Transforms::globalise_PhasedX().apply(c2);
      THEN("There are three global NPhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 3);
        BGL_FORALL_VERTICES(v, c2.dag, DAG) {
          OpType type = c2.get_OpType_from_Vertex(v);
          if (type == OpType::NPhasedX) {
            REQUIRE(is_global(v, c2));
          } else if (c2.detect_singleq_unitary_op(v)) {
            // single-qb gates are Rz
            REQUIRE(type == OpType::Rz);
          }
        }
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
    WHEN("Applying the transform (without squashing)") {
      Circuit c2 = c1;
      Transforms::globalise_PhasedX(false).apply(c2);
      THEN("There are six global NPhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 6);
        BGL_FORALL_VERTICES(v, c2.dag, DAG) {
          OpType type = c2.get_OpType_from_Vertex(v);
          if (type == OpType::NPhasedX) {
            REQUIRE(is_global(v, c2));
          } else if (c2.detect_singleq_unitary_op(v)) {
            // single-qb gates are Rz
            REQUIRE(type == OpType::Rz);
          }
        }
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
  }
  GIVEN("A serious-ish 4-qb circuit") {
    Circuit c1(4);
    c1.add_op<unsigned>(OpType::NPhasedX, {0.4, 0.3}, {0, 1, 2});
    c1.add_op<unsigned>(OpType::NPhasedX, {0.4, 0.3}, {0, 1, 2});
    c1.add_op<unsigned>(OpType::NPhasedX, {1.2, 0.5}, {0, 1, 2, 3});
    c1.add_op<unsigned>(OpType::CX, {0, 1});
    c1.add_op<unsigned>(OpType::NPhasedX, {0.4, 0.3}, {0, 1, 2});
    c1.add_op<unsigned>(OpType::NPhasedX, {1.4, 0.8}, {0, 1});
    c1.add_op<unsigned>(OpType::CX, {1, 2});
    c1.add_op<unsigned>(OpType::CX, {2, 3});
    c1.add_op<unsigned>(OpType::NPhasedX, {0.2, 0.5}, {0, 1});
    Circuit c2 = c1;

    WHEN("Applying the transform") {
      Transforms::globalise_PhasedX().apply(c2);
      THEN("There are five global NPhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 5);
        BGL_FORALL_VERTICES(v, c2.dag, DAG) {
          OpType type = c2.get_OpType_from_Vertex(v);
          if (type == OpType::NPhasedX) {
            REQUIRE(is_global(v, c2));
          } else if (c2.detect_singleq_unitary_op(v)) {
            // single-qb gates are Rz
            REQUIRE(type == OpType::Rz);
          }
        }
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
    WHEN("Applying the transform (without squashing)") {
      Transforms::globalise_PhasedX(false).apply(c2);
      THEN("There are three global NPhasedX gates") {
        REQUIRE(c2.count_gates(OpType::NPhasedX) == 11);
        BGL_FORALL_VERTICES(v, c2.dag, DAG) {
          OpType type = c2.get_OpType_from_Vertex(v);
          if (type == OpType::NPhasedX) {
            REQUIRE(is_global(v, c2));
          } else if (c2.detect_singleq_unitary_op(v)) {
            // single-qb gates are Rz
            REQUIRE(type == OpType::Rz);
          }
        }
      }
      THEN("The resulting unitaries are equal") {
        auto u1 = tket_sim::get_unitary(c1);
        auto u2 = tket_sim::get_unitary(c2);
        REQUIRE(u1.isApprox(u2));
      }
    }
  }
}

}  // namespace test_GlobalisePhasedX
}  // namespace tket
