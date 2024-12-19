// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <fstream>

#include "../testutil.hpp"
#include "tket/Circuit/Multiplexor.hpp"
#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/Circuit/Simulation/CircuitSimulator.hpp"
#include "tket/PauliGraphRefactor/Converters.hpp"
#include "tket/PauliGraphRefactor/PauliGraph.hpp"
#include "tket/Transformations/CliffordReductionPass.hpp"

namespace tket {
namespace test_PauliGraph3 {

using namespace pg;

bool comp_seqs(const std::list<PGOp_ptr>& seq1, std::list<PGOp_ptr> seq2) {
  if (seq1.size() != seq2.size()) return false;
  // If operations commute, they may appear in a different order. We can
  // identify the reordering by searching for each element of seq1 in seq2 and
  // checking they can be reordered through any skipped elements to get to the
  // position in seq1. We iterate through seq1 and for each element we search
  // through seq2 until we find it. If it is found, we need it to commute with
  // every preceding unmatched operation in seq2 to bring it to the front. We
  // pop it from seq2 so it only contains the unmatched operations.
  for (const PGOp_ptr& op1 : seq1) {
    auto it2 = seq2.begin();
    bool found = false;
    while (it2 != seq2.end()) {
      if (*op1 == *(*it2)) {
        found = true;
        seq2.erase(it2);
        break;
      } else if (!op1->commutes_with(*(*it2))) {
        return false;
      }
      ++it2;
    }
    if (!found) return false;
  }
  return true;
}

SCENARIO("Correct creation of refactored PauliGraphs") {
  GIVEN("A Clifford circuit") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Vdg, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A 1qb circuit") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.6, {0});
    circ.add_op<unsigned>(OpType::Ry, 1.2, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A 2qb circuit with no interaction") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.6, {0});
    circ.add_op<unsigned>(OpType::Ry, 1.2, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.2, {1});
    circ.add_op<unsigned>(OpType::Rx, 1.6, {1});
    circ.add_op<unsigned>(OpType::Rz, 1.3, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A 2qb circuit with some anti-commuting interaction") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.2, {1});
    circ.add_op<unsigned>(OpType::XXPhase, 1.1, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A 2qb circuit with some commuting interaction") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::ZZPhase, 1.1, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A 2qb circuit a Clifford-angled ZZPhase") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::ZZPhase, 0.5, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A 1qb circuit with stuff to merge") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.3, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.6, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A 2qb circuit with stuff to merge") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::ZZPhase, 1.1, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.8, {0});
    circ.add_op<unsigned>(OpType::ZZPhase, 1.6, {1, 0});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A circuit with Cliffords and non-Cliffords") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.4, {0});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 1.8, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A dense example") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {1});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {2});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {1});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {2});
    circ.add_op<unsigned>(OpType::Ry, 0.3, {3});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A more interesting example (tof_3)") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::T, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::T, {0});
    circ.add_op<unsigned>(OpType::Tdg, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::CX, {4, 3});
    circ.add_op<unsigned>(OpType::Tdg, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::T, {3});
    circ.add_op<unsigned>(OpType::CX, {4, 3});
    circ.add_op<unsigned>(OpType::Tdg, {3});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::T, {3});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::CX, {2, 4});
    circ.add_op<unsigned>(OpType::T, {2});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {2, 4});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::CX, {1, 4});
    circ.add_op<unsigned>(OpType::Tdg, {4});
    circ.add_op<unsigned>(OpType::CX, {0, 4});
    circ.add_op<unsigned>(OpType::T, {4});
    circ.add_op<unsigned>(OpType::T, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::T, {0});
    circ.add_op<unsigned>(OpType::Tdg, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    circ.add_op<unsigned>(OpType::H, {4});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A circuit with a PauliExpBox") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::ZZPhase, 0.2, {0, 1});
    circ.add_op<unsigned>(OpType::Vdg, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    PauliExpBox peb(SymPauliTensor({Pauli::Y, Pauli::X}, 0.333));
    circ.add_box(peb, {0, 1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("Teleportation") {
    Circuit circ(3, 2);
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    ChoiAPState correct_ap_in(3);
    correct_ap_in.ap_.phase = 0.;
    ChoiAPState correct_ap = circuit_to_choi_apstate(circ);
    correct_ap.normal_form();
    circ.add_measure(0, 0);
    circ.add_measure(1, 1);
    circ.add_conditional_gate<unsigned>(OpType::X, {}, uvec{2}, {1}, 1);
    circ.add_conditional_gate<unsigned>(OpType::Z, {}, uvec{2}, {0}, 1);
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    PGOp_ptr out_tab = pg.get_vertex_PGOp_ptr(*pg.get_output_tableau());
    PGOutputTableau& out_tab_op = dynamic_cast<PGOutputTableau&>(*out_tab);
    out_tab_op.normal_form();
    std::list<PGOp_ptr> sequence = pg.pgop_sequence();
    ChoiMixTableau correct_out_tab = ChoiMixTableau({
        {SpPauliStabiliser({Pauli::X, Pauli::Z, Pauli::X}),
         SpPauliStabiliser(Qubit(0), Pauli::Z)},
        {SpPauliStabiliser({Pauli::Z, Pauli::I, Pauli::I}),
         SpPauliStabiliser(Qubit(0), Pauli::X)},
        {SpPauliStabiliser({Pauli::Z, Pauli::X, Pauli::I}),
         SpPauliStabiliser(Qubit(1), Pauli::Z)},
        {SpPauliStabiliser({Pauli::I, Pauli::Z, Pauli::X}),
         SpPauliStabiliser(Qubit(1), Pauli::X)},
        {SpPauliStabiliser({Pauli::I, Pauli::X, Pauli::Z}),
         SpPauliStabiliser(Qubit(2), Pauli::Z)},
        {SpPauliStabiliser({Pauli::I, Pauli::I, Pauli::X}),
         SpPauliStabiliser(Qubit(2), Pauli::X)},
    });
    std::list<PGOp_ptr> correct_sequence{
        std::make_shared<PGInputTableau>(ChoiMixTableau(3), correct_ap_in),
        std::make_shared<PGMeasure>(
            SpPauliStabiliser({Pauli::Z, Pauli::X, Pauli::I}), Bit(1)),
        std::make_shared<PGMeasure>(
            SpPauliStabiliser({Pauli::X, Pauli::Z, Pauli::X}), Bit(0)),
        std::make_shared<PGConditional>(
            std::make_shared<PGCliffordRot>(
                SpPauliStabiliser({Pauli::I, Pauli::I, Pauli::X}), 2),
            bit_vector_t{Bit(1)}, 1),
        std::make_shared<PGConditional>(
            std::make_shared<PGCliffordRot>(
                SpPauliStabiliser({Pauli::I, Pauli::X, Pauli::Z}), 2),
            bit_vector_t{Bit(0)}, 1),
        std::make_shared<PGOutputTableau>(correct_out_tab, correct_ap)};
    REQUIRE_NOTHROW(pg.verify());
    CHECK(comp_seqs(sequence, correct_sequence));
    THEN("Print diagram to file") {
      std::ofstream dot_file("pauligraph.dot");
      pg.to_graphviz(dot_file);
      dot_file.close();
      remove("pauligraph.dot");
    }
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    PauliGraph res_pg = circuit_to_pauli_graph3(res);
    PGOp_ptr res_out_tab =
        res_pg.get_vertex_PGOp_ptr(*res_pg.get_output_tableau());
    PGOutputTableau& res_out_tab_op =
        dynamic_cast<PGOutputTableau&>(*res_out_tab);
    res_out_tab_op.normal_form();
    std::list<PGOp_ptr> res_sequence = res_pg.pgop_sequence();
    CHECK(comp_seqs(res_sequence, correct_sequence));
    Circuit res_sets = pauli_graph3_to_circuit_sets(pg);
    PauliGraph res_sets_pg = circuit_to_pauli_graph3(res_sets);
    PGOp_ptr res_sets_out_tab =
        res_pg.get_vertex_PGOp_ptr(*res_sets_pg.get_output_tableau());
    PGOutputTableau& res_sets_out_tab_op =
        dynamic_cast<PGOutputTableau&>(*res_sets_out_tab);
    res_sets_out_tab_op.normal_form();
    std::list<PGOp_ptr> res_sets_sequence = res_pg.pgop_sequence();
    CHECK(comp_seqs(res_sets_sequence, correct_sequence));
  }
  GIVEN("A conjugated Reset and Collapse") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    circ.add_op<unsigned>(OpType::Reset, {1});
    circ.add_op<unsigned>(OpType::CY, {0, 2});
    circ.add_op<unsigned>(OpType::Collapse, {2});
    circ.add_op<unsigned>(OpType::ZZMax, {1, 2});
    circ.add_op<unsigned>(OpType::V, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    REQUIRE(res.count_gates(OpType::Reset) == 1);
    ChoiMixTableau circ_tab = circuit_to_cm_tableau(circ);
    ChoiMixTableau res_tab = circuit_to_cm_tableau(res);
    REQUIRE(circ_tab == res_tab);
    Circuit res_sets = pauli_graph3_to_circuit_sets(pg);
    res_sets.decompose_boxes_recursively();
    REQUIRE(res_sets.count_gates(OpType::Reset) == 1);
    ChoiMixTableau res_sets_tab = circuit_to_cm_tableau(res_sets);
    REQUIRE(circ_tab == res_sets_tab);
  }
  GIVEN("Conjugated QControlBox and Multiplexors") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    circ.add_op<unsigned>(OpType::CY, {1, 3});
    circ.add_op<unsigned>(
        std::make_shared<QControlBox>(
            get_op_ptr(OpType::ISWAPMax), 2, std::vector<bool>{0, 1}),
        {2, 1, 0, 3});
    ctrl_op_map_t op_map = {
        {{0, 0}, get_op_ptr(OpType::CX)},
        {{0, 1}, get_op_ptr(OpType::Sycamore)},
        {{1, 1}, get_op_ptr(OpType::XXPhase, 0.34)}};
    circ.add_op<unsigned>(
        std::make_shared<MultiplexorBox>(op_map), {0, 3, 1, 2});
    ctrl_op_map_t rot_map = {
        {{0, 0, 0}, get_op_ptr(OpType::Rx, 0.14)},
        {{1, 0, 1}, get_op_ptr(OpType::Rx, -1.45)}};
    circ.add_op<unsigned>(
        std::make_shared<MultiplexedRotationBox>(rot_map), {3, 0, 2, 1});
    ctrl_op_map_t u2_map = {
        {{0}, get_op_ptr(OpType::Rz, -0.87)},
        {{1}, get_op_ptr(OpType::TK1, {0.98, -0.12, 1.2})}};
    circ.add_op<unsigned>(std::make_shared<MultiplexedU2Box>(u2_map), {1, 3});
    ctrl_tensored_op_map_t tensor_map = {
        {{0}, {get_op_ptr(OpType::Ry, 0.98), get_op_ptr(OpType::H)}},
        {{1}, {get_op_ptr(OpType::Rx, -0.87), get_op_ptr(OpType::Vdg)}}};
    circ.add_op<unsigned>(
        std::make_shared<MultiplexedTensoredU2Box>(tensor_map), {3, 2, 1});
    circ.add_op<unsigned>(OpType::CY, {0, 2});
    circ.add_op<unsigned>(OpType::ZZMax, {1, 2});
    circ.add_op<unsigned>(OpType::V, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(res.count_gates(OpType::QControlBox) == 1);
      REQUIRE(res.count_gates(OpType::MultiplexorBox) == 1);
      REQUIRE(res.count_gates(OpType::MultiplexedRotationBox) == 1);
      REQUIRE(res.count_gates(OpType::MultiplexedU2Box) == 1);
      REQUIRE(res.count_gates(OpType::MultiplexedTensoredU2Box) == 1);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      res.decompose_boxes_recursively(
          {OpType::QControlBox, OpType::MultiplexorBox,
           OpType::MultiplexedRotationBox, OpType::MultiplexedU2Box,
           OpType::MultiplexedTensoredU2Box});
      REQUIRE(res.count_gates(OpType::QControlBox) == 1);
      REQUIRE(res.count_gates(OpType::MultiplexorBox) == 1);
      REQUIRE(res.count_gates(OpType::MultiplexedRotationBox) == 1);
      REQUIRE(res.count_gates(OpType::MultiplexedU2Box) == 1);
      REQUIRE(res.count_gates(OpType::MultiplexedTensoredU2Box) == 1);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A conjugated Box") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    circ.add_op<unsigned>(OpType::Sycamore, {1, 2});
    circ.add_op<unsigned>(OpType::CY, {0, 2});
    circ.add_op<unsigned>(OpType::ZZMax, {1, 2});
    circ.add_op<unsigned>(OpType::V, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(res.count_gates(OpType::Sycamore) == 1);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(res.count_gates(OpType::Sycamore) == 1);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("Some end-of-circuit measurements") {
    Circuit circ(3, 2);
    circ.add_op<unsigned>(OpType::U1, 1.35, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_measure(0, 0);
    circ.add_measure(2, 1);
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    PGOp_ptr out_tab = pg.get_vertex_PGOp_ptr(*pg.get_output_tableau());
    PGOutputTableau& out_tab_op = dynamic_cast<PGOutputTableau&>(*out_tab);
    out_tab_op.normal_form();
    REQUIRE_NOTHROW(pg.verify());
    SpPauliStabiliser zzi(DensePauliMap{Pauli::Z, Pauli::Z});
    ChoiMixTableau correct_out_tab = ChoiMixTableau({
        {zzi, SpPauliStabiliser(Qubit(0), Pauli::Z)},
        {SpPauliStabiliser(Qubit(0), Pauli::X),
         SpPauliStabiliser(Qubit(0), Pauli::X)},
        {SpPauliStabiliser(Qubit(1), Pauli::Z),
         SpPauliStabiliser(Qubit(1), Pauli::Z)},
        {SpPauliStabiliser({Pauli::X, Pauli::X, Pauli::I}),
         SpPauliStabiliser(Qubit(1), Pauli::X)},
        {SpPauliStabiliser(Qubit(2), Pauli::Z),
         SpPauliStabiliser(Qubit(2), Pauli::Z)},
        {SpPauliStabiliser(Qubit(2), Pauli::X),
         SpPauliStabiliser(Qubit(2), Pauli::X)},
    });
    Circuit cliff_circ(3);
    cliff_circ.add_op<unsigned>(OpType::CX, {1, 0});
    ChoiAPState correct_ap = circuit_to_choi_apstate(cliff_circ);
    correct_ap.normal_form();
    ChoiAPState correct_ap_in(3);
    correct_ap_in.ap_.phase = 0.;
    std::list<PGOp_ptr> correct_sequence{
        std::make_shared<PGInputTableau>(ChoiMixTableau(3), correct_ap_in),
        std::make_shared<PGMeasure>(
            SpPauliStabiliser(Qubit(2), Pauli::Z), Bit(1)),
        std::make_shared<PGRotation>(
            SpPauliStabiliser(Qubit(0), Pauli::Z), 1.35),
        std::make_shared<PGMeasure>(zzi, Bit(0)),
        std::make_shared<PGOutputTableau>(correct_out_tab, correct_ap)};
    CHECK(comp_seqs(pg.pgop_sequence(), correct_sequence));
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      CHECK(res.count_gates(OpType::Measure) == 2);
      PauliGraph res_pg = circuit_to_pauli_graph3(res);
      PGOp_ptr res_out_tab =
          res_pg.get_vertex_PGOp_ptr(*res_pg.get_output_tableau());
      PGOutputTableau& res_out_tab_op =
          dynamic_cast<PGOutputTableau&>(*res_out_tab);
      res_out_tab_op.normal_form();
      CHECK(comp_seqs(res_pg.pgop_sequence(), correct_sequence));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      CHECK(res.count_gates(OpType::Measure) == 2);
      PauliGraph res_pg = circuit_to_pauli_graph3(res);
      PGOp_ptr res_out_tab =
          res_pg.get_vertex_PGOp_ptr(*res_pg.get_output_tableau());
      PGOutputTableau& res_out_tab_op =
          dynamic_cast<PGOutputTableau&>(*res_out_tab);
      res_out_tab_op.normal_form();
      CHECK(comp_seqs(res_pg.pgop_sequence(), correct_sequence));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      CHECK(res.count_gates(OpType::Measure) == 2);
      PauliGraph res_pg = circuit_to_pauli_graph3(res);
      PGOp_ptr res_out_tab =
          res_pg.get_vertex_PGOp_ptr(*res_pg.get_output_tableau());
      PGOutputTableau& res_out_tab_op =
          dynamic_cast<PGOutputTableau&>(*res_out_tab);
      res_out_tab_op.normal_form();
      CHECK(comp_seqs(res_pg.pgop_sequence(), correct_sequence));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      CHECK(res.count_gates(OpType::Measure) == 2);
      PauliGraph res_pg = circuit_to_pauli_graph3(res);
      PGOp_ptr res_out_tab =
          res_pg.get_vertex_PGOp_ptr(*res_pg.get_output_tableau());
      PGOutputTableau& res_out_tab_op =
          dynamic_cast<PGOutputTableau&>(*res_out_tab);
      res_out_tab_op.normal_form();
      CHECK(comp_seqs(res_pg.pgop_sequence(), correct_sequence));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      res.decompose_boxes_recursively();
      CHECK(res.count_gates(OpType::Measure) == 2);
      PauliGraph res_pg = circuit_to_pauli_graph3(res);
      PGOp_ptr res_out_tab =
          res_pg.get_vertex_PGOp_ptr(*res_pg.get_output_tableau());
      PGOutputTableau& res_out_tab_op =
          dynamic_cast<PGOutputTableau&>(*res_out_tab);
      res_out_tab_op.normal_form();
      CHECK(comp_seqs(res_pg.pgop_sequence(), correct_sequence));
    }
  }
  GIVEN("Some stabiliser assertions") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::Rz, 1.5, {0});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PauliStabiliser pauli1 = {{Pauli::X, Pauli::X}, 0};
    PauliStabiliser pauli2 = {{Pauli::Z, Pauli::Z}, 0};
    PauliStabiliser pauli3 = {{Pauli::Y, Pauli::Y}, 2};
    PauliStabiliserVec stabilisers = {pauli1, pauli2, pauli3};
    StabiliserAssertionBox box(stabilisers);
    circ.add_assertion(box, {Qubit(0), Qubit(2)}, Qubit(1));
    circ.add_assertion(box, {Qubit(0), Qubit(2)}, Qubit(1));
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    PGOp_ptr out_tab = pg.get_vertex_PGOp_ptr(*pg.get_output_tableau());
    PGOutputTableau& out_tab_op = dynamic_cast<PGOutputTableau&>(*out_tab);
    out_tab_op.normal_form();
    REQUIRE_NOTHROW(pg.verify());
    std::list<PGOp_ptr> sequence = pg.pgop_sequence();
    SpPauliStabiliser anc_z(Qubit(1), Pauli::Z);
    SpPauliStabiliser anc_x(DensePauliMap{Pauli::X, Pauli::X});
    ChoiMixTableau correct_out_tab = ChoiMixTableau({
        {SpPauliStabiliser({Pauli::Z, Pauli::Z, Pauli::I}),
         SpPauliStabiliser(Qubit(0), Pauli::Z)},
        {SpPauliStabiliser(Qubit(0), Pauli::X),
         SpPauliStabiliser(Qubit(0), Pauli::X)},
        {SpPauliStabiliser(Qubit(1), Pauli::Z),
         SpPauliStabiliser(Qubit(1), Pauli::Z)},
        {SpPauliStabiliser({Pauli::X, Pauli::X, Pauli::I}),
         SpPauliStabiliser(Qubit(1), Pauli::X)},
        {SpPauliStabiliser(Qubit(2), Pauli::Z),
         SpPauliStabiliser(Qubit(2), Pauli::Z)},
        {SpPauliStabiliser(Qubit(2), Pauli::X),
         SpPauliStabiliser(Qubit(2), Pauli::X)},
    });
    Circuit cliff_circ(3);
    cliff_circ.add_op<unsigned>(OpType::CX, {1, 0});
    ChoiAPState correct_ap = circuit_to_choi_apstate(cliff_circ);
    correct_ap.normal_form();
    // Phase from the Rz gate gets added here
    correct_ap.ap_.phase -= 0.75;
    ChoiAPState correct_ap_in(3);
    std::list<PGOp_ptr> correct_sequence{
        std::make_shared<PGInputTableau>(ChoiMixTableau(3), correct_ap_in),
        std::make_shared<PGRotation>(
            SpPauliStabiliser(Qubit(0), Pauli::Z), 1.5),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::X, Pauli::I, Pauli::X}), anc_z, anc_x,
            Bit(c_debug_zero_prefix() + "_" + c_debug_default_name(), 0)),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::Z, Pauli::Z, Pauli::Z}), anc_z, anc_x,
            Bit(c_debug_zero_prefix() + "_" + c_debug_default_name(), 1)),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::Y, Pauli::Z, Pauli::Y}, 2), anc_z, anc_x,
            Bit(c_debug_one_prefix() + "_" + c_debug_default_name(), 0)),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::X, Pauli::I, Pauli::X}), anc_z, anc_x,
            Bit(c_debug_zero_prefix() + "_" + c_debug_default_name() + "(1)",
                0)),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::Z, Pauli::Z, Pauli::Z}), anc_z, anc_x,
            Bit(c_debug_zero_prefix() + "_" + c_debug_default_name() + "(1)",
                1)),
        std::make_shared<PGStabAssertion>(
            SpPauliStabiliser({Pauli::Y, Pauli::Z, Pauli::Y}, 2), anc_z, anc_x,
            Bit(c_debug_one_prefix() + "_" + c_debug_default_name() + "(1)",
                0)),
        std::make_shared<PGOutputTableau>(correct_out_tab, correct_ap)};
    CHECK(comp_seqs(sequence, correct_sequence));
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(res.count_gates(OpType::StabiliserAssertionBox) == 6);
      PauliGraph res_pg = circuit_to_pauli_graph3(res);
      PGOp_ptr res_out_tab =
          res_pg.get_vertex_PGOp_ptr(*res_pg.get_output_tableau());
      PGOutputTableau& res_out_tab_op =
          dynamic_cast<PGOutputTableau&>(*res_out_tab);
      res_out_tab_op.normal_form();
      std::list<PGOp_ptr> res_sequence = res_pg.pgop_sequence();
      CHECK(comp_seqs(res_sequence, correct_sequence));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      res.decompose_boxes_recursively({OpType::StabiliserAssertionBox});
      REQUIRE(res.count_gates(OpType::StabiliserAssertionBox) == 6);
      PauliGraph res_pg = circuit_to_pauli_graph3(res);
      PGOp_ptr res_out_tab =
          res_pg.get_vertex_PGOp_ptr(*res_pg.get_output_tableau());
      PGOutputTableau& res_out_tab_op =
          dynamic_cast<PGOutputTableau&>(*res_out_tab);
      res_out_tab_op.normal_form();
      std::list<PGOp_ptr> res_sequence = res_pg.pgop_sequence();
      CHECK(comp_seqs(res_sequence, correct_sequence));
    }
  }
  GIVEN("A circuit with initialisations, discards, and implicit permutations") {
    Circuit circ(4);
    circ.qubit_create(Qubit(0));
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {3});
    circ.add_op<unsigned>(OpType::CX, {1, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 1});
    circ.qubit_discard(Qubit(2));
    REQUIRE(Transforms::clifford_reduction(true).apply(circ));
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      res.decompose_boxes_recursively();
      CHECK(res.created_qubits() == qubit_vector_t{Qubit(0)});
      CHECK(res.discarded_qubits() == qubit_vector_t{Qubit(2)});
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      res.decompose_boxes_recursively();
      CHECK(res.created_qubits() == qubit_vector_t{Qubit(0)});
      CHECK(res.discarded_qubits() == qubit_vector_t{Qubit(2)});
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("A symbolic circuit") {
    SymEngine::map_basic_basic sub_map;
    Sym a = SymEngine::symbol("a");
    sub_map[a] = Expr(0.8);
    Sym b = SymEngine::symbol("b");
    sub_map[b] = Expr(1.3);
    Sym c = SymEngine::symbol("c");
    sub_map[c] = Expr(-0.2);
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    circ.add_op<unsigned>(OpType::FSim, {Expr(a), Expr(b)}, {1, 2});
    circ.add_op<unsigned>(OpType::CY, {0, 2});
    circ.add_op<unsigned>(OpType::ZZPhase, {2 * Expr(c)}, {1, 2});
    circ.add_op<unsigned>(OpType::V, {1});
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    REQUIRE_NOTHROW(pg.verify());
    CHECK(pg.free_symbols() == SymSet{a, b, c});
    pg.symbol_substitution(sub_map);
    REQUIRE_NOTHROW(pg.verify());
    CHECK_FALSE(pg.is_symbolic());
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(res.count_gates(OpType::FSim) == 1);
      circ.symbol_substitution(sub_map);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(res.count_gates(OpType::FSim) == 1);
      circ.symbol_substitution(sub_map);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
  GIVEN("Don't collect cliffords") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::Y, {0});
    circ.add_op<unsigned>(OpType::Sdg, {1});
    circ.add_op<unsigned>(OpType::V, {2});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CY, {2, 0});
    circ.add_op<unsigned>(OpType::PhaseGadget, 0.198, {0, 1, 2});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::Vdg, {2});
    circ.add_op<unsigned>(OpType::CZ, {1, 2});
    circ.add_op<unsigned>(OpType::ZZMax, {1, 2});
    circ.add_op<unsigned>(OpType::SWAP, {0, 2});
    circ.add_op<unsigned>(OpType::YYPhase, 1.387, {0, 1});
    circ.add_op<unsigned>(OpType::TK1, {0.98, 0.2, 1.87}, {1});
    circ.add_op<unsigned>(OpType::TK2, {1.34, 0.23, 1.42}, {1, 0});
    PauliGraph pg = circuit_to_pauli_graph3(circ, false);
    REQUIRE_NOTHROW(pg.verify());
    WHEN("Legacy individual synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_individually(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy pairwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("Legacy setwise synthesis") {
      Circuit res = pauli_graph3_to_pauli_exp_box_circuit_sets(pg);
      res.decompose_boxes_recursively();
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General individual synthesis") {
      Circuit res = pauli_graph3_to_circuit_individual(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
    WHEN("General setwise synthesis") {
      Circuit res = pauli_graph3_to_circuit_sets(pg);
      REQUIRE(test_unitary_comparison(circ, res, false));
    }
  }
}

SCENARIO("Check global phase in converters") {
  std::list<std::pair<OpType, std::vector<unsigned>>> test_gates = {
      {OpType::Z, {0}},       {OpType::X, {0}},
      {OpType::Y, {0}},       {OpType::S, {0}},
      {OpType::Sdg, {0}},     {OpType::V, {0}},
      {OpType::Vdg, {0}},     {OpType::SX, {0}},
      {OpType::SXdg, {0}},    {OpType::H, {0}},
      {OpType::CX, {0, 1}},   {OpType::CY, {0, 1}},
      {OpType::CZ, {0, 1}},   {OpType::ZZMax, {0, 1}},
      {OpType::ECR, {0, 1}},  {OpType::ISWAPMax, {0, 1}},
      {OpType::SWAP, {0, 1}}, {OpType::BRIDGE, {0, 1, 2}},
      {OpType::noop, {0}},    {OpType::T, {0}},
      {OpType::Tdg, {0}},
  };
  for (const auto& com : test_gates) {
    Circuit circ(3);
    circ.add_op<unsigned>(com.first, com.second);
    PauliGraph pg = circuit_to_pauli_graph3(circ, false);
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    CHECK(test_unitary_comparison(circ, res, false));
  }
  std::list<std::tuple<OpType, std::vector<Expr>, std::vector<unsigned>>>
      test_param_gates = {
          {OpType::Rz, {0.35}, {0}},
          {OpType::U1, {0.35}, {0}},
          {OpType::Rx, {0.35}, {0}},
          {OpType::Ry, {0.35}, {0}},
          {OpType::TK1, {0.35, 0.98, 1.72}, {0}},
          {OpType::PhaseGadget, {0.35}, {0, 1, 2}},
          {OpType::ZZPhase, {0.35}, {0, 1}},
          {OpType::XXPhase, {0.35}, {0, 1}},
          {OpType::YYPhase, {0.35}, {0, 1}},
          {OpType::TK2, {0.35, 0.98, 1.72}, {0, 1}},
      };
  for (const auto& com : test_param_gates) {
    Circuit circ(3);
    circ.add_op<unsigned>(std::get<0>(com), std::get<1>(com), std::get<2>(com));
    PauliGraph pg = circuit_to_pauli_graph3(circ, false);
    Circuit res = pauli_graph3_to_circuit_individual(pg);
    CHECK(test_unitary_comparison(circ, res, false));
  }
}

}  // namespace test_PauliGraph3
}  // namespace tket
