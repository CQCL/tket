#include <catch2/catch.hpp>
#include <numeric>
#include <optional>

#include "Circuit/CircPool.hpp"
#include "Circuit/CircUtils.hpp"
#include "CircuitsForTesting.hpp"
#include "Gate/Rotation.hpp"
#include "OpType/OpType.hpp"

namespace tket {

namespace test_SubcircuitAnalysis {
SCENARIO("Find large subcircuits") {
  GIVEN("A CX circuit") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::CX, {0, 1}, "1");
    circ.add_op<unsigned>(OpType::CX, {2, 3}, "2");
    circ.add_op<unsigned>(OpType::CX, {3, 4}, "3");
    circ.add_op<unsigned>(OpType::CX, {1, 2}, "4");
    circ.add_op<unsigned>(OpType::CX, {2, 3}, "5");
    circ.add_op<unsigned>(OpType::CX, {1, 2}, "6");
    circ.add_op<unsigned>(OpType::CX, {2, 4}, "7");
    circ.add_op<unsigned>(OpType::CX, {2, 4}, "8");
    circ.add_op<unsigned>(OpType::CX, {2, 4}, "9");
    std::vector<Subcircuit> subs = circ.get_subcircuits(4, 3);
    REQUIRE(subs.size() == 2);
    REQUIRE(subs[0].verts.size() == 4);
    REQUIRE(subs[1].verts.size() == 3);
  }
  GIVEN("A multi-q circuit") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::CCX, {0, 1, 2}, "1");
    circ.add_op<unsigned>(OpType::CX, {2, 3}, "2");
    circ.add_op<unsigned>(OpType::CX, {2, 1}, "3");
    circ.add_op<unsigned>(OpType::CX, {0, 1}, "4");
    circ.add_op<unsigned>(OpType::CX, {1, 3}, "5");
    circ.add_op<unsigned>(OpType::CX, {3, 4}, "6");
    circ.add_op<unsigned>(OpType::CCX, {1, 3, 4}, "7");
    circ.add_op<unsigned>(OpType::CX, {2, 4}, "8");
    circ.add_op<unsigned>(OpType::CX, {2, 4}, "9");
    std::vector<Subcircuit> subs = circ.get_subcircuits(4, 3);
    // for (Subcircuit & sub : subs) {
    //     std::cout<<"\nCircuit:";
    //     for (const Vertex & v : sub.verts) {
    //         std::cout<<circ.get_opgroup_from_Vertex(v).value();
    //     }

    // }
    REQUIRE(subs.size() == 2);
    REQUIRE(subs[0].verts.size() == 5);
    REQUIRE(subs[1].verts.size() == 4);
  }
}

}  // namespace test_SubcircuitAnalysis
}  // namespace tket