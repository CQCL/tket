#include <catch2/catch_test_macros.hpp>

#include "Gate/GatePtr.hpp"
#include "OpType/OpTypeFunctions.hpp"
#include "OpType/OpTypeInfo.hpp"

namespace tket {
namespace test_OpTypeFunctions {

static Gate_ptr get_gate_ptr_from_optype(OpType ot) {
  std::vector<Expr> params(optypeinfo().at(ot).n_params(), 0.);
  return std::make_shared<Gate>(ot, params, 1);
}

SCENARIO("is_single_qubit_unitary_type <=> get_tk1_angles implemented") {
  for (OpType ot : all_single_qubit_types()) {
    Gate_ptr gp = get_gate_ptr_from_optype(ot);
    bool is_single_qb_unitary = true;
    try {
      gp->get_tk1_angles();
    } catch (const NotImplemented&) {
      is_single_qb_unitary = false;
    }
    if (is_single_qb_unitary) {
      REQUIRE(is_single_qubit_unitary_type(ot));
    } else {
      REQUIRE(!is_single_qubit_unitary_type(ot));
    }
  }
}

}  // namespace test_OpTypeFunctions
}  // namespace tket
