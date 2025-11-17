#include <cassert>
#include <iostream>
#include <string>

#include "tket-c-api.h"

int main() {
  std::string circ_json =
      R"({"bits": [], "commands": [{"args": [["q", [0]], ["q", [1]]], "op": {"type": "CX"}}, {"args": [["q", [1]], ["q", [0]]], "op": {"type": "CX"}}], "created_qubits": [], "discarded_qubits": [], "implicit_permutation": [[["q", [0]], ["q", [0]]], [["q", [1]], ["q", [1]]]], "phase": "0.0", "qubits": [["q", [0]], ["q", [1]]]})";
  TketCircuit *circ = tket_circuit_from_json(circ_json.c_str());
  assert(circ != nullptr);

  std::string pass_json =
      R"({"StandardPass": {"allow_swaps": false, "basis_allowed": ["H", "Rz", "CZ"], "name": "AutoRebase"}, "pass_class": "StandardPass"})";
  TketPass *pass = tket_pass_from_json(pass_json.c_str());
  assert(pass != nullptr);

  TketError rv = tket_apply_pass(circ, pass);
  if (rv != TKET_SUCCESS) {
    std::cerr << "Error applying pass: " << rv << std::endl;
    return -1;
  }

  char *circ1_json = nullptr;
  tket_circuit_to_json(circ, &circ1_json);
  assert(strstr(circ1_json, "CZ"));

  tket_free_circuit(circ);
  tket_free_pass(pass);
  free(circ1_json);

  std::cout << "Success" << std::endl;
  return 0;
}
