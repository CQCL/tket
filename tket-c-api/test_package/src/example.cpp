#include <cassert>
#include <iostream>
#include <string>

#include "tket-c-api.h"

int main() {
  std::string json_str =
      R"({"bits": [], "commands": [{"args": [["q", [0]], ["q", [1]]], "op": {"type": "CX"}}, {"args": [["q", [0]], ["q", [1]]], "op": {"type": "CX"}}], "created_qubits": [], "discarded_qubits": [], "implicit_permutation": [[["q", [0]], ["q", [0]]], [["q", [1]], ["q", [1]]]], "phase": "0.0", "qubits": [["q", [0]], ["q", [1]]]})";
  TketCircuit *circ = tket_circuit_from_json(json_str.c_str());
  assert(circ != nullptr);
  tket_free_circuit(circ);

  std::cout << "Success" << std::endl;
}
