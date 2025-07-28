extern "C" {
#include "tket-c-api.h"
}

#include <cstring>
#include <system_error>

#include "tket/Circuit/Circuit.hpp"
#include "tket/Transformations/BasicOptimisation.hpp"
#include "tket/Transformations/OptimisationPass.hpp"

using namespace tket;
using json = nlohmann::json;

struct TketCircuit {
  Circuit circuit;
};

OpType convert_target_gate(TketTargetGate target_gate) {
  switch (target_gate) {
    case TKET_TARGET_CX:
      return OpType::CX;
    case TKET_TARGET_TK2:
      return OpType::TK2;
    default:
      std::cerr << "Invalid target gate\n";
      std::exit(EXIT_FAILURE);
  }
}

TketCircuit *tket_circuit_from_json(const char *json_str) {
  if (!json_str) return nullptr;

  TketCircuit *tc = nullptr;

  // Parse JSON and create circuit
  try {
    tc = new TketCircuit;
    tc->circuit = json::parse(json_str);
  } catch (const json::parse_error &e) {
    std::cerr << "Invalid JSON in tket_circuit_from_json: " << e.what() << "\n";
    if (tc) tket_free_circuit(tc);
    tc = nullptr;
  } catch (...) {
    // Clean up memory, print error, and exit
    if (tc) tket_free_circuit(tc);
    std::cerr << "Unknown error in tket_circuit_from_json\n";
    std::exit(EXIT_FAILURE);
  }

  return tc;
}

TketError tket_circuit_to_json(const TketCircuit *tc, char **json_str) {
  if (!tc || !json_str) return TKET_ERROR_NULL_POINTER;

  std::string s;

  // Convert circuit to JSON
  try {
    json j = tc->circuit;
    s = j.dump();
  } catch (const json::exception &e) {
    // Something went wrong with reading the circuit into JSON
    return TKET_ERROR_CIRCUIT_INVALID;
  }

  // Allocate memory and copy JSON to C string
  try {
    *json_str = (char *)malloc(s.size() + 1);
    std::strcpy(*json_str, s.c_str());
  } catch (...) {
    // Clean up memory, print error, and exit
    if (*json_str) free(*json_str);
    *json_str = nullptr;
    std::cerr << "Unknown error in tket_circuit_from_json\n";
    std::exit(EXIT_FAILURE);
  }

  return TKET_SUCCESS;
}

void tket_free_circuit(TketCircuit *tc) { delete tc; }

void tket_free_string(char *str) { free(str); }

TketError tket_two_qubit_squash(
    TketCircuit *tc, TketTargetGate target_gate, double cx_fidelity,
    bool allow_swaps) {
  if (!tc) return TKET_ERROR_NULL_POINTER;

  Transforms::two_qubit_squash(
      convert_target_gate(target_gate), cx_fidelity, allow_swaps)
      .apply(tc->circuit);

  return TKET_SUCCESS;
}

TketError tket_clifford_simp(
    TketCircuit *tc, TketTargetGate target_gate, bool allow_swaps) {
  if (!tc) return TKET_ERROR_NULL_POINTER;

  Transforms::clifford_simp(allow_swaps, convert_target_gate(target_gate))
      .apply(tc->circuit);

  return TKET_SUCCESS;
}

const char *tket_error_string(TketError error) {
  switch (error) {
    case TKET_SUCCESS:
      return "Success";
    case TKET_ERROR_NULL_POINTER:
      return "Invalid NULL pointer in arguments";
    case TKET_ERROR_CIRCUIT_INVALID:
      return "Invalid circuit: could not convert to JSON";
    default:
      return "Unknown error";
  }
}
