from typing import Dict, List, Set

import pytket._tket.circuit
import pytket._tket.pauli

class FrameRandomisation:
    def __init__(self, arg0: Set[pytket._tket.circuit.OpType], arg1: Set[pytket._tket.circuit.OpType], arg2: Dict[pytket._tket.circuit.OpType,Dict[tuple,tuple]]) -> None: ...
    def get_all_circuits(self, circuit: pytket._tket.circuit.Circuit) -> List[pytket._tket.circuit.Circuit]: ...
    def sample_circuits(self, circuit: pytket._tket.circuit.Circuit, samples: int) -> List[pytket._tket.circuit.Circuit]: ...

class PauliFrameRandomisation:
    def __init__(self) -> None: ...
    def get_all_circuits(self, circuit: pytket._tket.circuit.Circuit) -> List[pytket._tket.circuit.Circuit]: ...
    def sample_circuits(self, circuit: pytket._tket.circuit.Circuit, samples: int) -> List[pytket._tket.circuit.Circuit]: ...

class UniversalFrameRandomisation:
    def __init__(self) -> None: ...
    def get_all_circuits(self, circuit: pytket._tket.circuit.Circuit) -> List[pytket._tket.circuit.Circuit]: ...
    def sample_circuits(self, circuit: pytket._tket.circuit.Circuit, samples: int) -> List[pytket._tket.circuit.Circuit]: ...

def apply_clifford_basis_change(pauli: pytket._tket.pauli.QubitPauliString, circuit: pytket._tket.circuit.Circuit) -> pytket._tket.pauli.QubitPauliString: ...
