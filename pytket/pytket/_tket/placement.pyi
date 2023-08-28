from typing import Any
from typing import Dict, List, Tuple

import pytket._tket.architecture
import pytket._tket.circuit
import pytket._tket.unit_id

class GraphPlacement(Placement):
    def __init__(
        self,
        arc: pytket._tket.architecture.Architecture,
        maximum_matches: int = ...,
        timeout: int = ...,
        maximum_pattern_gates: int = ...,
        maximum_pattern_depth: int = ...,
    ) -> None: ...
    def modify_config(self, **kwargs: Any) -> None: ...

class LinePlacement(Placement):
    def __init__(
        self,
        arc: pytket._tket.architecture.Architecture,
        maximum_line_gates: int = ...,
        maximum_line_depth: int = ...,
    ) -> None: ...

class NoiseAwarePlacement(Placement):
    def __init__(
        self,
        arc: pytket._tket.architecture.Architecture,
        node_errors: Dict[pytket._tket.unit_id.Node, float] = ...,
        link_errors: Dict[
            Tuple[pytket._tket.unit_id.Node, pytket._tket.unit_id.Node], float
        ] = ...,
        readout_errors: Dict[pytket._tket.unit_id.Node, float] = ...,
        maximum_matches: int = ...,
        timeout: int = ...,
        maximum_pattern_gates: int = ...,
        maximum_pattern_depth: int = ...,
    ) -> None: ...
    def modify_config(self, **kwargs: Any) -> None: ...

class Placement:
    def __init__(self, arc: pytket._tket.architecture.Architecture) -> None: ...
    @classmethod
    def from_dict(cls, arg0: dict) -> Placement: ...
    def get_placement_map(
        self, circuit: pytket._tket.circuit.Circuit
    ) -> Dict[pytket._tket.unit_id.Qubit, pytket._tket.unit_id.Node]: ...
    def get_placement_maps(
        self, circuit: pytket._tket.circuit.Circuit, matches: int = ...
    ) -> List[Dict[pytket._tket.unit_id.Qubit, pytket._tket.unit_id.Node]]: ...
    def place(self, circuit: pytket._tket.circuit.Circuit) -> bool: ...
    @classmethod
    def place_with_map(
        cls,
        circuit: pytket._tket.circuit.Circuit,
        qmap: Dict[pytket._tket.unit_id.Qubit, pytket._tket.unit_id.Node],
    ) -> bool: ...
    def to_dict(self) -> object: ...

def place_fully_connected(
    circuit: pytket._tket.circuit.Circuit,
    fully_connected: pytket._tket.architecture.FullyConnected,
) -> None: ...
def place_with_map(
    circuit: pytket._tket.circuit.Circuit,
    qmap: Dict[pytket._tket.unit_id.Qubit, pytket._tket.unit_id.Node],
) -> None: ...
