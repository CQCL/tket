from typing import Callable, Dict, List, Tuple

import pytket._tket.architecture
import pytket._tket.circuit
import pytket._tket.unit_id

class AASLabellingMethod(RoutingMethod):
    def __init__(self) -> None: ...

class AASRouteRoutingMethod(RoutingMethod):
    def __init__(self, aaslookahead: int) -> None: ...

class BoxDecompositionRoutingMethod(RoutingMethod):
    def __init__(self) -> None: ...

class LexiLabellingMethod(RoutingMethod):
    def __init__(self) -> None: ...

class LexiRouteRoutingMethod(RoutingMethod):
    def __init__(self, lookahead: int = ...) -> None: ...

class MappingManager:
    def __init__(
        self, architecture: pytket._tket.architecture.Architecture
    ) -> None: ...
    def route_circuit(
        self,
        circuit: pytket._tket.circuit.Circuit,
        routing_methods: List[RoutingMethod],
    ) -> bool: ...

class MultiGateReorderRoutingMethod(RoutingMethod):
    def __init__(self, max_depth: int = ..., max_size: int = ...) -> None: ...

class RoutingMethod:
    def __init__(self) -> None: ...

class RoutingMethodCircuit(RoutingMethod):
    def __init__(
        self,
        route_subcircuit: Callable[
            [pytket._tket.circuit.Circuit, pytket._tket.architecture.Architecture],
            Tuple[
                bool,
                pytket._tket.circuit.Circuit,
                Dict[pytket._tket.unit_id.UnitID, pytket._tket.unit_id.UnitID],
                Dict[pytket._tket.unit_id.UnitID, pytket._tket.unit_id.UnitID],
            ],
        ],
        max_size: int,
        max_depth: int,
    ) -> None: ...
