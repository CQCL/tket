from numpy.typing import NDArray
from typing import Any, List

from typing import overload
import numpy
import pytket._tket.circuit
import pytket._tket.unit_id

class UnitaryTableau:
    @overload
    def __init__(self, nqb: int) -> None: ...
    @overload
    def __init__(
        self,
        xx: NDArray[numpy.bool_],
        xz: NDArray[numpy.bool_],
        xph: NDArray[numpy.bool_],
        zx: NDArray[numpy.bool_],
        zz: NDArray[numpy.bool_],
        zph: NDArray[numpy.bool_],
    ) -> None: ...
    def apply_gate_at_end(
        self, arg0: pytket._tket.circuit.OpType, arg1: List[pytket._tket.unit_id.Qubit]
    ) -> None: ...
    def apply_gate_at_front(
        self, arg0: pytket._tket.circuit.OpType, arg1: List[pytket._tket.unit_id.Qubit]
    ) -> None: ...
    def get_row_product(self, *args: Any, **kwargs: Any) -> Any: ...
    def get_xrow(self, *args: Any, **kwargs: Any) -> Any: ...
    def get_zrow(self, *args: Any, **kwargs: Any) -> Any: ...

class UnitaryTableauBox(pytket._tket.circuit.Op):
    @overload
    def __init__(self, tab: UnitaryTableau) -> None: ...
    @overload
    def __init__(
        self,
        xx: NDArray[numpy.bool_],
        xz: NDArray[numpy.bool_],
        xph: NDArray[numpy.bool_],
        zx: NDArray[numpy.bool_],
        zz: NDArray[numpy.bool_],
        zph: NDArray[numpy.bool_],
    ) -> None: ...
    def get_circuit(self) -> pytket._tket.circuit.Circuit: ...
    def get_tableau(self) -> UnitaryTableau: ...
