from pytket.circuit import Circuit
from pytket.circuit.logic_exp import RegOr, RegWiseOp


def test_som() -> None:
    circ = Circuit(3)
    a = circ.add_c_register("a", 5)
    b = circ.add_c_register("b", 5)
    c = circ.add_c_register("c", 5)
    m = a | b
    assert isinstance(m, RegOr)
    circ.add_classicalexpbox_register(a | b, c.to_list())
    for g in circ:
        assert isinstance(g.op.get_exp(), RegOr)