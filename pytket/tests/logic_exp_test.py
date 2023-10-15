from pytket._tket.unit_id import BitRegister
from pytket.circuit.logic_exp import (
    RegOr,
    RegAnd,
    BitOr,
    BitAnd,
    BitXor,
    RegXor,
    RegAdd,
    RegSub,
    RegMul,
    RegDiv,
    RegPow,
    RegLsh,
    RegRsh,
)


def test_types() -> None:
    a = BitRegister("a", 5)
    b = BitRegister("b", 5)
    a0 = a[0]
    b0 = b[0]
    assert isinstance(a0 | b0, BitOr)
    assert isinstance(a0 & b0, BitAnd)
    assert isinstance(a0 ^ b0, BitXor)

    assert isinstance(a | b, RegOr)
    assert isinstance(a & b, RegAnd)
    assert isinstance(a ^ b, RegXor)
    assert isinstance(a + b, RegAdd)
    assert isinstance(a - b, RegSub)
    assert isinstance(a * b, RegMul)
    assert isinstance(a // b, RegDiv)
    assert isinstance(a**b, RegPow)
    assert isinstance(a << b, RegLsh)
    assert isinstance(a >> b, RegRsh)
