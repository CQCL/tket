# Copyright 2019-2022 Cambridge Quantum Computing
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from pytket.circuit import Circuit, OpType  # type: ignore
from pytket.architecture import Architecture  # type: ignore
from pytket.passes import AASRouting, CNotSynthType, ComposePhasePolyBoxes  # type: ignore
from pytket.predicates import CompilationUnit  # type: ignore


def test_AAS() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc, lookahead=2)
    assert pass1.apply(circ)


def test_AAS_2() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc)
    assert pass1.apply(circ)


def test_AAS_3() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc, lookahead=2)
    assert pass1.apply(circ)


def test_AAS_4() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc)
    assert pass1.apply(circ)


def test_AAS_5() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc, lookahead=2)
    assert pass1.apply(circ)


def test_AAS_6() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc)
    assert pass1.apply(circ)


def test_AAS_7() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.H(0).H(2)
    circ.CX(0, 1).CX(1, 2).CX(3, 4)
    circ.Rz(0, 1)
    pass1 = AASRouting(arc, lookahead=2)
    assert pass1.apply(circ)


def test_AAS_8() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    circ = Circuit(5)
    circ.CX(0, 1)
    circ.H(0)
    circ.Z(1)
    circ.CX(0, 3)
    circ.Rx(1.5, 3)
    circ.CX(2, 4)
    circ.X(2)
    circ.CX(1, 4)
    circ.CX(0, 4)
    pass1 = AASRouting(arc, lookahead=2)
    assert pass1.apply(circ)


def test_AAS_9() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8]])
    circ = Circuit(9)
    circ.CX(0, 8).CX(8, 1).CX(1, 7).CX(7, 2).CX(2, 6).CX(6, 3).CX(3, 5).CX(5, 4)
    circ.Rz(0.5, 4)
    pass1 = AASRouting(arc, lookahead=2)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() < 56


def test_AAS_10() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])
    circ = Circuit(7)
    circ.CX(0, 6).CX(6, 1).CX(1, 5).CX(5, 2).CX(2, 4).CX(4, 3)
    circ.Rz(0.5, 3)
    pass1 = AASRouting(arc, lookahead=2)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() < 33


def test_AAS_11() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])
    circ = Circuit(7)
    circ.CX(0, 6).CX(6, 1).CX(1, 5).CX(5, 2).CX(2, 4).CX(4, 3)
    circ.Rz(0.5, 3)
    pass1 = AASRouting(arc, lookahead=1, cnotsynthtype=CNotSynthType.SWAP)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() == 119


def test_AAS_12() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])
    circ = Circuit(7)
    circ.CX(0, 6).CX(6, 1).CX(1, 5).CX(5, 2).CX(2, 4).CX(4, 3)
    circ.Rz(0.5, 3)
    pass1 = AASRouting(arc, lookahead=1, cnotsynthtype=CNotSynthType.HamPath)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() == 36


def test_AAS_13() -> None:
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6]])
    circ = Circuit(7)
    circ.CX(0, 6).CX(6, 1).CX(1, 5).CX(5, 2).CX(2, 4).CX(4, 3)
    circ.Rz(0.5, 3)
    pass1 = AASRouting(arc, lookahead=1, cnotsynthtype=CNotSynthType.Rec)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() == 28


def test_AAS_14() -> None:
    arc = Architecture([[0, 1], [1, 0], [1, 2], [2, 1]])
    circ = Circuit(3).CZ(0, 1)
    pass1 = AASRouting(arc, lookahead=1, cnotsynthtype=CNotSynthType.Rec)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() == 3


def test_AAS_15() -> None:
    arc = Architecture([[0, 1], [1, 0], [1, 2], [2, 1]])
    circ = Circuit(2).CZ(0, 1)
    pass1 = AASRouting(arc, lookahead=1, cnotsynthtype=CNotSynthType.Rec)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.valid_connectivity(arc, False, True)
    assert out_circ.depth() == 3


def test_noncontiguous_arc_phase_poly() -> None:
    # testing non-contiguous ascending named nodes
    arc = Architecture([[0, 2]])
    pass1 = AASRouting(arc, lookahead=1)
    c = Circuit(2).H(0).H(1)
    pass1.apply(c)
    assert c.n_gates_of_type(OpType.H) == 2
    assert c.n_gates_of_type(OpType.CX) == 0
    assert c.n_gates_of_type(OpType.CX) == 0


def test_compose_ppb() -> None:
    circ = Circuit(5).CZ(0, 1).CZ(1, 2).CX(2, 3).CX(3, 4)
    pass1 = ComposePhasePolyBoxes(min_size=2)
    cu = CompilationUnit(circ)
    assert pass1.apply(cu)
    out_circ = cu.circuit
    assert out_circ.depth() == 6


if __name__ == "__main__":
    test_AAS()
    test_AAS_2()
    test_AAS_3()
    test_AAS_4()
    test_AAS_5()
    test_AAS_6()
    test_AAS_7()
    test_AAS_8()
    test_AAS_9()
    test_AAS_10()
    test_AAS_11()
    test_AAS_12()
    test_AAS_13()
    test_AAS_14()
    test_AAS_15()
    test_noncontiguous_arc_phase_poly()
    test_compose_ppb()
