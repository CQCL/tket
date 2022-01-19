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

from pytket import OpType, Qubit, Bit
from pytket.program import Program  # type: ignore

import pytest  # type: ignore


def test_straight_line_gen() -> None:
    p = Program(2, 2)
    p.add_gate(OpType.X, [0])
    p.add_gate(OpType.Y, [1], condition_bits=[0, 1], condition_value=3)
    p.add_gate(OpType.Measure, [0, 0])

    commands = p.get_commands()
    assert len(commands) == 4
    assert str(commands[0]) == "X q[0];"
    assert str(commands[1]) == "IF ([c[0], c[1]] == 3) THEN Y q[1];"
    assert str(commands[2]) == "Measure q[0] --> c[0];"
    assert str(commands[3]) == "Stop;"


def test_append() -> None:
    p = Program()
    qr = p.add_q_register("qr", 3)
    cr = p.add_c_register("cr", 2)
    p.add_gate(OpType.X, [], [qr[0]])
    p.add_gate(OpType.Measure, [], [qr[2], cr[0]])
    assert len(p.qubits) == 3
    assert len(p.bits) == 2

    p2 = Program()
    anc = Qubit("anc")
    err = Bit("err")
    p2.add_qubit(anc)
    assert len(p2.qubits) == 1
    p2.add_bit(err)
    assert len(p2.bits) == 1

    p2.add_gate(OpType.Measure, [], [anc, err])
    p.append(p2)
    commands = p.get_commands()
    assert len(commands) == 4
    assert str(commands[2]) == "Measure anc --> err;"

    bmap = p.bit_readout
    assert len(bmap) == 3
    assert bmap[err] == 2
    qmap = p.qubit_readout
    assert len(qmap) == 1  # Only the anc-->err measure is in the last block
    assert qmap[anc] == 2


def test_append_if() -> None:
    p = Program(2, 2)
    p.add_gate(OpType.X, [0])

    p2 = Program(2)
    p2.add_gate(OpType.CX, [0, 1])

    p.append_if(Bit(0), p2)
    commands = p.get_commands()
    assert len(commands) == 7
    assert str(commands[0]) == "X q[0];"
    assert str(commands[1]) == "Branch lab_0 c[0];"
    assert str(commands[2]) == "Goto lab_1;"
    assert str(commands[3]) == "Label lab_0;"
    assert str(commands[4]) == "CX q[0], q[1];"
    assert str(commands[5]) == "Label lab_1;"
    assert str(commands[6]) == "Stop;"


def test_append_if_else() -> None:
    p = Program(2, 2)
    p.add_gate(OpType.X, [0])

    p2 = Program(2)
    p2.add_gate(OpType.CX, [0, 1])

    p3 = Program(2)
    p3.add_gate(OpType.CZ, [1, 0])

    p.append_if_else(Bit(0), p2, p3)
    commands = p.get_commands()
    assert len(commands) == 8
    assert str(commands[0]) == "X q[0];"
    assert str(commands[1]) == "Branch lab_0 c[0];"
    assert str(commands[2]) == "CZ q[1], q[0];"
    assert str(commands[3]) == "Goto lab_1;"
    assert str(commands[4]) == "Label lab_0;"
    assert str(commands[5]) == "CX q[0], q[1];"
    assert str(commands[6]) == "Label lab_1;"
    assert str(commands[7]) == "Stop;"


def test_append_while() -> None:
    p = Program(2, 2)
    p.add_gate(OpType.X, [0])

    p2 = Program(2)
    p2.add_gate(OpType.CX, [0, 1])

    p.append_while(Bit(0), p2)
    commands = p.get_commands()
    assert len(commands) == 9
    assert str(commands[0]) == "X q[0];"
    assert str(commands[1]) == "Label lab_0;"
    assert str(commands[2]) == "Branch lab_1 c[0];"
    assert str(commands[3]) == "Goto lab_2;"
    assert str(commands[4]) == "Label lab_1;"
    assert str(commands[5]) == "CX q[0], q[1];"
    assert str(commands[6]) == "Goto lab_0;"
    assert str(commands[7]) == "Label lab_2;"
    assert str(commands[8]) == "Stop;"


if __name__ == "__main__":
    test_straight_line_gen()
    test_append()
    test_append_if()
    test_append_if_else()
    test_append_while()
