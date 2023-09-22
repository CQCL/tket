# Copyright 2019-2023 Cambridge Quantum Computing
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

from pytket.circuit import Circuit, OpType
from pytket.passes import FullPeepholeOptimise, superpass


def test_compilation() -> None:
    # TKET-2335
    c = Circuit.from_dict(
        {
            "bits": [["c", [0]], ["c", [1]], ["c", [2]]],
            "commands": [
                {"args": [["q", [2]]], "op": {"type": "X"}},
                {"args": [["q", [3]]], "op": {"type": "SX"}},
                {"args": [["q", [4]]], "op": {"type": "X"}},
                {"args": [["q", [4]], ["q", [0]]], "op": {"type": "CX"}},
                {
                    "args": [["q", [2]]],
                    "op": {"params": ["1.92159003147011"], "type": "Rz"},
                },
                {
                    "args": [["q", [3]]],
                    "op": {"params": ["2.27434149111356"], "type": "Rz"},
                },
                {
                    "args": [["q", [0]]],
                    "op": {"params": ["3.6441046525849"], "type": "Rz"},
                },
                {"args": [["q", [2]]], "op": {"type": "X"}},
                {
                    "args": [["q", [3]]],
                    "op": {"params": ["3.89662172760735"], "type": "Rz"},
                },
                {"args": [["q", [4]]], "op": {"type": "X"}},
                {
                    "args": [["q", [0]]],
                    "op": {"params": ["1.16126164227384"], "type": "Rz"},
                },
                {"args": [["q", [2]]], "op": {"type": "SX"}},
                {"args": [["q", [3]]], "op": {"type": "X"}},
                {"args": [["q", [4]]], "op": {"type": "SX"}},
                {"args": [["q", [2]], ["q", [0]]], "op": {"type": "CX"}},
                {"args": [["q", [3]]], "op": {"type": "SX"}},
                {"args": [["q", [4]]], "op": {"type": "X"}},
                {"args": [["q", [0]]], "op": {"type": "X"}},
                {"args": [["q", [1]], ["q", [2]]], "op": {"type": "CX"}},
                {"args": [["q", [4]]], "op": {"type": "X"}},
                {"args": [["q", [0]]], "op": {"type": "X"}},
                {
                    "args": [["q", [1]]],
                    "op": {"params": ["0.803110840391647"], "type": "Rz"},
                },
                {"args": [["q", [2]]], "op": {"type": "SX"}},
                {
                    "args": [["q", [4]]],
                    "op": {"params": ["2.1139514637727"], "type": "Rz"},
                },
                {
                    "args": [["q", [0]]],
                    "op": {"params": ["2.91558747399915"], "type": "Rz"},
                },
                {
                    "args": [["q", [1]]],
                    "op": {"params": ["3.43170320594144"], "type": "Rz"},
                },
                {
                    "args": [["q", [2]]],
                    "op": {"params": ["0.00131066898193488"], "type": "Rz"},
                },
                {"args": [["q", [4]]], "op": {"type": "X"}},
                {"args": [["q", [2]], ["q", [0]]], "op": {"type": "CX"}},
                {
                    "args": [["q", [1]]],
                    "op": {"params": ["0.792771074284277"], "type": "Rz"},
                },
                {"args": [["q", [0]]], "op": {"type": "SX"}},
                {
                    "args": [["q", [1]]],
                    "op": {"params": ["3.5432014179236"], "type": "Rz"},
                },
                {"args": [["q", [2]]], "op": {"type": "SX"}},
                {
                    "args": [["q", [0]]],
                    "op": {"params": ["3.15399322190124"], "type": "Rz"},
                },
                {"args": [["q", [1]]], "op": {"type": "X"}},
                {"args": [["q", [0]], ["q", [4]]], "op": {"type": "CX"}},
                {"args": [["q", [3]], ["q", [1]]], "op": {"type": "CX"}},
                {"args": [["q", [0]]], "op": {"type": "SX"}},
                {
                    "args": [["q", [1]]],
                    "op": {"params": ["0.993390217293483"], "type": "Rz"},
                },
                {"args": [["q", [4]], ["q", [2]]], "op": {"type": "CX"}},
                {"args": [["q", [3]]], "op": {"type": "SX"}},
                {"args": [["q", [0]]], "op": {"type": "X"}},
                {"args": [["q", [2]]], "op": {"type": "X"}},
                {"args": [["q", [4]]], "op": {"type": "SX"}},
                {"args": [["q", [2]], ["q", [1]]], "op": {"type": "CX"}},
                {"args": [["q", [1]]], "op": {"type": "SX"}},
                {"args": [["q", [4]], ["q", [2]]], "op": {"type": "CX"}},
                {"args": [["q", [1]]], "op": {"type": "SX"}},
                {"args": [["q", [2]]], "op": {"type": "X"}},
                {"args": [["q", [4]], ["q", [3]]], "op": {"type": "CX"}},
                {"args": [["q", [2]], ["q", [3]]], "op": {"type": "CX"}},
                {"args": [["q", [4]]], "op": {"type": "X"}},
                {"args": [["q", [3]], ["c", [1]]], "op": {"type": "Measure"}},
                {
                    "args": [["q", [4]]],
                    "op": {"params": ["2.12070321579766"], "type": "Rz"},
                },
                {"args": [["q", [0]], ["q", [4]]], "op": {"type": "CX"}},
                {"args": [["q", [0]], ["c", [2]]], "op": {"type": "Measure"}},
                {
                    "args": [["q", [4]]],
                    "op": {"params": ["2.08463952265192"], "type": "Rz"},
                },
                {"args": [["q", [4]], ["q", [2]]], "op": {"type": "CX"}},
                {
                    "args": [["q", [2]]],
                    "op": {"params": ["3.23214826870156"], "type": "Rz"},
                },
                {"args": [["q", [4]], ["q", [2]]], "op": {"type": "CX"}},
                {"args": [["q", [2]]], "op": {"type": "X"}},
                {
                    "args": [["q", [4]]],
                    "op": {"params": ["2.02237956848157"], "type": "Rz"},
                },
                {"args": [["q", [4]], ["c", [0]]], "op": {"type": "Measure"}},
            ],
            "implicit_permutation": [
                [["q", [0]], ["q", [0]]],
                [["q", [1]], ["q", [1]]],
                [["q", [2]], ["q", [2]]],
                [["q", [3]], ["q", [3]]],
                [["q", [4]], ["q", [4]]],
            ],
            "phase": "0.0",
            "qubits": [["q", [0]], ["q", [1]], ["q", [2]], ["q", [3]], ["q", [4]]],
        }
    )
    assert FullPeepholeOptimise().apply(c)


def test_superpass() -> None:
    fp = FullPeepholeOptimise()
    fp2 = FullPeepholeOptimise(allow_swaps=False)

    sp = superpass([fp, fp2])

    circ = Circuit(2).H(1).H(0).H(1).H(0).X(1).CX(1, 0).CX(0, 1).CX(1, 0)

    result = sp.apply(circ)

    assert sp.get_result_size() == [1, 4]

    assert result.depth() == min(sp.get_result_size())


def test_superpass_ii() -> None:
    fp = FullPeepholeOptimise()
    fp2 = FullPeepholeOptimise(allow_swaps=False)

    sp = superpass([fp, fp2])

    circ = Circuit(2).H(1).H(0).H(1).H(0).X(1).CX(1, 0).CX(0, 1).CX(1, 0)

    result = sp.apply(circ)

    assert sp.get_result_size() == [1, 4]

    assert result.depth() == min(sp.get_result_size())


def test_superpass_iii() -> None:
    fp = FullPeepholeOptimise()
    fp2 = FullPeepholeOptimise(allow_swaps=False)

    def count_gates(circ: Circuit) -> int:
        return circ.n_gates_of_type(OpType.CX)

    sp = superpass([fp, fp2], count_gates)

    circ = Circuit(2).H(1).H(0).H(1).H(0).X(1).CX(1, 0).CX(0, 1).CX(1, 0)

    result = sp.apply(circ)

    assert sp.get_result_size() == [0, 3]

    assert count_gates(result) == min(sp.get_result_size())


def test_superpass_iv() -> None:
    fp = FullPeepholeOptimise()
    fp2 = FullPeepholeOptimise(allow_swaps=False)

    def count_gates(circ: Circuit) -> int:
        return circ.n_gates_of_type(OpType.X)

    sp = superpass([fp, fp2], count_gates)

    circ = Circuit(2).H(1).H(0).H(1).H(0).X(1).CX(1, 0).CX(0, 1).CX(1, 0)

    result = sp.apply(circ)

    assert sp.get_result_size() == [0, 0]

    assert count_gates(result) == min(sp.get_result_size())
