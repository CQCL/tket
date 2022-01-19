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

from lark import Lark
import pytest
from pytket.circuit import OpType  # type: ignore
from pytket.passes import compilation_pass_from_script, compilation_pass_grammar  # type: ignore
from pytket.passes import (  # type: ignore
    CliffordSimp,
    ContextSimp,
    DecomposeBoxes,
    EulerAngleReduction,
    FullPeepholeOptimise,
    OptimisePhaseGadgets,
    PauliSimp,
    PauliSquash,
    SimplifyInitial,
)
from pytket.passes import BasePass, RepeatPass, SequencePass  # type: ignore
from pytket.transform import CXConfigType, PauliSynthStrat  # type: ignore


def test_grammar() -> None:
    g = compilation_pass_grammar()
    Lark(g)  # will raise exception if grammar invalid


@pytest.mark.parametrize(
    "script,comp_pass",
    [
        (
            "DecomposeBoxes",
            DecomposeBoxes(),
        ),
        (
            "[DecomposeBoxes, FullPeepholeOptimise]",
            SequencePass([DecomposeBoxes(), FullPeepholeOptimise()]),
        ),
        (
            "[DecomposeBoxes, repeat([DecomposeBoxes, DecomposeBoxes]), [DecomposeBoxes]]",
            SequencePass(
                [
                    DecomposeBoxes(),
                    RepeatPass(SequencePass([DecomposeBoxes(), DecomposeBoxes()])),
                    SequencePass([DecomposeBoxes()]),
                ]
            ),
        ),
        (
            "\n \t[ EulerAngleReduction(\nRy , Rz),\r\nFullPeepholeOptimise ]\t",
            SequencePass(
                [EulerAngleReduction(OpType.Ry, OpType.Rz), FullPeepholeOptimise()]
            ),
        ),
        (
            "CliffordSimp",
            CliffordSimp(allow_swaps=True),
        ),
        (
            "CliffordSimpNoSwaps",
            CliffordSimp(allow_swaps=False),
        ),
        (
            "OptimisePhaseGadgets(Tree)",
            OptimisePhaseGadgets(cx_config=CXConfigType.Tree),
        ),
        (
            "OptimisePhaseGadgets",
            OptimisePhaseGadgets(),
        ),
        (
            "PauliSimp(Sets, MultiQGate)",
            PauliSimp(strat=PauliSynthStrat.Sets, cx_config=CXConfigType.MultiQGate),
        ),
        (
            "PauliSquash",
            PauliSquash(),
        ),
        (
            "SimplifyInitial",
            SimplifyInitial(allow_classical=True),
        ),
        (
            "SimplifyInitialNoClassical",
            SimplifyInitial(allow_classical=False),
        ),
        (
            "ContextSimp",
            ContextSimp(),
        ),
        (
            "ContextSimpNoClassical",
            ContextSimp(allow_classical=False),
        ),
        (
            "FullPeepholeOptimise",
            FullPeepholeOptimise(allow_swaps=True),
        ),
        (
            "FullPeepholeOptimiseNoSwaps",
            FullPeepholeOptimise(allow_swaps=False),
        ),
    ],
)
def test_script(script: str, comp_pass: BasePass) -> None:
    p = compilation_pass_from_script(script)
    assert isinstance(p, BasePass)
    assert p.to_dict() == comp_pass.to_dict()


def test_invalid_script() -> None:
    invalid_scripts = [
        "Farewell, cruel world!",
        "EulerAngleReduction(Rz, Rw)",
        "[]",
    ]
    for script in invalid_scripts:
        with pytest.raises(Exception) as errorinfo:
            compilation_pass_from_script(script)
            assert "UnexpectedCharacters" in str(errorinfo.value)
