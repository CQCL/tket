# Copyright 2019-2024 Cambridge Quantum Computing
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

from typing import Optional, Callable

from pytket.circuit import Circuit
from .._tket.passes import BasePass


class PassSelector:
    """
    Collection of pytket compilation passes which are
    all applied to the same circuit. The result of the
    compilation is the best circuit as selected by a given metric.
    """

    def __init__(
        self,
        passlist: list[BasePass],
        score_func: Callable[[Circuit], int],
    ):
        """
        Constructs a PassSelector

        :param passlist: list of pytket compilation passes
        :param score_func: function to score the
          results of the compilation (lower scores are preferred)
        """
        self._passlist = passlist
        self._score_func = score_func
        if len(self._passlist) < 1:
            raise ValueError("passlist needs to contain at least one pass")

    def apply(self, circ: Circuit) -> Circuit:
        """
        Compiles the given circuit with the best of the given passes.

        :param circ: Circuit that should be compiled
        :return: compiled circuit
        """
        circ_list = [circ.copy() for _ in self._passlist]

        self._scores: list[Optional[int]] = []

        for p, c in zip(self._passlist, circ_list):
            try:
                p.apply(c)
                self._scores.append(self._score_func(c))
            except:  # in case of any error the pass should be ignored
                self._scores.append(None)

        try:
            return circ_list[
                self._scores.index(min(x for x in self._scores if x is not None))
            ]
        except ValueError:
            raise RuntimeError("No passes have successfully run on this circuit")

    def get_scores(self) -> list[Optional[int]]:
        """
        :return: scores of the circuit after compiling
          for each of the compilations passes
        """
        return self._scores
