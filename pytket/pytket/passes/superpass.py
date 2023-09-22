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

import copy
from typing import Optional, Callable

from pytket.circuit import Circuit
from .._tket.passes import BasePass


class superpass:
    """
    collection of different pytket compilation passes which are
    all applied on the same circuit. The result of the
    compilation is best circuit selected by a given metric.
    """

    def __init__(
        self,
        passlist: list[BasePass],
        score_func: Optional[Callable[[Circuit], int]] = None,
    ):
        """
        Constructs a superpass

        :param passlist: list of pytket compilation passes
        :param score_func: function to score the
          results of the compilation, this is optional.
          The default function is the circuit depth.
        """
        self._passlist = passlist
        self._score_func = score_func

    def apply(self, circ: Circuit) -> Circuit:
        """
        Compiles the given circuit with the best of the given passes.

        :param circ: Circuit that should be compiled
        :return: compiled circuit
        """
        circ_list = [copy.copy(circ) for _ in self._passlist]
        for p, c in zip(self._passlist, circ_list):
            p.apply(c)

        if self._score_func is None:
            self._result_size = [c.depth() for c in circ_list]
        else:
            self._result_size = [self._score_func(c) for c in circ_list]

        return circ_list[self._result_size.index(min(self._result_size))]

    def get_result_size(self) -> list[int]:
        """
        :return: scores of the circuit after compiling
          for each of the compilations passes
        """
        return self._result_size
