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

from splitting import Splitter
from knitting import Knitter
from mpattern import MPattern
from pytket.circuit import Circuit, Qubit  # type: ignore
from pytket.architecture import Architecture  # type: ignore
from typing import Dict, Tuple


def format_circ(circ: Circuit) -> Tuple[Circuit, Dict[str, Dict[Qubit, int]]]:
    """
    Takes a pytket circuit object as input and returns a tuple containing the circuit as
    well as a dictionary mapping the original input and output qubits to their indices
    in the qubit register.
    
    :param g:       A zx diagram with some remaining simple edges we want to remove.
    :param type:    GraphS
    
    :param io_map:  A dictionary containing the current i/o mapping.
    :param type:    dict
    
    :returns:        A tuple containing the circuit as well as its corresponding i/o map.
    :rtype:          Tuple[Circuit,Dict[str, Dict[Qubit, int]]]
    """
    new_map: Dict[str, Dict[Qubit, int]] = {"i": {}, "o": {}}
    for qubit in circ.qubits:
        new_map["i"][qubit] = circ.qubits.index(qubit)
        new_map["o"][qubit] = circ.qubits.index(qubit)
    new_tuple = (circ.copy(), new_map)
    return new_tuple


def repeated_mbqc_conversion(
    circ: Circuit, arch: Architecture, splits: int, splitter: Splitter, knitter: Knitter
) -> Tuple[Circuit, Dict[str, Dict[Qubit, int]]]:
    """
    Takes a pytket circuit object, an architecture, a number of mbqc segments to split into,
    a splitter method and a knitter method as inputs. It splits the original circuit into
    the specified number of segments, converts each to a measurement pattern and joins them
    up in an architecture aware manner. Outputs a tuple containing the final circuit and
    its corresponding i/o map.
    
    :param circ:    The original pytket circuit.
    :param type:    Circuit
    
    :param arch:    The architecture to route onto.
    :param type:    Architecture
    
    :param splits:  Maximum number of segments to split into.
    :param type:    int
    
    :param splitter:  A splitter which determines how the circuit will be split.
    :param type:      Splitter
    
    :param knitter:   A knitter which determines how the converted segments will be joined up.
    :param type:      Knitter
    
    :returns:        A tuple containing the new circuit as well as its corresponding i/o map.
    :rtype:          Tuple[Circuit,Dict[str, Dict[Qubit, int]]]
    """
    split_circuits = splitter(circ, splits)  # type: ignore
    pattern_list = []
    for sc in split_circuits:
        if sc[1]:
            mp = MPattern(sc[0])
            pattern_list.append(mp.single_conversion())
        else:
            pattern_list.append(format_circ(sc[0]))
    final_circuit = knitter(pattern_list, arch)  # type: ignore
    return final_circuit


"""
def single_mbqc_conversion(circ: Circuit, arc: Architecture, mpattern: MPattern) -> Circuit:

    does the thing but once calls the other thing
    return repeated_mbqc_conversion(circ, arc, Splitter.SetSubcircuits(1), mpattern, Knitter.Routed)
"""
