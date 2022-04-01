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


def format_circ(circ: Circuit) -> Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]:
    """
    Takes a pytket circuit object as input and returns a tuple containing the circuit as
    well as two dictionaries, the first mapping the original input qubits to the indices of the corresponding
    qubits in the new circuit and the second mapping the orginal output qubits to the indices of the corresponding
    qubits in the new circuit.
    
    :param circ:    The pytket circuit which requires formatting.
    :param type:    Circuit
    
    :returns:        A tuple containing the circuit as well as its corresponding input/output maps.
                        These map the input/output qubits of the original circuit to their index
                        in the qubit register.
    :rtype:          Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]
    """
    new_inputs: Dict[Qubit, int] = {}
    new_outputs: Dict[Qubit, int] = {}
    for qubit in circ.qubits:
        new_inputs[qubit] = circ.qubits.index(qubit)
        new_outputs[qubit] = circ.qubits.index(qubit)
    return (circ, new_inputs, new_outputs)


def repeated_mbqc_conversion(
    circ: Circuit,
    arch: Architecture,
    splits: int,
    splitter: Splitter,
    knitter: Knitter,
    add_barriers: bool = False,
) -> Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]:
    """
    Takes a pytket circuit object, an architecture, a number of mbqc segments to split into,
    a splitter method and a knitter method as inputs. It splits the original circuit into
    the specified number of segments, converts each to a measurement pattern and joins them
    up in an architecture aware manner. Outputs a tuple containing the final circuit and
    its corresponding input/output maps. These are dictionaries mapping the original input/
    output qubits to the corresponding inputs/outputs of the new circuit.
    
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
    
    :param add_barriers:   Option to add barriers between measurement layers.
    :param type:           bool
    
    :returns:        A tuple containing the new circuit as well as its corresponding input/output maps.
                        These map the inputs/outputs of the original circuit to their corresponding qubits
                        in the new circuit.
    :rtype:          Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]
    """
    split_circuits = splitter(circ, splits)  # type: ignore
    pattern_list = []
    for sc in split_circuits:
        if sc[1]:
            mp = MPattern(sc[0])
            pattern_list.append(mp.single_conversion(add_barriers))
        else:
            pattern_list.append(format_circ(sc[0]))
    final_circuit = knitter(pattern_list, arch)  # type: ignore
    return final_circuit


def single_mbqc_conversion(
    circ: Circuit, arch: Architecture, add_barriers: bool = False
) -> Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]:
    """
    Takes a pytket circuit object and an architecture, converts the circuit to a
    measurement pattern, extracts a new circuit from the pattern and routes it on
    the architecture. Finally, it returns a tuple containing the final circuit and
    its corresponding input/output maps.
    
    :param circ:    The original pytket circuit.
    :param type:    Circuit
    
    :param arch:    The architecture to route onto.
    :param type:    Architecture
    
    :param add_barriers:   Option to add barriers between measurement layers.
    :param type:           bool
    
    :returns:       A tuple containing the new circuit as well as its corresponding input/output maps.
                        These map the input/output qubits of the original circuit to the corresponding
                        qubits in the new circuit.
    :rtype:         Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]
    """
    return repeated_mbqc_conversion(
        circ, arch, 1, Splitter.depth_split, Knitter.sequential, add_barriers
    )
