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

import numpy as np
from pytket.circuit import Command, Circuit, OpType  # type: ignore
from typing import List


def depth_structure(c: Circuit) -> List[List[Command]]:
    """
    Converts a pytket circuit to a list containing 'x' lists, each containing
    'y' gates, where 'x' is the depth of the circuit and 'y' is the number
    of gates acting in a given timestep. Essentially we split the circuit
    into a list of 'timeslices'.
    
    :returns:       A list containing lists of gates.
    :rtype:         List[List[Command]]
    """
    qubits = c.qubits
    current_frontiers = [0] * len(qubits)
    depth_slices: List[List[int]] = [[] for _ in range(c.depth())]
    for gate in c.get_commands():
        qubit_indices = [qubits.index(q) for q in gate.qubits]
        max_frontier = max([current_frontiers[qid] for qid in qubit_indices])
        for qid in qubit_indices:
            current_frontiers[qid] = max_frontier + 1
        depth_slices[max_frontier].append(gate)
    return depth_slices

def is_mbqc_clifford(aGate: Command) -> bool:
    """
    Check if a gate is a valid operation for splitting Circuits into Clifford components.
    
    :param aGate:    Command to check.
    :param type:     Command
    
    :returns:        The result of the check.
    :rtype:          bool
    """
    if aGate.op.type in {
        OpType.Z,
        OpType.X,
        OpType.Y,
        OpType.H,
        OpType.S,
        OpType.V,
        OpType.Sdg,
        OpType.Vdg,
        OpType.SX,
        OpType.SXdg,
        OpType.CX,
        OpType.CY,
        OpType.CZ,
        OpType.CH,
        OpType.CV,
        OpType.CVdg,
        OpType.CSX,
        OpType.CSXdg,
        OpType.CCX,
        OpType.SWAP,
        OpType.CSWAP,
        OpType.noop,
        OpType.BRIDGE,
        OpType.Reset,
    }:
        return True
    elif aGate.op.type in {OpType.Rx, OpType.Rz, OpType.Ry, OpType.CRx, OpType.CRy, OpType.CRz}:
        if aGate.op.params[0] in {0, 1 / 2, 1, 3 / 2, 2}:
            return True
    return False

def count_nCliffords(c: Circuit) -> int:
    """
    Returns number of non-Clifford gates in a circuit.
    
    :param c:        Circuit to check.
    :param type:     Circuit
    
    :returns:        The number of non-Clifford gates.
    :rtype:          int
    """
    return sum([int(not is_mbqc_clifford(cmd)) for cmd in c.get_commands()])

def is_worthwhile(
    self,
    depth_focus_only: bool = True,
    maxWidth: int = None,
    maxDepth: int = None,
    strictness=0.5,
) -> bool:
    """
    Check if a circuit is worth converting to MBQC.
    
    :param depth_focus_only:    True if we only care about depth. False if we care about both depth and width.
    :param type:                bool
    
    :param maxWidth:     Provide an upper limit to the width of the new circuit. If exceeded, disregard the new circuit entirely.
    :param type:         int
    
    :param maxDepth:     Provide an upper limit to the depth of the new circuit. If exceeded, disregard the new circuit entirely.
    :param type:         int
    
    :param strictness:   Parameter to control how strictly predicted circuits are evaluated.
    :param type:         float
    
    :returns:            Returns true if there is at least one expected circuit which meets the specifications.
    :rtype:              bool
    """

    # Numerical averages extracted from random Clifford+T circuit samples.
    # gw ~= 2.56cw + 1.1t
    # gd ~= 13.9 + 0.51gw − 0.37cw

    # hw(n) ~= (2.56cw + 1.1t/n)*1.1
    # hd(n) ~= n*(13.9 + 0.51hw(n)/1.1 - 0.37cw)

    cw = self.c.n_qubits
    cd = self.c.depth()
    t = count_nCliffords(self.c)
    n = np.array(range(1, cd))
    hw = (2.56 * cw + 1.1 * t / n) * 1.1
    hd = n * (13.9 + 0.51 * hw / 1.1 - 0.37 * cw)
    result_array = np.vstack((n, hw, hd, hw * hd))
    delete_columns = set()
    if not (maxWidth == None):
        for i in range(result_array.shape[1]):
            if result_array[1, i] * strictness > maxWidth:
                delete_columns |= {i}
    if not (maxDepth == None):
        for i in range(result_array.shape[1]):
            if result_array[2, i] * strictness > maxDepth:
                delete_columns |= {i}
    keep_columns = set(n - 1) - delete_columns
    interesting_row = 2
    if depth_focus_only:
        interesting_row = 2
    else:
        interesting_row = 3
    benchmarks = [None, cw, cd, cw * cd]
    for c in keep_columns:
        if result_array[interesting_row, c] * strictness < benchmarks[interesting_row]:
            return True
    return False
