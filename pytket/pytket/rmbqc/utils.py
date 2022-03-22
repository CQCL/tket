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


def is_Clifford(aGate: Command) -> bool:
    """
    Check if a gate is Clifford.
    
    :param aGate:    Command to check.
    :param type:     Command
    
    :returns:        The result of the check.
    :rtype:          bool
    """
    if aGate.op.get_name() in {'Z','X','Y','H','S','V','Sdg','Vdg','SX','SXdg','CX','CY','CZ','CH','CV','CVdg','CSX','CSXdg','CCX','SWAP','CSWAP','noop','BRIDGE','Reset'}:
        return True
    elif aGate.op.get_name() in {'T','Tdg'}:
        return False
    elif aGate.op.get_name() in {'Rx','Rz','Ry','CRx','CRy','CRz'}:
        if aGate.op.params[0] in {0,1/2,1,3/2,2}:
            return True
    else:
        return False


    
#OTHER----------------------------------------------------------------
#These also depend on SPLITTING.is_Clifford()
def count_nCliffords(c: Circuit) -> int:
    """
    Returns number of non-Clifford gates in a circuit.
    
    :param c:        Circuit to check.
    :param type:     Circuit
    
    :returns:        The number of non-Clifford gates.
    :rtype:          int
    """
    count = 0
    for g in c.get_commands():
        if MPattern.is_Clifford(g):
            count += 1
    return count
    
def is_worthwhile(self, improveOn: str = "Depth", maxWidth: int = None, maxDepth: int = None, strictness = 0.5) -> bool:
    """
    Check if a circuit is worth converting to MBQC.
    
    :param improveOn:    The resource we want to lower - can be "Depth", "Width" or "Both".
    :param type:         str
    
    :param maxWidth:     Provide an upper limit to the width of the new circuit. If exceeded, disregard the new circuit entirely.
    :param type:         int
    
    :param maxDepth:     Provide an upper limit to the depth of the new circuit. If exceeded, disregard the new circuit entirely.
    :param type:         int
    
    :param strictness:   Parameter to control how strictly predicted circuits are evaluated.
    :param type:         float
    
    :returns:            Returns true if there is at least one expected circuit which meets the specifications.
    :rtype:              bool
    """
    
    #Numerical averages extracted from random Clifford+T circuit samples.
    #gw ~= 2.56cw + 1.1t
    #gd ~= 13.9 + 0.51gw − 0.37cw
    
    #hw(n) ~= (2.56cw + 1.1t/n)*1.1
    #hd(n) ~= n*(13.9 + 0.51hw(n)/1.1 - 0.37cw)
    
    cw = self.c.n_qubits
    cd = self.c.depth()
    t = MPattern.count_nCliffords(self.c)
    n = np.array(range(1,cd))
    hw = (2.56*cw + 1.1*t/n)*1.1
    hd = n*(13.9 + 0.51*hw/1.1 - 0.37*cw)
    result_array = np.vstack((n,hw,hd,hw*hd))
    delete_columns = set()
    if not (maxWidth == None):
        for i in range(result_array.shape[1]):
            if result_array[1,i]*strictness > maxWidth:
                delete_columns |= {i}
    if not (maxDepth == None):
        for i in range(result_array.shape[1]):
            if result_array[2,i]*strictness > maxDepth:
                delete_columns |= {i}
    keep_columns = set(n-1) - delete_columns
    interesting_row = 2
    if improveOn == "Width":
        interesting_row = 1
    if improveOn == "Depth":
        interesting_row = 2
    elif improveOn == "Both":
        interesting_row = 3
    benchmarks = [None,cw,cd,cw*cd]
    for c in keep_columns:
        if result_array[interesting_row,c]*strictness < benchmarks[interesting_row]:
            return True
    return False