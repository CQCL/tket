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

from .splitting import Splitter
from .knitting import Knitter
from .mpattern import MPattern


def repeated_mbqc_conversion(Circuit: circ, Architecture: arc, Splitter: splitter, MPattern: mpattern, Knitter: knitter) -> Circuit:
    """
    does the thing
    
    """
    split_circuits = splitter(circ)
    mbqc_subcircuits = mpattern(split_circuits)
    final_circuit = knitter(mbqc_circuits, arc)
    return final_circuit

def single_mbqc_conversion(Circuit: circ, Architecture: arc, MPattern: mpattern) -> Circuit:
    """
    does the thing but once calls the other thing
    
    """
    return repeated_mbqc_conversion(circ, arc, Splitter.SetSubcircuits(1), mpattern, Knitter.Routed)