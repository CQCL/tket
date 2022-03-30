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

from pytket.rmbqc import MPattern,Splitter,Knitter
from pytket.circuit import Circuit, Qubit
from pyzx.gflow import gflow  # type: ignore
from pytket.architecture import SquareGrid  # type: ignore

def test_zx_diagram() -> None:
    c = Circuit(4)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).CZ(1,3).T(0).X(1).T(2).T(3)
    mp = MPattern(c)
    (zx_d,io_map) = mp.zx_diagram()
    assert len(zx_d.vertices()) == 10
    assert len(zx_d.edge_set()) == 9
    assert len(zx_d.inputs) == 4
    assert len(zx_d.outputs) == 4
    
def test_label_squish() -> None:
    c = Circuit(4)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).CZ(1,3).T(0).X(1).T(2).T(3)
    mp = MPattern(c)
    (zx_d,io_map) = mp.zx_diagram()
    vn = len(zx_d.vertices())
    for v in zx_d.vertices():
        assert v < vn
    
def test_split_subgraphs() -> None:
    c = Circuit(4)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).CZ(1,3).T(0).X(1).T(2).T(3)
    mp = MPattern(c)
    (zx_d,io_map) = mp.zx_diagram()
    g_list = mp.split_subgraphs(zx_d,io_map)
    assert len(g_list) == 2
    assert len(g_list[0].vertices()) in {2,8}
    assert len(g_list[1].vertices()) in {2,8}
    assert not len(g_list[1].vertices()) == len(g_list[0].vertices()) 

def test_entangle() -> None:
    c = Circuit(4)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).CZ(1,3).T(0).X(1).T(2).T(3)
    mp = MPattern(c)
    (zx_d,io_map) = mp.zx_diagram()
    cz_c = MPattern.entangle(zx_d)
    assert cz_c.n_qubits == 10
    assert len(cz_c.get_commands()) == 15
    
def test_layer_list() -> None:
    c = Circuit(4)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).CZ(1,3).T(0).X(1).T(2).T(3)
    mp = MPattern(c)
    (zx_d,io_map) = mp.zx_diagram()
    g_list = mp.split_subgraphs(zx_d,io_map)
    gf_list = [gflow(g)[0] for g in g_list]
    assert len(MPattern.layer_list(gf_list[0])) in {2,5}
    assert len(MPattern.layer_list(gf_list[1])) in {2,5}
    assert not len(MPattern.layer_list(gf_list[0])) == len(MPattern.layer_list(gf_list[1])) 
    
def test_correct() -> None:
    c = Circuit(4)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).CZ(1,3).T(0).X(1).T(2).T(3)
    mp = MPattern(c)
    (zx_d,io_map) = mp.zx_diagram()
    g_list = mp.split_subgraphs(zx_d,io_map)
    correct_c = MPattern.correct(g_list)
    assert correct_c.n_qubits == 10
    assert len(correct_c.get_commands()) == 40

def test_single_conversion() -> None:
    c = Circuit(4)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).CZ(1,3).T(0).X(1).T(2).T(3)
    mp = MPattern(c)
    new_c,io_map = mp.single_conversion()
    assert len(io_map["i"].keys()) == 4
    assert len(io_map["o"].keys()) == 4
    assert list(io_map["i"].values()) == [5,0,1,7]
    assert list(io_map["o"].values()) == [6,9,2,7]
    
def test_depth_structure() -> None:
    c = Circuit(4)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).CZ(1,3).T(0).X(1).T(2).T(3)
    d_struct = Splitter.depth_structure(c)
    assert len(d_struct) == c.depth()
    
def test_depth_split() -> None:
    c = Circuit(4)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).CZ(1,3).T(0).X(1).T(2).T(3)
    max_n = 0
    for n in range(1,11):
        output = Splitter.depth_split(c.copy(),n)
        max_n = max(max_n,len(output))
        assert len(output) >= 1 and len(output) <= n
        total_depth = 0
        for (subcirc,some_bool) in output:
            assert some_bool
            total_depth += subcirc.depth()
        assert total_depth == c.depth()
    assert max_n == 10
        
def test_gate_split() -> None:
    c = Circuit(4)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).CZ(1,3).T(0).X(1).T(2).T(3)
    max_n = 0
    for n in range(1,11):
        output = Splitter.gates_split(c.copy(),n)
        max_n = max(max_n,len(output))
        assert len(output) >= 1 and len(output) <= n and len(output) <= 5
        total_depth = 0
        for (subcirc,some_bool) in output:
            assert some_bool
            total_depth += subcirc.depth()
        assert total_depth == c.depth()
    assert max_n == 4
    
def test_clifford_split() -> None:
    c = Circuit(4)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).CZ(1,3).T(0).X(1).T(2).T(3)
    output = Splitter.clifford_split(c.copy())
    assert len(output) == 7
    curr_bool = output[0][1]
    for (subcircuit,some_bool) in output:
        assert curr_bool == some_bool
        curr_bool = not curr_bool
        
def test_unrouted() -> None:
    c = Circuit(3)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).T(0).X(1).T(2)
    output = Splitter.depth_split(c,3)
    mp_list = [MPattern(_c).single_conversion() for (_c,bool) in output]
    (final_c,final_map) = Knitter.unrouted(mp_list)
    assert final_map["i"] == {Qubit(0): Qubit(0), Qubit(1): Qubit(3), Qubit(2): Qubit(2)}
    assert final_map["o"] == {Qubit(0): Qubit(5), Qubit(1): Qubit(3), Qubit(2): Qubit(4)}
    assert final_c.depth() == 23
    assert final_c.n_qubits == 7
    
def test_sequential() -> None:
    c = Circuit(3)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).T(0).X(1).T(2)
    output = Splitter.depth_split(c,3)
    mp_list = [MPattern(_c).single_conversion() for (_c,bool) in output]
    arch = SquareGrid(3,3)
    (final_c,final_map) = Knitter.sequential(mp_list, arch)
    assert final_map["i"] == {Qubit(0): arch.nodes[5], Qubit(1): arch.nodes[0], Qubit(2): arch.nodes[2]}
    assert final_map["o"] == {Qubit(0): arch.nodes[7], Qubit(1): arch.nodes[3], Qubit(2): arch.nodes[4]}
    assert final_c.depth() == 35
    assert final_c.n_qubits == 9
    
def test_separate() -> None:
    c = Circuit(3)
    c.H(0).T(2).CZ(0,1).CX(0,2).X(1).H(2).T(2).CZ(1,2).CX(1,0).T(0).T(1).CZ(0,2).T(0).X(1).T(2)
    output = Splitter.depth_split(c,3)
    mp_list = [MPattern(_c).single_conversion() for (_c,bool) in output]
    arch = SquareGrid(3,3)
    (final_c,final_map) = Knitter.separate(mp_list, arch)
    assert final_map["i"] == {Qubit(0): arch.nodes[5], Qubit(1): arch.nodes[0], Qubit(2): arch.nodes[2]}
    assert final_map["o"] == {Qubit(0): arch.nodes[0], Qubit(1): arch.nodes[1], Qubit(2): arch.nodes[3]}
    assert final_c.depth() == 35
    assert final_c.n_qubits == 9