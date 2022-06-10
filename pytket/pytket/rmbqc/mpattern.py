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

from pytket.transform import Transform  # type: ignore
from pytket.circuit import Circuit, OpType, Qubit  # type: ignore
from pytket.extensions.pyzx import tk_to_pyzx  # type: ignore
from pyzx.simplify import interior_clifford_simp  # type: ignore
from pyzx.graph.graph_s import GraphS  # type: ignore
from pyzx.circuit.graphparser import circuit_to_graph  # type: ignore
from pyzx.gflow import gflow  # type: ignore
from typing import Dict, Tuple, List, Set, Union


class MPattern:
    """
    Class with tools to convert a pytket circuit into a new pytket circuit
    with lower depth and higher width, by using Measurement-Based Quantum
    Computation (MBQC) techniques.
    """

    def __init__(self, c: Circuit) -> None:
        """
        Initialises the MPattern object by giving it a pytket circuit.
        
        :param c:       An arbitrary pytket circuit that we want converted to 
                            a measurement pattern.
        :param type:    Circuit
        """
        self.c = c

    def label_squish(
        self, g: GraphS, input_map: Dict[Qubit, int], output_map: Dict[Qubit, int]
    ) -> None:
        """
        Updates the input/output labels of the MPattern to matched a squished
        graph caused by g.copy(). The reason why we use g.copy() is because the
        easiest way to split a zx diagram into multiple disjoint zx.diagrams is
        to first identify X disjoint clusters of qubits, the copy the original
        graph into X graphs, and then for each of the copies remove all the
        vertices outside the corresponding cluster.
        
        :param g:       A pyzx graph representing a zx diagram.
        :param type:    GraphS
        
        :param input_map:  A dictionary containing the current input mappings
                            which map the input qubits of the original circuit
                            to the corresponding vertices in the ZX diagram.
        :param type:    Dict[Qubit,int]
        
        :param output_map:  A dictionary containing the current output mappings
                            which map the output qubits of the original circuit
                            to the corresponding vertices in the ZX diagram.
        :param type:    Dict[Qubit,int]
        """
        original_labels = sorted(
            list(g.vertices())
        )  # Sort the original vertices in ascending order.
        # The following code will update the i/o maps with the new labels of the
        # input/output qubits after the squish.
        for i in input_map.keys():
            for v in range(len(original_labels)):
                if input_map[i] == original_labels[v]:
                    input_map[i] = v
                    break
        for o in output_map.keys():
            for v in range(len(original_labels)):
                if output_map[o] == original_labels[v]:
                    output_map[o] = v
                    break

    def remove_redundant(
        self, g: GraphS, input_map: Dict[Qubit, int], output_map: Dict[Qubit, int]
    ) -> None:
        """
        Removes simple edges from a zx diagram by merging the connected
        vertices.
        
        :param g:       A zx diagram with some remaining simple edges we want to remove.
        :param type:    GraphS
        
        :param input_map:   A dictionary containing the current input mapping
                            to map the input qubits of the original circuit
                            to their corresponding vertices in the ZX diagram.
        :param type:    Dict[Qubit,int]
        
        :param output_map:  A dictionary containing the current output mapping
                            to map the output qubits of the original circuit
                            to their corresponding vertices in the ZX diagram.
        :param type:    Dict[Qubit,int]
        """
        simple_edges = set()
        for edge in g.edge_set():
            if g.edge_type(edge) == 1:
                simple_edges.add(edge)
        for edge in simple_edges:
            v1, v2 = edge
            g.remove_edge(edge)
            is_boundary = True
            removing_input = True
            remove_vertex = v1
            keep_vertex = v2
            if v1 in g.inputs:
                remove_vertex = v1
                keep_vertex = v2
            elif v2 in g.inputs:
                remove_vertex = v2
                keep_vertex = v1
            elif v1 in g.outputs:
                removing_input = False
            elif v2 in g.outputs:
                removing_input = False
                remove_vertex = v2
                keep_vertex = v1
            else:
                is_boundary = False
            g.set_row(keep_vertex, g.row(remove_vertex))
            g.set_qubit(keep_vertex, g.qubit(remove_vertex))
            neighbors = g.neighbors(remove_vertex) - [keep_vertex]
            for neighbor in neighbors:
                this_type = g.edge_type((remove_vertex, neighbor))
                if g.connected(keep_vertex, neighbor):
                    if g.edge_type((keep_vertex, neighbor)) == this_type:
                        g.remove_edge((keep_vertex, neighbor))
                    else:
                        g.add_edge((keep_vertex, neighbor), this_type)
                else:
                    g.add_edge((keep_vertex, neighbor), this_type)
            g.add_to_phase(keep_vertex, g.phase(remove_vertex))
            g.remove_vertex(remove_vertex)
            if is_boundary:
                if removing_input:
                    g.inputs.remove(remove_vertex)
                    g.inputs.append(keep_vertex)
                    for i in input_map.keys():
                        if input_map[i] == remove_vertex:
                            input_map[i] = keep_vertex
                            break
                else:
                    g.outputs.remove(remove_vertex)
                    g.outputs.append(keep_vertex)
                    for o in output_map.keys():
                        if output_map[o] == remove_vertex:
                            output_map[o] = keep_vertex
                            break
                g.set_type(keep_vertex, 0)
        self.identity_cleanup(g, input_map, output_map)

    def identity_cleanup(
        self, g: GraphS, input_map: Dict[Qubit, int], output_map: Dict[Qubit, int]
    ) -> None:
        """
        Removes identity vertices from a zx diagram if any exist.
        
        :param g:       A zx diagram with a few possible remaining identity vertices.
        :param type:    GraphS
        
        :param input_map:  A dictionary mapping the input qubits of the original
                            circuit to their corresponding vertices in the ZX
                            diagram.
        :param type:    Dict[Qubit,int]
        
        :param output_map:  A dictionary mapping the output qubits of the original
                            circuit to their corresponding vertices in the ZX
                            diagram.
        :param type:    Dict[Qubit,int]
        """
        v_list = []
        for v in g.vertices():
            if (
                (len(g.neighbors(v)) == 2)
                and (g.phase(v) == 0)
                and not (v in g.inputs + g.outputs)
            ):
                v_list.append(v)
                neigh = g.neighbors(v)
                new_edge = (neigh[0], neigh[1])
                g.add_edge(new_edge, 1)
        for v in v_list:
            g.remove_vertex(v)
        if len(v_list) > 0:
            self.remove_redundant(g, input_map, output_map)

    def zx_diagram(self) -> Tuple[GraphS, Dict[Qubit, int], Dict[Qubit, int]]:
        """
        Converts a pytket circuit to a zx diagram.
        
        :returns:       A tuple containing the zx diagram and new input/output maps
                            which map the input/output qubits of the original circuit
                            to their corresponding vertices in the zx diagram.
        :rtype:         Tuple[GraphS, Dict[Qubit,int], Dict[Qubit,int]]
        """
        rebased_c = (
            self.c.copy()
        )  # Creates copy of original circuit on which we will work.
        rebased_c.flatten_registers()  # Pyzx can't handle pytket qubit labelling so we have to flatten the register.
        Transform.RebaseToPyZX().apply(
            rebased_c
        )  # Rebasing the gates into something pyzx can read.
        pyzxc = tk_to_pyzx(rebased_c)  # Convert to pyzx circuit.
        g = circuit_to_graph(pyzxc)  # Get graph (zx diagram) from pyzx circuit.
        interior_clifford_simp(g, quiet=True)  # Remove interior clifford spiders.
        new_inputs = {}
        new_outputs = {}
        sorted_vertices = sorted(
            list(g.vertices())
        )  # Probably already sorted but don't want to take the chance.
        for q in range(self.c.n_qubits):
            new_inputs[self.c.qubits[q]] = sorted_vertices[q]
            new_outputs[self.c.qubits[q]] = sorted_vertices[q - self.c.n_qubits]
        self.remove_redundant(
            g, new_inputs, new_outputs
        )  # Removes the last few remaining simple edges in the diagram.
        # We assume that g.copy() will squash the vertex
        # labels and thus we keep track of the new input/output vertices. If
        # pyzx is updated such that graph.copy() no longer changes vertex labels
        # then comment out the next line (label_squish(g)).
        self.label_squish(
            g, new_inputs, new_outputs
        )  # Squishes the labels of the graph vertices to fill in empty vertex indices.
        return (g.copy(), new_inputs, new_outputs)

    @staticmethod
    def entangle(g: GraphS, add_barriers: bool = False) -> Circuit:
        """
        Creates a pytket circuit which implements the edges of a zx diagram
        via CZ gates on pairs of qubits.
        
        :param g:       A zx diagram the edges of which we want to implement.
        :param type:    GraphS
        
        :param add_barriers: Set to True to add barrier at the end of the circuit.
        :param type:         bool
        
        :returns:       A pytket circuit.
        :rtype:         Circuit
        """
        vlist = list(g.vertices())
        vnum = len(vlist)
        c = Circuit(vnum, vnum)
        dlist = [
            len(g.neighbors(v)) for v in vlist
        ]  # Tracks number of remaining edges for each vertex.
        for v in vlist:
            if not v in g.inputs:
                c.H(v)
        sorted_vlist = [
            x for _, x in sorted(zip(dlist, vlist))
        ]  # Sort the vertices by descending number of edges.
        sorted_vlist.reverse()
        edge_pool = set(g.edge_set())
        finished_edge_pool: Dict[int, int] = {}
        doneround = [False] * len(
            sorted_vlist
        )  # Which qubits have had CZ gates placed on them during this round of application.
        while (
            len(edge_pool) > 0
        ):  # Do this while there are still un-implemented edges pending.
            for vid1 in range(len(sorted_vlist)):  # Look through all the vertices.
                if not doneround[
                    vid1
                ]:  # Only interested in the ones that haven't been involved in a CZ during the current timestep.
                    for vid2 in range(
                        vid1 + 1, len(sorted_vlist)
                    ):  # Look for another vertex which has a pending edge with the current vertex and has also not been involved yet.
                        if not doneround[vid2]:
                            if (
                                (sorted_vlist[vid1], sorted_vlist[vid2]) in edge_pool
                            ) or (
                                (sorted_vlist[vid2], sorted_vlist[vid1]) in edge_pool
                            ):
                                c.CZ(
                                    sorted_vlist[vid1], sorted_vlist[vid2]
                                )  # Implement the CZ gate corresponding to the edge.
                                doneround[vid1] = True
                                doneround[vid2] = True
                                edge_pool -= set(
                                    [(sorted_vlist[vid1], sorted_vlist[vid2])]
                                )  # Remove the edge from the edge pool.
                                edge_pool -= set(
                                    [(sorted_vlist[vid2], sorted_vlist[vid1])]
                                )
                                # Following if statements increment the values in a dictionary which correspond to the the number of implemented edges for each vertex.
                                if sorted_vlist[vid1] in finished_edge_pool.keys():
                                    finished_edge_pool[sorted_vlist[vid1]] += 1
                                else:
                                    finished_edge_pool[sorted_vlist[vid1]] = 1
                                if sorted_vlist[vid2] in finished_edge_pool.keys():
                                    finished_edge_pool[sorted_vlist[vid2]] += 1
                                else:
                                    finished_edge_pool[sorted_vlist[vid2]] = 1
                                break
                        else:
                            continue
                else:
                    continue
            for v in range(len(doneround)):
                doneround[v] = False
            dlist = []
            # After every timestep we need to sort the list again to make sure that the vertices with the highest vertex degree are always prioritized.
            for v in vlist:
                if v in finished_edge_pool.keys():
                    dlist.append(
                        len(g.neighbors(v)) - finished_edge_pool[v]
                    )  # Remaining vertex degree = original vertex degree minus the number of edges already implemented.
                else:
                    dlist.append(len(g.neighbors(v)))
                    finished_edge_pool[v] = 0
            sorted_vlist = [x for _, x in sorted(zip(dlist, vlist))]
            sorted_vlist.reverse()
        if add_barriers:
            active_vertices = []
            for v in g.vertices():
                if len(g.neighbors(v)) > 0:
                    active_vertices.append(v)
            c.add_barrier(active_vertices, active_vertices)
        return c

    def split_subgraphs(
        self, g: GraphS, input_map: Dict[Qubit, int], output_map: Dict[Qubit, int]
    ) -> List[GraphS]:
        """
        If a zx diagram contains sub-diagrams which are not connected to each
        other, it splits them into multiple zx diagrams. It returns a list of
        all the irreducible zx diagrams contained by the original. The easiest
        way to do this is to first identify all the disjoint clusters of vertices
        in the original diagram, then produce a copy of the original for each
        cluster and remove all the vertices that do not belong to that cluster
        from the corresponding copy. However, as a result of using graph.copy()
        the labelling in the zx diagram changes which is why label_squish is
        needed.
        
        :param g:       A zx diagram which may contain disjointed sub-diagrams.
        :param type:    GraphS
        
        :param input_map:  A dictionary containing the current input map, which
                            maps the input qubits of the original circuit to
                            their corresponding vertices in the ZX diagram.
        :param type:    Dict[Qubit, int]
        
        :param output_map:  A dictionary containing the current output map, which
                             maps the output qubits of the original circuit to
                             their corresponding vertices in the ZX diagram.
        :param type:    Dict[Qubit, int]
        
        :returns:       A list of zx diagrams.
        :rtype:         List[GraphS]
        """
        #'label_squish()' is ran before 'g.copy()' to keep track of input/
        # output qubit labels.
        self.label_squish(
            g, input_map, output_map
        )  # Must squish labels again because we are going to use graph.copy()
        g1 = g.copy()  # Create copy of the graph to work on.
        cluster_list: List[
            Set[int]
        ] = []  # Will contain all the sets of interconnected vertices.
        for v in g1.vertices():
            found = False
            for cluster in cluster_list:
                if v in cluster:
                    found = (
                        True  # Adds a flag if the current vertex exists in any cluster.
                    )
            if (
                not found
            ):  # If the current vertex isn't in a cluster, create new cluster.
                new_set = set([v])  # Current state of new cluster.
                new_nodes = set([v])  # The latest additions to the new cluster.
                while (
                    len(new_nodes) > 0
                ):  # If there are no new nodes to add to the cluster the cluster is complete.
                    temp = set()
                    for v2 in new_nodes:
                        temp |= set(g1.neighbors(v2))  # Get neighbors of new additions.
                    new_nodes = (
                        temp - new_set
                    )  # If they have already been added to the cluster they are not new additions.
                    new_set |= new_nodes  # Add new additions to the cluster.
                cluster_list.append(new_set)  # Add new cluster to the list.
        graph_list = []  # This is a list of all the disjoint subgraphs.
        for cl in range(
            len(cluster_list)
        ):  # We will extract one subgraph for each cluster.
            curr_cluster = cluster_list[cl]
            new_g = g1.copy()  # The subgraph can start as a copy of the full graph.
            for v in set(
                new_g.vertices()
            ):  # Remove each vertex not in the current cluster from the current subgraph.
                if v not in curr_cluster:
                    new_g.remove_edges(new_g.incident_edges(v))
                    if v in new_g.inputs:
                        new_g.inputs.remove(v)
                    if v in new_g.outputs:
                        new_g.outputs.remove(v)
                    new_g.remove_vertex(v)
            graph_list.append(new_g)  # Add new subgraph to the list.
        return graph_list

    @staticmethod
    def layer_list(layers: Dict[int, int]) -> List[Set[int]]:
        """
        This method takes a dictionary which maps each vertex in a zx diagram
        to a correction layer as input. It then produces a list of sets of 
        vertices as output. The sets in the list represent the layers of
        measurements for the diagram in ascending order.
        
        :param layers:  Dictionary mapping vertices to their correction layer.
        :param type:    Dict[int,int]
        
        :returns:       A list of sets containing integers.
        :rtype:         List[Set[int]]
        """
        new_list: List[Set[int]] = []
        for l in range(max(list(layers.values())) + 1):
            new_list.append(set())
        for vertex in layers.keys():
            new_list[layers[vertex]] |= {vertex}
        return new_list

    @staticmethod
    def correct(glist: List[GraphS], add_barriers: bool = False) -> Circuit:
        """
        This method takes a list of subgraphs as input and produces a circuit
        of measurements and corrections which ensures that the underlying
        graphs are implemented deterministically.
        
        :param glist:   A list of unconnected graphs.
        :param type:    List[GraphS]
        
        :param add_barriers:   Option to add barriers between measurement layers.
        :param type:           bool
        
        :returns:       A circuit containing measurements and conditional gates.
        :rtype:         Circuit
        """
        signals: Dict[int, Dict[str, Set[Union[int, bool]]]] = {}
        total_v = 0
        for g in glist:
            total_v += len(g.vertices())
            for v in g.vertices():
                signals[v] = {"x": set(), "z": set()}
        new_c = Circuit(total_v, total_v)
        for g in glist:
            gf = gflow(g)
            if gflow == None:
                raise TypeError(
                    "This graph doesn't have gflow (which shouldn't be the case if it's based on a circuit)."
                )
            else:
                l_list = MPattern.layer_list(gf[0])
                reset_list: List[int] = []
                for corr_layer in range(len(l_list) - 1):
                    if corr_layer > 0:
                        isClifford = True
                        for v in l_list[
                            -1 - corr_layer
                        ]:  # If all the measurements in the next layer are Clifford then it can be performed in parallel to the current layer.
                            if g.phase(v) not in {
                                0,
                                1 / 2,
                                1,
                                3 / 2,
                            }:  # Check if the phase of the measurement is non-Clifford.
                                isClifford = False
                                break
                        if not isClifford:
                            if add_barriers:
                                new_c.add_barrier(
                                    list(g.vertices()), list(g.vertices())
                                )
                            for v in reset_list:
                                new_c.add_gate(OpType.Reset, [v])
                            reset_list = []
                    for v in l_list[
                        -1 - corr_layer
                    ]:  # Iterate through the layers in reverse order (layer 0 are the outputs and layer 1 is the last layer to be measured)
                        my_result = {
                            v
                        }  # Measurement of qubit 'v' is initially assumed to only depend on itself.
                        if g.phase(v) in {0, 1 / 2, 1, 3 / 2}:
                            my_result ^= signals[v][
                                "z"
                            ]  # Depending on the phase of qubit v its outcome will be affected by Z corrections.
                        if g.phase(v) in {1 / 2, 1, 3 / 2}:
                            my_result ^= signals[v][
                                "x"
                            ]  # Depending on the phase of qubit v its outcome will be affected by X corrections.
                        if g.phase(v) in {1 / 2, 3 / 2}:
                            my_result ^= {
                                True
                            }  # Depending on the phase of qubit v its outcome may be flipped.
                        for u in gf[1][v] - {
                            v
                        }:  # Adds corrections to other qubits which carry over the dependencies of current qubit.
                            signals[u]["x"] ^= my_result
                        for u in g.vertices() - {v}:
                            Nu = set(g.neighbors(u))
                            if (len(Nu & gf[1][v]) % 2) == 1:
                                signals[u][
                                    "z"
                                ] ^= my_result  # Adds corections to other qubits which carry over the dependencies of current qubit.
                        # Measurement plane/angle of current qubit is affected by its phase and dependencies.
                        # For more details on how these measurements and corrections work refer to
                        # https://arxiv.org/abs/quant-ph/0702212
                        if g.phase(v) in {0, 1}:
                            new_c.H(v)
                        elif g.phase(v) in {1 / 2, 3 / 2}:
                            new_c.Rx(-g.phase(v), v)
                        else:
                            new_c.H(v)
                            # theta = zi-(((-1)**xi)*g.phase(v))
                            zi = False
                            for val in signals[v]["z"]:
                                if type(val) == bool:
                                    zi ^= val
                                else:
                                    zi ^= new_c.bits[val]
                            xi = False
                            for val in signals[v]["x"]:
                                if type(val) == bool:
                                    xi ^= val
                                else:
                                    xi ^= new_c.bits[val]
                            if type(zi) == bool:
                                if zi:
                                    new_c.X(v)
                            else:
                                new_c.X(v, condition=zi)
                            if type(xi) == bool:
                                if xi:
                                    new_c.Rx(g.phase(v), v)
                                else:
                                    new_c.Rx(-g.phase(v), v)
                            else:
                                new_c.Rx(g.phase(v), v, condition=xi)
                                new_c.Rx(-g.phase(v), v, condition=(xi ^ True))
                        new_c.Measure(v, v)
                        reset_list.append(v)
                if add_barriers and len(l_list) > 1:
                    new_c.add_barrier(list(g.vertices()), list(g.vertices()))
                for v in reset_list:
                    new_c.add_gate(OpType.Reset, [v])
                for v in l_list[0]:
                    zi = False
                    for val in signals[v]["z"]:
                        if type(val) == bool:
                            zi ^= val
                        else:
                            zi ^= new_c.bits[val]
                    xi = False
                    for val in signals[v]["x"]:
                        if type(val) == bool:
                            xi ^= val
                        else:
                            xi ^= new_c.bits[val]
                    if type(xi) == bool:
                        if xi:
                            new_c.X(v)
                    else:
                        new_c.X(v, condition=xi)
                    if type(zi) == bool:
                        if zi:
                            new_c.Z(v)
                    else:
                        new_c.Z(v, condition=zi)
                    if g.phase(v) == 1:
                        new_c.Z(v)
                    elif g.phase(v) == 0:
                        pass
                    else:
                        new_c.Rz(-g.phase(v), v)
        return new_c

    def single_conversion(
        self, add_barriers: bool = False
    ) -> Tuple[Circuit, Dict[Qubit, int], Dict[Qubit, int]]:
        """
        Converts a pytket circuit to another with reduced depth and higher width.
        
        :param add_barriers:   Option to add barriers between measurement layers.
        :param type:           bool
        
        :returns:       A tuple containing the new circuit and the input/output maps
                            which are dictionaries mapping the inputs and outputs
                            of the original circuit to their new locations in the
                            qubit register.
        :rtype:         Tuple[Circuit, Dict[Qubit,int], Dict[Qubit,int]]
        """
        (
            g,
            input_map,
            output_map,
        ) = self.zx_diagram()  # Creates a simplified ZX diagram.
        subs = self.split_subgraphs(
            g, input_map, output_map
        )  # Returns list of disjoint subgraphs.
        cz_circ = self.entangle(
            g, add_barriers
        )  # Creates the CZ circuit part of the pattern.
        m_circ = self.correct(
            subs, add_barriers
        )  # Circuit implementing measurements/corrections.
        cz_circ.add_circuit(m_circ, [])
        return (cz_circ, input_map, output_map)
