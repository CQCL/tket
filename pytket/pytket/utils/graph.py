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

from collections import defaultdict
from itertools import combinations
from tempfile import NamedTemporaryFile
from typing import Optional

import networkx as nx  # type: ignore
import graphviz as gv  # type: ignore

from pytket.circuit import Circuit


class Graph:
    def __init__(self, c: Circuit):
        """
        A class for visualising a circuit as a directed acyclic graph (DAG).

        Note: in order to use graph-rendering methods, such as
        :py:meth:`Graph.save_DAG`, it is necessary to have the Graphviz tools installed
        and on your path. See the `Graphviz website <https://graphviz.org/download/>`_
        for instructions on how to install them.

        :param      c:    Circuit
        :type       c:    pytket.Circuit
        """
        (
            q_inputs,
            c_inputs,
            w_inputs,
            q_outputs,
            c_outputs,
            w_outputs,
            input_names,
            output_names,
            node_data,
            edge_data,
        ) = c._dag_data
        self.q_inputs = q_inputs
        self.c_inputs = c_inputs
        self.w_inputs = w_inputs
        self.q_outputs = q_outputs
        self.c_outputs = c_outputs
        self.w_outputs = w_outputs
        self.input_names = input_names
        self.output_names = output_names
        self.node_data = node_data
        self.Gnx: Optional[nx.MultiDiGraph] = None
        self.G: Optional[gv.Digraph] = None
        self.Gqc: Optional[gv.Graph] = None
        self.edge_data: dict[tuple[int, int], list[tuple[int, int, str]]] = defaultdict(
            list
        )
        self.port_counts: dict = defaultdict(int)
        for src_node, tgt_node, src_port, tgt_port, edge_type in edge_data:
            self.edge_data[(src_node, tgt_node)].append((src_port, tgt_port, edge_type))
            self.port_counts[(src_node, src_port)] += 1

    def as_nx(self) -> nx.MultiDiGraph:
        """
        Return a logical representation of the circuit as a DAG.

        :returns:   Representation of the DAG
        :rtype:     networkx.MultiDiGraph
        """
        if self.Gnx is not None:
            return self.Gnx
        Gnx = nx.MultiDiGraph()
        for node, desc in self.node_data.items():
            Gnx.add_node(node, desc=desc)
        for nodepair, portpairlist in self.edge_data.items():
            src_node, tgt_node = nodepair
            for src_port, tgt_port, edge_type in portpairlist:
                Gnx.add_edge(
                    src_node,
                    tgt_node,
                    src_port=src_port,
                    tgt_port=tgt_port,
                    edge_type=edge_type,
                )

        # Add node IDs to edges
        for edge in nx.topological_sort(nx.line_graph(Gnx)):
            src_node, tgt_node, _ = edge
            # List parent edges with matching port number
            src_port = Gnx.edges[edge]["src_port"]
            prev_edges = [
                e
                for e in Gnx.in_edges(src_node, keys=True)
                if Gnx.edges[e]["tgt_port"] == src_port
            ]
            if not prev_edges:
                # The source must be an input node
                unit_id = src_node
                nx.set_edge_attributes(Gnx, {edge: {"unit_id": unit_id}})
            else:
                # The parent must be unique
                assert len(prev_edges) == 1
                prev_edge = prev_edges[0]
                unit_id = Gnx.edges[prev_edge]["unit_id"]
                nx.set_edge_attributes(Gnx, {edge: {"unit_id": unit_id}})

        # Remove unnecessary port attributes to avoid clutter:
        for node in Gnx.nodes:
            if Gnx.in_degree(node) == 1:
                for edge in Gnx.in_edges(node, keys=True):
                    nx.set_edge_attributes(Gnx, {edge: {"tgt_port": None}})
                for edge in Gnx.out_edges(node, keys=True):
                    nx.set_edge_attributes(Gnx, {edge: {"src_port": None}})

        self.Gnx = Gnx
        return Gnx

    def get_DAG(self) -> gv.Digraph:
        """
        Return a visual representation of the DAG as a graphviz object.

        :returns:   Representation of the DAG
        :rtype:     graphviz.DiGraph
        """
        if self.G is not None:
            return self.G
        G = gv.Digraph(
            "Circuit",
            strict=True,
        )
        G.attr(rankdir="LR", ranksep="0.3", nodesep="0.15", margin="0")
        q_color = "blue"
        c_color = "slategray"
        b_color = "gray"
        w_color = "green"
        gate_color = "lightblue"
        boundary_cluster_attr = {
            "style": "rounded, filled",
            "color": "lightgrey",
            "margin": "5",
        }
        boundary_node_attr = {"fontname": "Courier", "fontsize": "8"}
        with G.subgraph(name="cluster_q_inputs") as c:
            c.attr(rank="source", **boundary_cluster_attr)
            c.node_attr.update(shape="point", color=q_color)
            for node in self.q_inputs:
                c.node(
                    str((node, 0)), xlabel=self.input_names[node], **boundary_node_attr
                )
        with G.subgraph(name="cluster_c_inputs") as c:
            c.attr(rank="source", **boundary_cluster_attr)
            c.node_attr.update(shape="point", color=c_color)
            for node in self.c_inputs:
                c.node(
                    str((node, 0)), xlabel=self.input_names[node], **boundary_node_attr
                )
        with G.subgraph(name="cluster_w_inputs") as c:
            c.attr(rank="source", **boundary_cluster_attr)
            c.node_attr.update(shape="point", color=w_color)
            for node in self.w_inputs:
                c.node(
                    str((node, 0)), xlabel=self.input_names[node], **boundary_node_attr
                )
        with G.subgraph(name="cluster_q_outputs") as c:
            c.attr(rank="sink", **boundary_cluster_attr)
            c.node_attr.update(shape="point", color=q_color)
            for node in self.q_outputs:
                c.node(
                    str((node, 0)), xlabel=self.output_names[node], **boundary_node_attr
                )
        with G.subgraph(name="cluster_c_outputs") as c:
            c.attr(rank="sink", **boundary_cluster_attr)
            c.node_attr.update(shape="point", color=c_color)
            for node in self.c_outputs:
                c.node(
                    str((node, 0)), xlabel=self.output_names[node], **boundary_node_attr
                )
        with G.subgraph(name="cluster_w_outputs") as c:
            c.attr(rank="sink", **boundary_cluster_attr)
            c.node_attr.update(shape="point", color=w_color)
            for node in self.w_outputs:
                c.node(
                    str((node, 0)), xlabel=self.output_names[node], **boundary_node_attr
                )
        boundary_nodes = (
            self.q_inputs
            | self.c_inputs
            | self.w_inputs
            | self.q_outputs
            | self.c_outputs
            | self.w_outputs
        )
        Gnx = self.as_nx()
        node_cluster_attr = {
            "style": "rounded, filled",
            "color": gate_color,
            "fontname": "Times-Roman",
            "fontsize": "10",
            "margin": "5",
            "lheight": "100",
        }
        port_node_attr = {
            "shape": "point",
            "weight": "2",
            "fontname": "Helvetica",
            "fontsize": "8",
        }
        for node, ndata in Gnx.nodes.items():
            if node not in boundary_nodes:
                with G.subgraph(name="cluster_" + str(node)) as c:
                    c.attr(label=ndata["desc"], **node_cluster_attr)
                    n_ports = Gnx.in_degree(node)
                    if n_ports == 1:
                        c.node(name=str((node, 0)), **port_node_attr)
                    else:
                        for i in range(n_ports):
                            c.node(name=str((node, i)), xlabel=str(i), **port_node_attr)
        edge_colors = {
            "Quantum": q_color,
            "Boolean": b_color,
            "Classical": c_color,
            "WASM": w_color,
        }
        edge_attr = {
            "weight": "2",
            "arrowhead": "vee",
            "arrowsize": "0.2",
            "headclip": "true",
            "tailclip": "true",
        }
        for edge, edata in Gnx.edges.items():
            src_node, tgt_node, _ = edge
            src_port = edata["src_port"] or 0
            tgt_port = edata["tgt_port"] or 0
            edge_type = edata["edge_type"]
            src_nodename = str((src_node, src_port))
            tgt_nodename = str((tgt_node, tgt_port))
            G.edge(
                src_nodename, tgt_nodename, color=edge_colors[edge_type], **edge_attr
            )
        self.G = G
        return G

    def save_DAG(self, name: str, fmt: str = "pdf") -> None:
        """
        Save an image of the DAG to a file.

        The actual filename will be "<name>.<fmt>". A wide range of formats is
        supported. See https://graphviz.org/doc/info/output.html.

        :param      name:   Prefix of file name
        :type       name:   str
        :param      fmt:    File format, e.g. "pdf", "png", ...
        :type       fmt:    str
        """
        G = self.get_DAG()
        G.render(name, cleanup=True, format=fmt, quiet=True)

    def view_DAG(self) -> str:
        """
        View the DAG.

        This method creates a temporary file, and returns its filename so that the
        caller may delete it afterwards.

        :returns: filename of temporary created file
        """
        G = self.get_DAG()
        filename = NamedTemporaryFile(delete=False).name
        G.view(filename, quiet=True)
        return filename

    def get_qubit_graph(self) -> gv.Graph:
        """
        Return a visual representation of the qubit connectivity graph as a graphviz
        object.

        :returns:   Representation of the qubit connectivity graph of the circuit
        :rtype:     graphviz.Graph
        """
        if self.Gqc is not None:
            return self.Gqc
        Gnx = self.as_nx()
        Gqcnx = nx.Graph()
        for node in Gnx.nodes():
            qubits = []
            for e in Gnx.in_edges(node, keys=True):
                unit_id = Gnx.edges[e]["unit_id"]
                if unit_id in self.q_inputs:
                    qubits.append(unit_id)

            Gqcnx.add_edges_from(combinations(qubits, 2))
        G = gv.Graph(
            "Qubit connectivity",
            node_attr={
                "shape": "circle",
                "color": "blue",
                "fontname": "Courier",
                "fontsize": "10",
            },
            engine="neato",
        )
        G.edges(
            (self.input_names[src], self.input_names[tgt]) for src, tgt in Gqcnx.edges()
        )
        self.Gqc = G
        return G

    def view_qubit_graph(self) -> str:
        """
        View the qubit connectivity graph.

        This method creates a temporary file, and returns its filename so that the
        caller may delete it afterwards.

        :returns: filename of temporary created file
        """
        G = self.get_qubit_graph()
        filename = NamedTemporaryFile(delete=False).name
        G.view(filename, quiet=True)
        return filename

    def save_qubit_graph(self, name: str, fmt: str = "pdf") -> None:
        """
        Save an image of the qubit connectivity graph to a file.

        The actual filename will be "<name>.<fmt>". A wide range of formats is
        supported. See https://graphviz.org/doc/info/output.html.

        :param      name:   Prefix of file name
        :type       name:   str
        :param      fmt:    File format, e.g. "pdf", "png", ...
        :type       fmt:    str
        """
        G = self.get_qubit_graph()
        G.render(name, cleanup=True, format=fmt, quiet=True)
