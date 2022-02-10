// Copyright 2019-2022 Cambridge Quantum Computing
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


BASIC PROBLEM:

We are interested in the standard subgraph monomorphism problem, also known as the non-induced subgraph isomorphism problem:

Given graphs P, T (the pattern graph and target graph), find a subgraph monomorphism f from P to T. Thus, f is an injective function

    f : V(P) -> V(T)

such that

    e in E(P) => f(e) in E(T).

In other words, f maps the vertices of P into the vertices of T, without clashes, such that if (v1,v2) are adjacent in P, then (f(v1), f(v2)) are adjacent in T.

(It is "non-induced" because we allow (f(v1), f(v2)) to be adjacent in T even if (v1,v2) are not adjacent in P).

Our graphs are always undirected, without loops or multiple edges.

In applications to qubit routing and placement problems, P is usually the interaction graph between the logical qubits (the qubits of a quantum circuit), and T is the interaction graph of an actual quantum computer with physical qubits.

We also consider a new problem, called Weighted Subgraph Monomorphism (WSM): each pattern edge and target edge additionally has a weight w >= 0 (which is an integer). For each subgraph monomorphism f, let

    S(f) = sum_(e in E(P)) w(e).w(f(e))

We call S(f) the "total weight" or "scalar product" of a monomorphism f.

Then, the WSM problem is to find f which minimises S(f).

(It seems unlikely that this is a new problem, but we have not found this in existing papers).

The point is that we can adjust the P-edge weights to account for less "important" qubit interactions: i.e., less often used, or occurring further in the future of the circuit.

We can also adjust the T-edge weights to account for different 2-qubit gate fidelities on the hardware, and also add extra target edges to account for the fact that it is unlikely that an f exists mapping EVERY P-edge to a T-edge on the hardware device: thus we expect that we will have to apply extra SWAP gates to move the qubits around; thus we still want adjacent P-vertices to be mapped to closeby T-vertices, even if not exactly adjacent.

(So, as time progresses, the locations of the logical qubits, i.e. the physical qubits to which they have been assigned, are expected to "drift" away from the initial positions, and so we adjust the P-edge weights with some kind of "time decay", to take account of the fact that far future assignments are not really worth taking into account very much, because the qubit locations may be quite different by then).


ALGORITHMIC COMMENTS:

This may be termed version 0.1. Many improvements are planned soon!

Some of the basic ideas are similar to the Glasgow Subgraph Solver, described for example in the paper "A Parallel, Backjumping Subgraph Isomorphism Algorithm using Supplemental Graphs" by Ciaran McCreesh and Patrick Prosser.

We do a standard backtracking search, using some simple graph theory to decide if a p-vertex could be mapped to a t-vertex, and simple heuristics to guide the search.

We aim to implement parallel searching with full nogood recording soon.

The data structures (composed of std::maps, std::sets, etc.) are far from optimised at present; much fancier data structures are planned once the basic algorithms and heuristics are nailed down.
