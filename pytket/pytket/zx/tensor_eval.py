# Copyright Quantinuum
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

"""Collection of methods to evaluate a ZXDiagram to a tensor. This uses the
numpy tensor features, in particular the einsum evaluation and optimisations."""

import warnings
from math import cos, floor, pi, sin, sqrt
from typing import Any

import numpy as np
import sympy

from pytket.zx import (
    CliffordGen,
    DirectedGen,
    PhasedGen,
    QuantumType,
    Rewrite,
    ZXBox,
    ZXDiagram,
    ZXGen,
    ZXType,
    ZXVert,
)

try:
    import quimb.tensor as qtn  # type: ignore
except ModuleNotFoundError:
    warnings.warn(
        'Missing package for tensor evaluation of ZX diagrams. Run "pip '
        "install 'pytket[ZX]'\" to install the optional dependencies.",
        stacklevel=2,
    )


def _gen_to_tensor(gen: ZXGen, rank: int) -> np.ndarray:
    if isinstance(gen, PhasedGen):
        return _spider_to_tensor(gen, rank)
    if isinstance(gen, CliffordGen):
        return _clifford_to_tensor(gen, rank)
    if isinstance(gen, DirectedGen):
        return _dir_gen_to_tensor(gen)
    if isinstance(gen, ZXBox):
        return _tensor_from_basic_diagram(gen.diagram)
    raise ValueError(f"Cannot convert generator of type {gen.type} to a tensor")


def _spider_to_tensor(gen: PhasedGen, rank: int) -> np.ndarray:
    try:
        if gen.type == ZXType.Hbox:
            param_c = complex(gen.param)
        else:
            param = float(gen.param)
    except TypeError as e:
        # If parameter is symbolic, we cannot evaluate the tensor
        raise ValueError(
            f"Evaluation of ZXDiagram failed due to symbolic expression {gen.param}"
        ) from e
    size = pow(2, rank)
    if gen.type == ZXType.ZSpider:
        x = param / 2.0
        modval = 2.0 * (x - floor(x))
        phase = np.exp(1j * modval * pi)
        t = np.zeros(size, dtype=complex)
        t[0] = 1.0
        t[size - 1] = phase
    elif gen.type == ZXType.XSpider:
        x = param / 2.0
        modval = 2.0 * (x - floor(x))
        phase = np.exp(1j * modval * pi)
        t = np.full(size, 1.0, dtype=complex)
        constant = pow(sqrt(0.5), rank)
        for i in range(size):
            parity = (i).bit_count()
            t[i] += phase if parity % 2 == 0 else -phase
            t[i] *= constant
    elif gen.type == ZXType.Hbox:
        t = np.full(size, 1.0, dtype=complex)
        t[size - 1] = param_c
    elif gen.type == ZXType.XY:
        x = param / 2.0
        modval = 2.0 * (x - floor(x))
        phase = np.exp(-1j * modval * pi)
        t = np.zeros(size, dtype=complex)
        t[0] = sqrt(0.5)
        t[size - 1] = sqrt(0.5) * phase
    elif gen.type == ZXType.XZ:
        x = param / 2.0
        modval = x - floor(x)
        t = np.zeros(size, dtype=complex)
        t[0] = cos(modval * pi)
        t[size - 1] = sin(modval * pi)
    elif gen.type == ZXType.YZ:
        x = param / 2.0
        modval = x - floor(x)
        t = np.zeros(size, dtype=complex)
        t[0] = cos(modval * pi)
        t[size - 1] = -1j * sin(modval * pi)
    else:
        raise ValueError(
            f"Cannot convert phased generator of type {gen.type} to a tensor"
        )
    return t.reshape(tuple([2] * rank))


def _clifford_to_tensor(gen: CliffordGen, rank: int) -> np.ndarray:
    size = pow(2, rank)
    t = np.zeros(size, dtype=complex)
    if gen.type == ZXType.PX:
        t[0] = sqrt(0.5)
        t[size - 1] = -sqrt(0.5) if gen.param else sqrt(0.5)
    elif gen.type == ZXType.PY:
        t[0] = sqrt(0.5)
        t[size - 1] = 1j * sqrt(0.5) if gen.param else -1j * sqrt(0.5)
    elif gen.type == ZXType.PZ:
        if gen.param:
            t[size - 1] = 1.0
        else:
            t[0] = 1.0
    else:
        raise ValueError(
            f"Cannot convert Clifford generator of type {gen.type} to a tensor"
        )
    return t.reshape(tuple([2] * rank))


def _dir_gen_to_tensor(gen: DirectedGen) -> np.ndarray:
    if gen.type == ZXType.Triangle:
        t = np.ones((2, 2), dtype=complex)
        t[1, 0] = 0.0
        return t
    raise ValueError(
        f"Cannot convert directed generator of type {gen.type} to a tensor"
    )


_id_tensor = np.asarray([[1, 0], [0, 1]], dtype=complex)

_boundary_types = [ZXType.Input, ZXType.Output, ZXType.Open]


def _tensor_from_basic_diagram(diag: ZXDiagram) -> np.ndarray:
    try:
        scalar = complex(diag.scalar)
    except TypeError as e:
        raise ValueError(
            f"Cannot evaluate a diagram with a symbolic scalar. Given scalar: "
            f"{diag.scalar}"
        ) from e
    all_wires = diag.wires
    indices = dict(zip(all_wires, range(len(all_wires)), strict=False))
    next_index = len(all_wires)
    tensor_list: list[Any]
    tensor_list = []
    id_wires = set()
    res_indices = []
    for b in diag.get_boundary():
        # Boundaries are handled separately to get the correct order for the
        # final indices
        bw = diag.adj_wires(b)[0]
        bwi = indices[bw]
        other = diag.other_end(bw, b)
        if diag.get_zxtype(other) in _boundary_types and bw not in id_wires:
            # Two boundaries are directly connected, so insert an id tensor for
            # this boundary
            id_ind = [bwi, next_index]
            qt = qtn.Tensor(data=_id_tensor, inds=id_ind)
            tensor_list.append(qt)
            res_indices.append(next_index)
            next_index += 1
            id_wires.add(bw)
        else:
            res_indices.append(bwi)
    for v in diag.vertices:
        gen = diag.get_vertex_ZXGen(v)
        if gen.type in _boundary_types:
            # Boundaries already handled above
            continue
        v_ind = []
        for w in diag.adj_wires(v):
            v_ind.append(indices[w])
            if diag.other_end(w, v) == v:
                v_ind.append(indices[w])
        t = _gen_to_tensor(gen, len(v_ind))
        qt = qtn.Tensor(data=t, inds=v_ind)
        tensor_list.append(qt)
    net = qtn.TensorNetwork(tensor_list)
    net.full_simplify_(seq="ADCR")
    res_ten = net.contract(output_inds=res_indices, optimize="greedy")
    result: np.ndarray
    if isinstance(res_ten, qtn.Tensor):  # noqa: SIM108
        result = res_ten.data
    else:
        # Scalar
        result = np.asarray(res_ten)
    return result * scalar


def tensor_from_quantum_diagram(diag: ZXDiagram) -> np.ndarray:
    """
    Evaluates a purely quantum :py:class:`~.ZXDiagram` as a tensor. Indices of
    the resulting tensor match the order of the boundary vertices from
    :py:meth:`~.ZXDiagram.get_boundary`.

    Throws an exception if the diagram contains any non-quantum vertex or wire,
    or if it contains any symbolic parameters.
    """
    for v in diag.vertices:
        if diag.get_qtype(v) != QuantumType.Quantum:
            raise ValueError(
                "Non-quantum vertex found. tensor_from_quantum_diagram only "
                "supports diagrams consisting of only quantum components"
            )
    for w in diag.wires:
        if diag.get_wire_qtype(w) != QuantumType.Quantum:
            raise ValueError(
                "Non-quantum wire found. tensor_from_quantum_diagram only "
                "supports diagrams consisting of only quantum components"
            )
    diag_copy = ZXDiagram(diag)
    diag_copy.multiply_scalar(1 / sympy.sqrt(diag.scalar))
    Rewrite.basic_wires().apply(diag_copy)
    return _tensor_from_basic_diagram(diag_copy)


def tensor_from_mixed_diagram(diag: ZXDiagram) -> np.ndarray:
    """
    Evaluates an arbitrary :py:class:`~.ZXDiagram` as a tensor in the doubled
    picture - that is, each quantum generator is treated as a pair of conjugate
    generators, whereas a classical generator is just itself.

    The indices of the resulting tensor match the order of the boundary
    vertices from :py:meth:`~.ZXDiagram.get_boundary`, with quantum boundaries
    split into two. For example, if the boundary is ``[qb1, cb1, qb2]``, the
    indices will match ``[qb1, qb1_conj, cb1, qb2, qb2_conj]``.

    Throws an exception if the diagram contains any symbolic parameters.
    """
    expanded = diag.to_doubled_diagram()
    Rewrite.basic_wires().apply(expanded)
    return _tensor_from_basic_diagram(expanded)


def _format_tensor_as_unitary(diag: ZXDiagram, tensor: np.ndarray) -> np.ndarray:
    in_ind = []
    out_ind = []
    boundary = diag.get_boundary()
    for i in range(len(boundary)):
        if diag.get_zxtype(boundary[i]) == ZXType.Input:
            in_ind.append(i)
        else:
            out_ind.append(i)
    shape = (pow(2, len(in_ind)), pow(2, len(out_ind)))
    all_ind = in_ind + out_ind
    reshaped = np.transpose(tensor, all_ind).reshape(shape)
    return reshaped.T


def unitary_from_quantum_diagram(diag: ZXDiagram) -> np.ndarray:
    """
    Evaluates a purely quantum :py:class:`~.ZXDiagram` as a matrix describing the
    linear map from inputs to outputs. Qubits are indexed according to ILO-BE
    convention based on relative position amongst inputs/outputs in
    :py:meth`~.ZXDiagram.get_boundary`.

    Throws an exception if the diagram contains any non-quantum vertex or wire,
    or if it contains any symbolic parameters.
    """
    tensor = tensor_from_quantum_diagram(diag)
    return _format_tensor_as_unitary(diag, tensor)


def unitary_from_classical_diagram(diag: ZXDiagram) -> np.ndarray:
    """
    Evaluates a purely classical :py:class:`~.ZXDiagram` as a matrix describing
    the linear map from inputs to outputs. Bits are indexed according to the
    ILO-BE convention based on relative position amongst inputs/outputs in
    :py:meth:`~.ZXDiagram.get_boundary`. Each quantum generator is treated as a
    pair of conjugate generators.

    Throws an exception if the diagram contains any non-classical boundary, or
    if it contains any symbolic parameters.
    """
    for b in diag.get_boundary():
        if diag.get_qtype(b) != QuantumType.Classical:
            raise ValueError(
                "Non-classical boundary vertex found. "
                "unitary_from_classical_diagram only supports diagrams with "
                "only classical boundaries"
            )
    tensor = tensor_from_mixed_diagram(diag)
    return _format_tensor_as_unitary(diag, tensor)


def density_matrix_from_cptp_diagram(diag: ZXDiagram) -> np.ndarray:
    """
    Evaluates a :py:class:`~.ZXDiagram` with quantum boundaries but possibly
    mixed quantum and classical generators as a density matrix. Inputs are
    treated identically to outputs, i.e. the result is the Choi-state of the
    diagram. Qubits are indexed according to the ILO-BE convention based on the
    ordering of boundary vertices in :py:meth:`~.ZXDiagram.get_boundary`.

    Throws an exception if the diagram contains any non-quantum boundary, or if
    it contains any symbolic parameters.
    """
    for b in diag.get_boundary():
        if diag.get_qtype(b) != QuantumType.Quantum:
            raise ValueError(
                "Non-quantum boundary vertex found. "
                "density_matrix_from_cptp_diagram only supports diagrams with "
                "only quantum boundaries"
            )
    tensor = tensor_from_mixed_diagram(diag)
    n_bounds = len(diag.get_boundary())
    shape = (pow(2, n_bounds), pow(2, n_bounds))
    # diag.to_doubled_diagram() in tensor_from_mixed_diagram will alternate
    # original boundary vertices and their conjugates
    indices = [2 * i for i in range(n_bounds)] + [2 * i + 1 for i in range(n_bounds)]
    reshaped = np.transpose(tensor, indices).reshape(shape)
    return reshaped.T


def fix_boundaries_to_binary_states(
    diag: ZXDiagram, vals: dict[ZXVert, int]
) -> ZXDiagram:
    """
    Fixes (a subset of) the boundary vertices of a :py:class:`~.ZXDiagram` to
    computational basis states/post-selection.
    """
    new_diag = ZXDiagram(diag)
    b_lookup = dict(zip(diag.get_boundary(), new_diag.get_boundary(), strict=False))
    for b, val in vals.items():
        if diag.get_zxtype(b) not in _boundary_types:
            raise ValueError("Can only set states of boundary vertices")
        if val not in [0, 1]:
            raise ValueError("Can only fix boundary states to |0> and |1>.")
        new_b = b_lookup[b]
        qtype = diag.get_qtype(b)
        assert qtype is not None
        fix_b = new_diag.add_vertex(ZXType.XSpider, float(val), qtype)
        bw = new_diag.adj_wires(new_b)[0]
        adj = new_diag.other_end(bw, new_b)
        adj_p = dict(new_diag.get_wire_ends(bw))[adj]
        new_diag.add_wire(
            u=fix_b, v=adj, v_port=adj_p, type=new_diag.get_wire_type(bw), qtype=qtype
        )
        new_diag.remove_vertex(new_b)
        new_diag.multiply_scalar(0.5 if qtype == QuantumType.Quantum else sqrt(0.5))
    return new_diag


def fix_inputs_to_binary_state(diag: ZXDiagram, vals: list[int]) -> ZXDiagram:
    """
    Fixes all input vertices of a :py:class:`~.ZXDiagram` to computational basis states.
    """
    inputs = diag.get_boundary(type=ZXType.Input)
    if len(inputs) != len(vals):
        raise ValueError(
            f"Gave {len(vals)} values for {len(inputs)} inputs of ZXDiagram"
        )
    val_dict = dict(zip(inputs, vals, strict=False))
    return fix_boundaries_to_binary_states(diag, val_dict)


def fix_outputs_to_binary_state(diag: ZXDiagram, vals: list[int]) -> ZXDiagram:
    """
    Fixes all output vertices of a :py:class:`~.ZXDiagram` to computational basis
    states/post-selection.
    """
    outputs = diag.get_boundary(type=ZXType.Output)
    if len(outputs) != len(vals):
        raise ValueError(
            f"Gave {len(vals)} values for {len(outputs)} outputs of ZXDiagram"
        )
    val_dict = dict(zip(outputs, vals, strict=False))
    return fix_boundaries_to_binary_states(diag, val_dict)
