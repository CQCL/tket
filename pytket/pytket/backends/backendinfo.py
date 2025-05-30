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

"""BackendInfo class: additional information on Backends"""

from dataclasses import asdict, dataclass, field
from typing import Any

from pytket.architecture import Architecture, FullyConnected
from pytket.circuit import Node, OpType

_OpTypeErrs = dict[OpType, float]
_Edge = tuple[Node, Node]


def _serialize_all_node_gate_errors(
    d: dict[Node, _OpTypeErrs] | None,
) -> list[list] | None:
    if d is None:
        return None
    return [
        [n.to_list(), {ot.name: err for ot, err in errs.items()}]
        for n, errs in d.items()
    ]


def _deserialize_all_node_gate_errors(
    l: list[list] | None,
) -> dict[Node, _OpTypeErrs] | None:
    if l is None:
        return None
    return {
        Node.from_list(n): {OpType.from_name(ot): err for ot, err in errs.items()}
        for n, errs in l
    }


def _serialize_all_edge_gate_errors(d: dict[_Edge, _OpTypeErrs] | None) -> list | None:
    if d is None:
        return None
    return [
        [[n0.to_list(), n1.to_list()], {ot.name: err for ot, err in errs.items()}]
        for (n0, n1), errs in d.items()
    ]


def _deserialize_all_edge_gate_errors(
    l: list | None,
) -> dict[_Edge, _OpTypeErrs] | None:
    if l is None:
        return None
    return {
        (Node.from_list(n0), Node.from_list(n1)): {
            OpType.from_name(ot): err for ot, err in errs.items()
        }
        for (n0, n1), errs in l
    }


def _serialize_all_readout_errors(
    d: dict[Node, list[list[float]]] | None,
) -> list[list] | None:
    if d is None:
        return None
    return [[n.to_list(), errs] for n, errs in d.items()]


def _deserialize_all_readout_errors(
    l: list[list] | None,
) -> dict[Node, list[list[float]]] | None:
    if l is None:
        return None
    return {Node.from_list(n): errs for n, errs in l}


def _serialize_averaged_node_gate_errors(
    d: dict[Node, float] | None,
) -> list[list] | None:
    if d is None:
        return None
    return [[n.to_list(), err] for n, err in d.items()]


def _deserialize_averaged_node_gate_errors(
    l: list[list] | None,
) -> dict[Node, float] | None:
    if l is None:
        return None
    return {Node.from_list(n): err for n, err in l}


def _serialize_averaged_edge_gate_errors(
    d: dict[_Edge, float] | None,
) -> list[list] | None:
    if d is None:
        return None
    return [[[n0.to_list(), n1.to_list()], err] for (n0, n1), err in d.items()]


def _deserialize_averaged_edge_gate_errors(
    l: list[list] | None,
) -> dict[tuple, float] | None:
    if l is None:
        return None
    return {(Node.from_list(n0), Node.from_list(n1)): err for (n0, n1), err in l}


def _serialize_averaged_readout_errors(
    d: dict[Node, float] | None,
) -> list[list] | None:
    if d is None:
        return None
    return [[n.to_list(), err] for n, err in d.items()]


def _deserialize_averaged_readout_errors(
    l: list[list] | None,
) -> dict[Node, float] | None:
    if l is None:
        return None
    return {Node.from_list(n): err for n, err in l}


@dataclass
class BackendInfo:
    """
    Stores various properties of a Backend.

    This provides all device information useful for compilation.

    :param name: Class name of the backend.
    :param device_name: Name of the device.
    :param version: Pytket-extension version installed when creating object.
    :param architecture: Optional device connectivity.
    :param gate_set: Set of supported gate types.
    :param n_cl_reg: number of classical registers supported.
    :param supports_fast_feedforward: Flag for hardware support of fast feedforward.
    :param supports_reset: Flag for hardware support of reset operation
    :param supports_midcircuit_meas: Flag for hardware support of midcircuit
        measurement.
    :param all_node_gate_errors: Dictionary between architecture Node and error rate
        for different single qubit operations.
    :param all_edge_gate_errors: Dictionary between architecture couplings and error
        rate for different two-qubit operations.
    :param all_readout_errors: Dictionary between architecture Node and uncorrelated
        single qubit readout errors (2x2 readout probability matrix).
    :param averaged_node_gate_errors: Dictionary between architecture Node and averaged
        error rate for all single qubit operations.
    :param averaged_edge_gate_errors: Dictionary between architecture couplings and
        averaged error rate for all two-qubit operations.
    :param averaged_readout_errors: Dictionary between architecture Node and averaged
        readout errors.
    :param misc: key-value map with further provider-specific information (must be
        JSON-serializable)
    """

    # identifying information
    name: str
    device_name: str | None
    version: str
    # hardware constraints
    architecture: Architecture | FullyConnected | None
    gate_set: set[OpType]
    n_cl_reg: int | None = None
    # additional feature support
    supports_fast_feedforward: bool = False
    supports_reset: bool = False
    supports_midcircuit_measurement: bool = False
    # additional basic device characterisation information
    all_node_gate_errors: dict[Node, dict[OpType, float]] | None = None
    all_edge_gate_errors: dict[tuple[Node, Node], dict[OpType, float]] | None = None
    all_readout_errors: dict[Node, list[list[float]]] | None = None
    averaged_node_gate_errors: dict[Node, float] | None = None
    averaged_edge_gate_errors: dict[tuple[Node, Node], float] | None = None
    averaged_readout_errors: dict[Node, float] | None = None

    # miscellaneous, eg additional noise characterisation and provider-supplied
    # information
    misc: dict[str, Any] = field(default_factory=dict)

    @property
    def nodes(self) -> list[Node]:
        """
        List of device nodes of the backend. Returns empty list
        if the `architecture` field is not provided.

        :return: List of nodes.
        """
        if self.architecture is None:
            return []
        return self.architecture.nodes

    @property
    def n_nodes(self) -> int:
        """
        Number of nodes in the architecture of the device. Returns 0
        if the `architecture` field is not provided.

        :return: Number of nodes.
        """
        return len(self.nodes)

    def add_misc(self, key: str, val: Any) -> None:
        """
        Add a new entry in BackendInfo's dictionary of additional information.

        :param key: Key to store and retrieve value.
        :param val: Value to be stored.
        """
        if key in self.misc:
            raise KeyError("Attempting to add an already existing entry to misc dict")
        self.misc[key] = val

    def get_misc(self, key: str) -> Any:
        """
        Retrieve information stored in Backend's additional information store

        :param key: Key to retrieve value.

        :raises KeyError: There is no value stored with the given key.

        :return: The value stored at the given key.
        """
        return self.misc[key]

    def to_dict(self) -> dict[str, Any]:
        """
        Generate a dictionary serialized representation of BackendInfo,
        suitable for writing to JSON.

        :return: JSON serializable dictionary.
        """

        self_dict = asdict(self)
        if self_dict["architecture"] is not None:
            self_dict["architecture"] = self_dict["architecture"].to_dict()
        self_dict["gate_set"] = [op.value for op in self_dict["gate_set"]]
        self_dict["all_node_gate_errors"] = _serialize_all_node_gate_errors(
            self_dict["all_node_gate_errors"]
        )
        self_dict["all_edge_gate_errors"] = _serialize_all_edge_gate_errors(
            self_dict["all_edge_gate_errors"]
        )
        self_dict["all_readout_errors"] = _serialize_all_readout_errors(
            self_dict["all_readout_errors"]
        )
        self_dict["averaged_node_gate_errors"] = _serialize_averaged_node_gate_errors(
            self_dict["averaged_node_gate_errors"]
        )
        self_dict["averaged_edge_gate_errors"] = _serialize_averaged_edge_gate_errors(
            self_dict["averaged_edge_gate_errors"]
        )
        self_dict["averaged_readout_errors"] = _serialize_averaged_readout_errors(
            self_dict["averaged_readout_errors"]
        )
        return self_dict

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "BackendInfo":
        """
        Construct BackendInfo object from JSON serializable dictionary
        representation, as generated by BackendInfo.to_dict.

        :return: Instance of BackendInfo constructed from dictionary.
        """
        args = dict(**d)
        arch = args["architecture"]
        if arch is not None:
            if "links" in arch:
                args["architecture"] = Architecture.from_dict(args["architecture"])
            else:
                args["architecture"] = FullyConnected.from_dict(args["architecture"])
        args["gate_set"] = {OpType(op) for op in args["gate_set"]}
        args["all_node_gate_errors"] = _deserialize_all_node_gate_errors(
            args["all_node_gate_errors"]
        )
        args["all_edge_gate_errors"] = _deserialize_all_edge_gate_errors(
            args["all_edge_gate_errors"]
        )
        args["all_readout_errors"] = _deserialize_all_readout_errors(
            args["all_readout_errors"]
        )
        args["averaged_node_gate_errors"] = _deserialize_averaged_node_gate_errors(
            args["averaged_node_gate_errors"]
        )
        args["averaged_edge_gate_errors"] = _deserialize_averaged_edge_gate_errors(
            args["averaged_edge_gate_errors"]
        )
        args["averaged_readout_errors"] = _deserialize_averaged_readout_errors(
            args["averaged_readout_errors"]
        )
        return cls(**args)


def fully_connected_backendinfo(  # type: ignore
    name: str,
    device_name: str | None,
    version: str,
    n_nodes: int,
    gate_set: set[OpType],
    **kwargs,
) -> BackendInfo:
    """
    Construct a BackendInfo with a FullyConnected architecture.

    :param name: Class name of the backend.
    :param device_name: Name of the device.
    :param version: Version of the pytket Backend.
    :param n_nodes: Number of nodes of the device.
    :param gate_set: Set of supported gate types.
    :param \\**kwargs: All further arguments are passed to the BackendInfo constructor.
    """
    return BackendInfo(
        name, device_name, version, FullyConnected(n_nodes), gate_set, **kwargs
    )
