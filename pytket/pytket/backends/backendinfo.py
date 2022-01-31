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

""" BackendInfo class: additional information on Backends """
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Set, cast, Tuple, Union

from pytket.routing import Architecture, FullyConnected  # type: ignore
from pytket.circuit import Node, OpType  # type: ignore


@dataclass
class BackendInfo:
    """
    Stores various properties of a Backend.

    This provides all device information useful for compilation.

    :param name: Class name of the backend.
    :param device_name: Name of the device.
    :param version: Pytket-extension version installed when creating object.
    :param architecture: Device connectivity.
    :param gate_set: Set of supported gate types.
    :param supports_fast_feedforward: Flag for hardware support of fast feedforward.
    :param supports_reset: Flag for hardware support of reset operation
    :param supports_midcircuit_meas: Flag for hardware support of midcircuit
        measurement.
    :param all_node_gate_errors: Dictionary between architecture Node and error rate
        for different single qubit operations.
    :param all_edge_gate_errors: Dictionary between architecture couplings and error
        rate for different two-qubit operations.
    :param all_readout_errors: Dictionary between architecture Node and uncorrelated
        single qubit readout errors.
    :param averaged_node_gate_errors: Dictionary between architecture Node and averaged
        error rate for all single qubit operations.
    :param averaged_edge_gate_errors: Dictionary between architecture couplings and
        averaged error rate for all two-qubit operations.
    :param averaged_readout_errors: Dictionary between architecture Node and averaged
        readout errors.
    :param misc: key-value map with further provider-specific information
    """

    # identifying information
    name: str
    device_name: Optional[str]
    version: str
    # hardware constraints
    architecture: Union[Architecture, FullyConnected]
    gate_set: Set[OpType]
    # additional feature support
    supports_fast_feedforward: bool = False
    supports_reset: bool = False
    supports_midcircuit_measurement: bool = False
    # additional basic device characterisation information
    all_node_gate_errors: Optional[Dict[Node, Dict[OpType, float]]] = None
    all_edge_gate_errors: Optional[Dict[Tuple[Node, Node], Dict[OpType, float]]] = None
    all_readout_errors: Optional[Dict[Node, List[List[float]]]] = None
    averaged_node_gate_errors: Optional[Dict[Node, float]] = None
    averaged_edge_gate_errors: Optional[Dict[Tuple[Node, Node], float]] = None
    averaged_readout_errors: Optional[Dict[Node, float]] = None

    # miscellaneous, eg additional noise characterisation and provider-supplied
    # information
    misc: Dict[str, Any] = field(default_factory=dict)

    @property
    def nodes(self) -> List[Node]:
        """
        List of device nodes of the backend.

        :return: List of nodes.
        :rtype: List[Node]
        """
        return cast(List[Node], self.architecture.nodes)

    @property
    def n_nodes(self) -> int:
        """
        Number of nodes in the architecture of the device.

        :return: Number of nodes.
        :rtype: int
        """
        return len(self.nodes)

    def add_misc(self, key: str, val: Any) -> None:
        """
        Add a new entry in BackendInfo's dictionary of additional information.

        :param key: Key to store and retrieve value.
        :type key: str
        :param val: Value to be stored.
        """
        if key in self.misc:
            raise KeyError("Attempting to add an already existing entry to misc dict")
        self.misc[key] = val

    def get_misc(self, key: str) -> Any:
        """
        Retrieve information stored in Backend's additional information store

        :param key: Key to retrieve value.
        :type key: str

        :raises KeyError: There is no value stored with the given key.

        :return: The value stored at the given key.
        """
        return self.misc[key]

    def to_dict(self) -> Dict[str, Any]:
        """
        Generate a dictionary serialized representation of BackendInfo,
        suitable for writing to JSON.

        :return: JSON serializable dictionary.
        :rtype: Dict[str, Any]
        """

        self_dict = asdict(self)
        self_dict["architecture"] = self_dict["architecture"].to_dict()
        self_dict["gate_set"] = [op.value for op in self_dict["gate_set"]]
        return self_dict

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "BackendInfo":
        """
        Construct BackendInfo object from JSON serializable dictionary
        representation, as generated by BackendInfo.to_dict.

        :return: Instance of BackendInfo constructed from dictionary.
        :rtype: BackendInfo
        """
        args = dict(**d)
        arch = args["architecture"]
        if "links" in arch:
            args["architecture"] = Architecture.from_dict(args["architecture"])
        else:
            args["architecture"] = FullyConnected.from_dict(args["architecture"])
        args["gate_set"] = {OpType(op) for op in args["gate_set"]}
        return cls(**args)


def fully_connected_backendinfo(  # type: ignore
    name: str,
    device_name: Optional[str],
    version: str,
    n_nodes: int,
    gate_set: Set[OpType],
    **kwargs
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
