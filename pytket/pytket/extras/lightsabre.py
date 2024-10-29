# Copyright 2019-2024 Quantinuum
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


"""Methods to route TKET circuits using Qiskit routing"""

from typing import (
    Callable,
    List,
    Tuple,
)

from qiskit import (
    QuantumCircuit,
)
from qiskit import qasm2

from qiskit.transpiler import PassManager, CouplingMap  # type: ignore
from qiskit.transpiler.preset_passmanagers.builtin_plugins import SabreLayoutPassManager  # type: ignore
from qiskit.transpiler.passmanager_config import PassManagerConfig  # type: ignore

from pytket.circuit import (
    Circuit,
    Node,
)

from pytket.qasm import circuit_to_qasm_str, circuit_from_qasm_str
from pytket.architecture import Architecture


def _architecture_to_couplingmap(architecture: Architecture) -> CouplingMap:
    """
    Converts a pytket Architecture object to a Qiskit CouplingMap object.

    :param architecture: Architecture to be converted
    """
    # we can make some assumptions from how the Architecture object is
    # originally constructed from the Qiskit CouplingMap:
    # 1) All nodes are single indexed
    # 2) All nodes are default register
    # 3) Node with index "i" corresponds to integer "i" in the original coupling map
    # We confirm assumption 1) and 2) while producing the coupling map
    coupling_map: List[Tuple[int, int]] = []
    for edge in architecture.coupling:
        assert len(edge[0].index) == 1
        assert len(edge[1].index) == 1
        assert edge[0].reg_name == "node"
        assert edge[1].reg_name == "node"
        coupling_map.append((edge[0].index[0], edge[1].index[0]))
    return CouplingMap(coupling_map)


def _gen_lightsabre_transformation(
    architecture: Architecture,
    optimisation_level: int = 2,
    seed=0,
) -> Callable[Circuit, Circuit]:
    """
    Generates a function that can be used in a Transform to make a PassPtr that
    uses LightSABRE routing.

    :param architecture: Architecture LightSABRE routes circuits to match
    :param optimisation_level: Corresponds to qiskit optmisation levels
    :param seed: LightSABRE routing is stochastic, with this parameter setting the seed
    """
    config: PassManagerConfig = PassManagerConfig(
        coupling_map=_architecture_to_couplingmap(architecture),
        routing_method="sabre",
        seed_transpiler=seed,
    )

    def lightsabre(circuit: Circuit) -> Circuit:
        sabre_pass: PassManager = SabreLayoutPassManager().pass_manager(
            config, optimisation_level=optimisation_level
        )
        c: Circuit = circuit_from_qasm_str(
            qasm2.dumps(
                sabre_pass.run(
                    QuantumCircuit.from_qasm_str(circuit_to_qasm_str(circuit))
                )
            )
        )
        c.remove_blank_wires()
        c.rename_units({q: Node(q.index[0]) for q in c.qubits})
        return c

    return lightsabre
