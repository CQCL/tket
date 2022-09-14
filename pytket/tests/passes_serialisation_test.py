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

import json
import pytest
from jsonschema import RefResolver, Draft7Validator, ValidationError  # type: ignore
from pathlib import Path
from typing import Any, Dict, List

from pytket.circuit import Node, Circuit, Qubit  # type: ignore
from pytket.passes import BasePass  # type: ignore
from pytket.predicates import Predicate  # type: ignore
from pytket.architecture import Architecture  # type: ignore
from pytket.placement import Placement  # type: ignore


def standard_pass_dict(content: Dict[str, Any]) -> Dict[str, Any]:
    return {"StandardPass": content, "pass_class": "StandardPass"}


def sequence_pass_dict(content: List[Dict[str, Any]]) -> Dict[str, Any]:
    return {"SequencePass": {"sequence": content}, "pass_class": "SequencePass"}


def repeat_pass_dict(content: Dict[str, Any]) -> Dict[str, Any]:
    return {"RepeatPass": {"body": content}, "pass_class": "RepeatPass"}


def repeat_until_satisfied_pass_dict(
    content: Dict[str, Any], pred: Dict[str, Any]
) -> Dict[str, Any]:
    return {
        "RepeatUntilSatisfiedPass": {"body": content, "predicate": pred},
        "pass_class": "RepeatUntilSatisfiedPass",
    }


def nonparam_pass_dict(name: str) -> Dict[str, Any]:
    return standard_pass_dict({"name": name})


def nonparam_predicate_dict(name: str) -> Dict[str, Any]:
    return {"type": name}


example_routing_config = [
    {"name": "LexiLabellingMethod"},
    {"depth": 100, "name": "LexiRouteRoutingMethod"},
]

_arch = Architecture(
    [
        [Node(0), Node(1)],
        [Node(1), Node(2)],
        [Node(2), Node(3)],
        [Node(3), Node(4)],
        [Node(4), Node(5)],
    ]
)

example_architecture = _arch.to_dict()

example_placement = Placement(_arch).to_dict()

example_qmap = [
    [Qubit(0).to_list(), Qubit(1).to_list()],
    [Qubit(1).to_list(), Qubit(0).to_list()],
]
example_2q_circuit = Circuit(2).CX(0, 1).to_dict()
example_1q_circuit = Circuit(1).X(0).to_dict()

PARAM_PREDICATES = {
    "GateSetPredicate": {
        "type": "GateSetPredicate",
        "allowed_types": ["CX", "Rx", "Rz"],
    },
    "PlacementPredicate": {
        "type": "PlacementPredicate",
        "node_set": [Node(0).to_list(), Node(1).to_list()],
    },
    "ConnectivityPredicate": {
        "type": "ConnectivityPredicate",
        "architecture": example_architecture,
    },
    "DirectednessPredicate": {
        "type": "DirectednessPredicate",
        "architecture": example_architecture,
    },
    "MaxNQubitsPredicate": {"type": "MaxNQubitsPredicate", "n_qubits": 10},
}

NONPARAM_PREDICATES = [
    "NoClassicalControlPredicate",
    "NoFastFeedforwardPredicate",
    "NoClassicalBitsPredicate",
    "NoWireSwapsPredicate",
    "MaxTwoQubitGatesPredicate",
    "CliffordCircuitPredicate",
    "DefaultRegisterPredicate",
    "NoBarriersPredicate",
    "NoMidMeasurePredicate",
    "NoSymbolsPredicate",
    "GlobalPhasedXPredicate",
    "NormalisedTK2Predicate",
]

PREDICATES = {name: nonparam_predicate_dict(name) for name in NONPARAM_PREDICATES}
PREDICATES.update(PARAM_PREDICATES)

# Parametrized passes that satisfy pass.from_dict(d).to_dict()==d
TWO_WAY_PARAM_PASSES = {
    "KAKDecomposition": standard_pass_dict(
        {
            "name": "KAKDecomposition",
            "allow_swaps": True,
            "fidelity": 0.99,
            "target_2qb_gate": "CX",
        }
    ),
    "ThreeQubitSquash": standard_pass_dict(
        {
            "name": "ThreeQubitSquash",
            "allow_swaps": True,
        }
    ),
    "FullPeepholeOptimise": standard_pass_dict(
        {
            "name": "FullPeepholeOptimise",
            "allow_swaps": True,
            "target_2qb_gate": "CX",
        }
    ),
    "ComposePhasePolyBoxes": standard_pass_dict(
        {
            "name": "ComposePhasePolyBoxes",
            "min_size": 10,
        }
    ),
    "EulerAngleReduction": standard_pass_dict(
        {
            "name": "EulerAngleReduction",
            "euler_p": "Rx",
            "euler_q": "Rz",
            "euler_strict": True,
        }
    ),
    "RoutingPass": standard_pass_dict(
        {
            "name": "RoutingPass",
            "architecture": example_architecture,
            "routing_config": example_routing_config,
        }
    ),
    "PlacementPass": standard_pass_dict(
        {
            "name": "PlacementPass",
            "placement": example_placement,
        }
    ),
    "NaivePlacementPass": standard_pass_dict(
        {
            "name": "NaivePlacementPass",
            "architecture": example_architecture,
        }
    ),
    "RenameQubitsPass": standard_pass_dict(
        {
            "name": "RenameQubitsPass",
            "qubit_map": example_qmap,
        }
    ),
    "CliffordSimp": standard_pass_dict(
        {
            "name": "CliffordSimp",
            "allow_swaps": True,
        }
    ),
    "DecomposeSwapsToCXs": standard_pass_dict(
        {
            "name": "DecomposeSwapsToCXs",
            "architecture": example_architecture,
            "directed": True,
        }
    ),
    "DecomposeSwapsToCircuit": standard_pass_dict(
        {
            "name": "DecomposeSwapsToCircuit",
            "swap_replacement": example_2q_circuit,
        }
    ),
    "OptimisePhaseGadgets": standard_pass_dict(
        {
            "name": "OptimisePhaseGadgets",
            "cx_config": "Snake",
        }
    ),
    "OptimisePairwiseGadgets": standard_pass_dict(
        {
            "name": "OptimisePairwiseGadgets",
            "cx_config": "Snake",
        }
    ),
    "PauliSimp": standard_pass_dict(
        {
            "name": "PauliSimp",
            "pauli_synth_strat": "Sets",
            "cx_config": "Snake",
        }
    ),
    "GuidedPauliSimp": standard_pass_dict(
        {
            "name": "GuidedPauliSimp",
            "pauli_synth_strat": "Sets",
            "cx_config": "Snake",
        }
    ),
    "SimplifyInitial": standard_pass_dict(
        {
            "name": "SimplifyInitial",
            "allow_classical": True,
            "create_all_qubits": True,
        }
    ),
    "SimplifyInitial with x_circuit set": standard_pass_dict(
        {
            "name": "SimplifyInitial",
            "allow_classical": True,
            "create_all_qubits": True,
            "x_circuit": example_1q_circuit,
        }
    ),
}

# non-parametrized passes that satisfy pass.from_dict(d).to_dict()==d
TWO_WAY_NONPARAM_PASSES = [
    "CommuteThroughMultis",
    "DecomposeArbitrarilyControlledGates",
    "DecomposeBoxes",
    "DecomposeMultiQubitsCX",
    "DecomposeSingleQubitsTK1",
    "PeepholeOptimise2Q",
    "RebaseTket",
    "RebaseUFR",
    "RemoveRedundancies",
    "SynthesiseHQS",
    "SynthesiseTK",
    "SynthesiseTket",
    "SynthesiseOQC",
    "SynthesiseUMD",
    "SquashTK1",
    "SquashRzPhasedX",
    "FlattenRegisters",
    "DelayMeasures",
    "ZZPhaseToRz",
    "RemoveDiscarded",
    "SimplifyMeasured",
    "RemoveBarriers",
    "DecomposeBridges",
]

TWO_WAY_PASSES = {name: nonparam_pass_dict(name) for name in TWO_WAY_NONPARAM_PASSES}
TWO_WAY_PASSES.update(TWO_WAY_PARAM_PASSES)

CUSTOM_TWO_WAY_PASSES = {
    "Simple SequencePass": sequence_pass_dict(
        [
            TWO_WAY_PASSES["ThreeQubitSquash"],
            TWO_WAY_PASSES["NaivePlacementPass"],
        ]
    ),
    "Simple RepeatPass": repeat_pass_dict(TWO_WAY_PASSES["DecomposeBoxes"]),
    "Simple RepeatUntilSatisfiedPass": repeat_until_satisfied_pass_dict(
        TWO_WAY_PASSES["DecomposeBoxes"], PREDICATES["NoBarriersPredicate"]
    ),
}

TWO_WAY_PASSES.update(CUSTOM_TWO_WAY_PASSES)


# Passes that don't satisfy pass.from_dict(d).to_dict()==d
ONE_WAY_PASSES = {
    "FullMappingPass": standard_pass_dict(
        {
            "name": "FullMappingPass",
            "architecture": example_architecture,
            "placement": example_placement,
            "routing_config": example_routing_config,
        }
    ),
    "DefaultMappingPass": standard_pass_dict(
        {
            "name": "DefaultMappingPass",
            "architecture": example_architecture,
            "delay_measures": True,
        }
    ),
    "CXMappingPass": standard_pass_dict(
        {
            "name": "CXMappingPass",
            "architecture": example_architecture,
            "placement": example_placement,
            "routing_config": example_routing_config,
            "directed": True,
            "delay_measures": True,
        }
    ),
    "PauliSquash": standard_pass_dict(
        {
            "name": "PauliSquash",
            "pauli_synth_strat": "Sets",
            "cx_config": "Snake",
        }
    ),
    "ContextSimp": standard_pass_dict(
        {
            "name": "ContextSimp",
            "allow_classical": True,
            "x_circuit": example_1q_circuit,
        }
    ),
    "DecomposeTK2": standard_pass_dict(
        {
            "name": "DecomposeTK2",
            "fidelities": {
                "CX": None,
                "ZZMax": 1.0,
            },
            "allow_swaps": True,
        }
    ),
}


# Load the json schemas for testing
# https://stackoverflow.com/a/61632081
curr_file_path = Path(__file__).resolve().parent
schema_dir = curr_file_path.parent.parent / "schemas"
with open(schema_dir / "compiler_pass_v1.json", "r") as f:
    pass_schema = json.load(f)
with open(schema_dir / "circuit_v1.json", "r") as f:
    circ_schema = json.load(f)
with open(schema_dir / "architecture_v1.json", "r") as f:
    arch_schema = json.load(f)
with open(schema_dir / "placement_v1.json", "r") as f:
    plact_schema = json.load(f)
with open(schema_dir / "predicate_v1.json", "r") as f:
    pred_schema = json.load(f)

schema_store = {
    pass_schema["$id"]: pass_schema,
    circ_schema["$id"]: circ_schema,
    arch_schema["$id"]: arch_schema,
    plact_schema["$id"]: plact_schema,
    pred_schema["$id"]: pred_schema,
}
pass_validator_resolver = RefResolver.from_schema(pass_schema, store=schema_store)
pass_validator = Draft7Validator(pass_schema, resolver=pass_validator_resolver)
predicate_validator_resolver = RefResolver.from_schema(pred_schema, store=schema_store)
predicate_validator = Draft7Validator(
    pred_schema, resolver=predicate_validator_resolver
)


def check_pass_serialisation(
    serialised_pass: Dict[str, Any], check_roundtrip: bool = True
) -> None:
    # Check the JSON is valid
    pass_validator.validate(serialised_pass)
    # Check the JSON can be deserialised
    tk_pass = BasePass.from_dict(serialised_pass)
    # Check the newly construct pass produces valid JSON
    new_serialised_pass = tk_pass.to_dict()
    pass_validator.validate(new_serialised_pass)
    # Optionally check for roundtrip
    if check_roundtrip:
        assert new_serialised_pass == serialised_pass


def check_predicate_serialisation(serialised_predicate: Dict[str, Any]) -> None:
    # Check the JSON is valid
    predicate_validator.validate(serialised_predicate)
    # Check the JSON can be deserialised
    tk_predicate = Predicate.from_dict(serialised_predicate)
    # Check the newly construct predicate produces valid JSON
    new_serialised_predicate = tk_predicate.to_dict()
    predicate_validator.validate(new_serialised_predicate)
    # Check for roundtrip
    assert new_serialised_predicate == serialised_predicate


def test_passes_roundtrip_serialisation() -> None:
    for k, p in TWO_WAY_PASSES.items():
        try:
            check_pass_serialisation(p)
        except ValidationError as e:
            print(p)
            raise ValueError(f"Pass {k} failed serialisation test.") from e


def test_passes_onw_way_serialisation() -> None:
    for k, p in ONE_WAY_PASSES.items():
        try:
            check_pass_serialisation(p, check_roundtrip=False)
        except ValidationError as e:
            raise ValueError(f"Pass {k} failed serialisation test.") from e


def test_predicate_serialisation() -> None:
    for k, p in PREDICATES.items():
        try:
            check_predicate_serialisation(p)
        except ValidationError as e:
            raise ValueError(f"Predicate {k} failed serialisation test.") from e


def test_invalid_pass_deserialisation() -> None:
    # Unknown pass
    p = standard_pass_dict(
        {
            "name": "UnknownPass",
            "architecture": example_architecture,
        }
    )
    with pytest.raises(ValidationError) as e:
        check_pass_serialisation(p)
    err_msg = "'UnknownPass' is not one of"
    assert err_msg in str(e.value)
    # Wrong enum value
    p = standard_pass_dict(
        {
            "name": "PauliSimp",
            "pauli_synth_strat": "Wrong Strat",
            "cx_config": "Snake",
        }
    )
    with pytest.raises(ValidationError) as e:
        check_pass_serialisation(p)
    err_msg = "'Wrong Strat' is not one of"
    assert err_msg in str(e.value)
    # Missing property
    p = standard_pass_dict(
        {
            "name": "PauliSimp",
            "cx_config": "Snake",
        }
    )
    with pytest.raises(ValidationError) as e:
        check_pass_serialisation(p)
    err_msg = "'pauli_synth_strat' is a required property"
    assert err_msg in str(e.value)
    # Unsupported property
    p = standard_pass_dict(
        {
            "name": "PauliSimp",
            "pauli_synth_strat": "Sets",
            "cx_config": "Snake",
            "architecture": example_architecture,
        }
    )
    with pytest.raises(ValidationError) as e:
        check_pass_serialisation(p)
    err_msg = "too many properties"
    assert err_msg in str(e.value)


def test_invalid_predicate_deserialisation() -> None:
    # Unknown predicate
    p = {
        "type": "MyPredicate",
        "architecture": example_architecture,
    }
    with pytest.raises(ValidationError) as e:
        check_predicate_serialisation(p)
    err_msg = "'MyPredicate' is not one of"
    assert err_msg in str(e.value)
    # Missing property
    p = {
        "type": "ConnectivityPredicate",
    }
    with pytest.raises(ValidationError) as e:
        check_predicate_serialisation(p)
    err_msg = "'architecture' is a required property"
    assert err_msg in str(e.value)
    # Unsupported property
    p = {
        "type": "ConnectivityPredicate",
        "architecture": example_architecture,
        "allowed_types": ["CX", "Rx", "Rz"],
    }
    with pytest.raises(ValidationError) as e:
        check_predicate_serialisation(p)
    err_msg = "too many properties"
    assert err_msg in str(e.value)
