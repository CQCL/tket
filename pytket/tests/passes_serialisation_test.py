# Copyright 2019-2023 Cambridge Quantum Computing
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

from pytket.circuit import Node, Circuit, Qubit, OpType  # type: ignore
from pytket.predicates import Predicate  # type: ignore
from pytket.architecture import Architecture  # type: ignore
from pytket.placement import Placement, GraphPlacement  # type: ignore
from pytket._tket.circuit import _library  # type: ignore

from pytket.passes import (  # type: ignore
    BasePass,
    SequencePass,
    RemoveRedundancies,
    RepeatUntilSatisfiedPass,
    CommuteThroughMultis,
    RepeatWithMetricPass,
    RebaseCustom,
    CXMappingPass,
    FullMappingPass,
    DefaultMappingPass,
    AASRouting,
    SquashCustom,
    RoundAngles,
)
from pytket.mapping import (  # type: ignore
    LexiLabellingMethod,
    LexiRouteRoutingMethod,
    MultiGateReorderRoutingMethod,
    BoxDecompositionRoutingMethod,
)


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
    "CommutableMeasuresPredicate",
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
    "PeepholeOptimise2Q": standard_pass_dict(
        {
            "name": "PeepholeOptimise2Q",
            "allow_swaps": True,
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
    "PauliExponentials": standard_pass_dict(
        {
            "name": "PauliExponentials",
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
    "DelayMeasures": standard_pass_dict(
        {
            "name": "DelayMeasures",
            "allow_partial": False,
        }
    ),
    "RoundAngles": standard_pass_dict(
        {"name": "RoundAngles", "n": 6, "only_zeros": False}
    ),
}

# non-parametrized passes that satisfy pass.from_dict(d).to_dict()==d
TWO_WAY_NONPARAM_PASSES = [
    "CommuteThroughMultis",
    "DecomposeArbitrarilyControlledGates",
    "DecomposeBoxes",
    "DecomposeMultiQubitsCX",
    "DecomposeSingleQubitsTK1",
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
    "ZZPhaseToRz",
    "RemoveDiscarded",
    "SimplifyMeasured",
    "RemoveBarriers",
    "DecomposeBridges",
    "CnXPairwiseDecomposition",
    "RemoveImplicitQubitPermutation",
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
    "PauliSimp": standard_pass_dict(
        {
            "name": "PauliSimp",
            "pauli_synth_strat": "Sets",
            "cx_config": "Snake",
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


def check_arc_dict(arc: Architecture, d: dict) -> bool:
    links = [
        {"link": [n1.to_list(), n2.to_list()], "weight": 1} for n1, n2 in arc.coupling
    ]

    if d["links"] != links:
        return False
    else:
        nodes = [Node(n[0], n[1]) for n in d["nodes"]]
        return set(nodes) == set(arc.nodes)


def test_pass_deserialisation_only() -> None:
    # SquashCustom
    def sq(a: float, b: float, c: float) -> Circuit:
        circ = Circuit(1)
        if c != 0:
            circ.Rz(c, 0)
        if b != 0:
            circ.Rx(b, 0)
        if a != 0:
            circ.Rz(a, 0)
        return circ

    squash_pass = SquashCustom({OpType.Rz, OpType.Rx, OpType.Ry}, sq)
    assert squash_pass.to_dict()["StandardPass"]["name"] == "SquashCustom"
    assert set(squash_pass.to_dict()["StandardPass"]["basis_singleqs"]) == {
        "Rz",
        "Rx",
        "Ry",
    }
    # RebaseCustom
    cx = Circuit(2)
    cx.CX(0, 1)
    pz_rebase = RebaseCustom(
        {OpType.CX, OpType.PhasedX, OpType.Rz}, cx, _library._TK1_to_TK1
    )
    assert pz_rebase.to_dict()["StandardPass"]["name"] == "RebaseCustom"
    assert set(pz_rebase.to_dict()["StandardPass"]["basis_allowed"]) == {
        "CX",
        "PhasedX",
        "Rz",
    }
    assert cx.to_dict() == pz_rebase.to_dict()["StandardPass"]["basis_cx_replacement"]

    # RebaseCustomViaTK2
    rebase = RebaseCustom(
        {OpType.XXPhase, OpType.YYPhase, OpType.ZZPhase, OpType.Rx, OpType.Rz},
        lambda a, b, c: Circuit(2).ZZPhase(c).YYPhase(b).XXPhase(a),
        lambda a, b, c: Circuit(1).Rz(c).Rx(b).Rz(a),
    )
    assert rebase.to_dict()["StandardPass"]["name"] == "RebaseCustomViaTK2"
    assert set(rebase.to_dict()["StandardPass"]["basis_allowed"]) == {
        "XXPhase",
        "YYPhase",
        "ZZPhase",
        "Rx",
        "Rz",
    }

    # FullMappingPass
    arc = Architecture([[0, 2], [1, 3], [2, 3], [2, 4]])
    placer = GraphPlacement(arc)
    fm_pass = FullMappingPass(
        arc,
        placer,
        config=[
            LexiLabellingMethod(),
            LexiRouteRoutingMethod(),
            MultiGateReorderRoutingMethod(),
            BoxDecompositionRoutingMethod(),
        ],
    )
    assert fm_pass.to_dict()["pass_class"] == "SequencePass"
    p_pass = fm_pass.get_sequence()[0]
    r_pass = fm_pass.get_sequence()[1]
    np_pass = fm_pass.get_sequence()[2]
    assert np_pass.to_dict()["StandardPass"]["name"] == "NaivePlacementPass"
    assert r_pass.to_dict()["StandardPass"]["name"] == "RoutingPass"
    assert p_pass.to_dict()["StandardPass"]["name"] == "PlacementPass"
    assert check_arc_dict(arc, r_pass.to_dict()["StandardPass"]["architecture"])
    assert p_pass.to_dict()["StandardPass"]["placement"]["type"] == "GraphPlacement"
    assert r_pass.to_dict()["StandardPass"]["routing_config"] == [
        {"name": "LexiLabellingMethod"},
        {
            "name": "LexiRouteRoutingMethod",
            "depth": 10,
        },
        {
            "name": "MultiGateReorderRoutingMethod",
            "depth": 10,
            "size": 10,
        },
        {"name": "BoxDecompositionRoutingMethod"},
    ]
    assert r_pass.to_dict()["StandardPass"]["routing_config"][3] == {
        "name": "BoxDecompositionRoutingMethod"
    }
    # DefaultMappingPass
    dm_pass = DefaultMappingPass(arc)
    assert dm_pass.to_dict()["pass_class"] == "SequencePass"
    p_pass = dm_pass.get_sequence()[0].get_sequence()[0]
    r_pass = dm_pass.get_sequence()[0].get_sequence()[1]
    np_pass = dm_pass.get_sequence()[0].get_sequence()[2]
    d_pass = dm_pass.get_sequence()[1]
    assert d_pass.to_dict()["StandardPass"]["name"] == "DelayMeasures"
    assert d_pass.to_dict()["StandardPass"]["allow_partial"] == False
    assert p_pass.to_dict()["StandardPass"]["name"] == "PlacementPass"
    assert np_pass.to_dict()["StandardPass"]["name"] == "NaivePlacementPass"
    assert r_pass.to_dict()["StandardPass"]["name"] == "RoutingPass"
    assert check_arc_dict(arc, r_pass.to_dict()["StandardPass"]["architecture"])
    assert p_pass.to_dict()["StandardPass"]["placement"]["type"] == "GraphPlacement"
    # DefaultMappingPass with delay_measures=False
    dm_pass = DefaultMappingPass(arc, False)
    assert dm_pass.to_dict()["pass_class"] == "SequencePass"
    assert len(dm_pass.get_sequence()) == 3
    p_pass = dm_pass.get_sequence()[0]
    r_pass = dm_pass.get_sequence()[1]
    np_pass = dm_pass.get_sequence()[2]
    assert p_pass.to_dict()["StandardPass"]["name"] == "PlacementPass"
    assert r_pass.to_dict()["StandardPass"]["name"] == "RoutingPass"
    assert np_pass.to_dict()["StandardPass"]["name"] == "NaivePlacementPass"
    assert check_arc_dict(arc, r_pass.to_dict()["StandardPass"]["architecture"])
    assert p_pass.to_dict()["StandardPass"]["placement"]["type"] == "GraphPlacement"
    # AASRouting
    aas_pass = AASRouting(arc, lookahead=2)
    assert aas_pass.to_dict()["pass_class"] == "SequencePass"
    comppba_plac_pass = aas_pass.get_sequence()[0]
    aasrou_pass = aas_pass.get_sequence()[1]
    assert aasrou_pass.to_dict()["StandardPass"]["name"] == "AASRoutingPass"
    assert check_arc_dict(arc, aasrou_pass.to_dict()["StandardPass"]["architecture"])
    assert (
        comppba_plac_pass.get_sequence()[0].to_dict()["StandardPass"]["name"]
        == "ComposePhasePolyBoxes"
    )
    assert (
        comppba_plac_pass.get_sequence()[1].to_dict()["StandardPass"]["name"]
        == "PlacementPass"
    )
    # CXMappingPass
    cxm_pass = CXMappingPass(arc, placer, directed_cx=True, delay_measures=True)
    assert cxm_pass.to_dict()["pass_class"] == "SequencePass"
    p0 = cxm_pass.get_sequence()[0]
    p1 = cxm_pass.get_sequence()[1]
    assert p0.to_dict()["pass_class"] == "SequencePass"
    assert p1.to_dict()["StandardPass"]["name"] == "DecomposeSwapsToCXs"
    assert p1.to_dict()["StandardPass"]["directed"] == True
    p00 = p0.get_sequence()[0]
    p01 = p0.get_sequence()[1]
    assert p00.to_dict()["pass_class"] == "SequencePass"
    assert p01.to_dict()["StandardPass"]["name"] == "RebaseCustom"
    assert p01.to_dict()["StandardPass"]["basis_cx_replacement"] == cx.to_dict()
    p000 = p00.get_sequence()[0]
    p001 = p00.get_sequence()[1]
    assert p000.to_dict()["pass_class"] == "SequencePass"
    assert p001.to_dict()["StandardPass"]["name"] == "DelayMeasures"
    assert p001.to_dict()["StandardPass"]["allow_partial"] == False
    p0000 = p000.get_sequence()[0]
    p0001 = p000.get_sequence()[1]
    assert p0000.to_dict()["StandardPass"]["name"] == "RebaseCustom"
    assert p0001.to_dict()["pass_class"] == "SequencePass"
    p00010 = p0001.get_sequence()[0]
    p00011 = p0001.get_sequence()[1]
    assert p00010.to_dict()["StandardPass"]["name"] == "PlacementPass"
    assert p00011.to_dict()["StandardPass"]["name"] == "RoutingPass"
    assert check_arc_dict(arc, p00011.to_dict()["StandardPass"]["architecture"])
    # RepeatWithMetricPass
    def number_of_CX(circ: Circuit) -> object:
        return circ.n_gates_of_type(OpType.CX)

    rp = RepeatWithMetricPass(
        SequencePass([CommuteThroughMultis(), RemoveRedundancies()]), number_of_CX
    )
    assert rp.to_dict()["pass_class"] == "RepeatWithMetricPass"
    sps = rp.get_pass().get_sequence()
    assert sps[0].to_dict()["StandardPass"]["name"] == "CommuteThroughMultis"
    assert sps[1].to_dict()["StandardPass"]["name"] == "RemoveRedundancies"
    cx = Circuit(2)
    cx.CX(0, 1)
    cx.CX(1, 0)
    assert number_of_CX(cx) == rp.get_metric()(cx)

    # RepeatUntilSatisfiedPass
    def no_CX(circ: Circuit) -> object:
        return circ.n_gates_of_type(OpType.CX) == 0

    rp = RepeatUntilSatisfiedPass(
        SequencePass([CommuteThroughMultis(), RemoveRedundancies()]), no_CX
    )
    assert rp.to_dict()["pass_class"] == "RepeatUntilSatisfiedPass"
    sps = rp.get_pass().get_sequence()
    assert sps[0].to_dict()["StandardPass"]["name"] == "CommuteThroughMultis"
    assert sps[1].to_dict()["StandardPass"]["name"] == "RemoveRedundancies"
    assert rp.get_predicate().__repr__() == "UserDefinedPredicate"
    assert (
        rp.to_dict()["RepeatUntilSatisfiedPass"]["predicate"]["type"]
        == "UserDefinedPredicate"
    )
