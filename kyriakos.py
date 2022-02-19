from pytket import Circuit
from pytket.predicates import CompilationUnit

circ_dict = {
    "bits": [
        ["c", [0]],
        ["c", [1]],
        ["c", [2]],
        ["c", [3]],
        ["c", [4]],
        ["c", [5]],
        ["tk_SCRATCH_BIT", [0]],
        ["tk_SCRATCH_BIT", [1]],
        ["tk_SCRATCH_BIT", [2]],
    ],
    "commands": [
        {"args": [["q", [1]], ["q", [3]]], "op": {"type": "CZ"}},
        {"args": [["q", [1]], ["q", [2]]], "op": {"type": "CZ"}},
        {"args": [["q", [4]], ["q", [3]]], "op": {"type": "CZ"}},
        {"args": [["q", [1]], ["q", [0]]], "op": {"type": "CZ"}},
        {
            "args": [
                ["q", [0]],
                ["q", [1]],
                ["q", [2]],
                ["q", [3]],
                ["q", [4]],
                ["q", [5]],
                ["c", [0]],
                ["c", [1]],
                ["c", [2]],
                ["c", [3]],
                ["c", [4]],
                ["c", [5]],
            ],
            "op": {
                "signature": [
                    "Q",
                    "Q",
                    "Q",
                    "Q",
                    "Q",
                    "Q",
                    "C",
                    "C",
                    "C",
                    "C",
                    "C",
                    "C",
                ],
                "type": "Barrier",
            },
        },
        {"args": [["q", [0]]], "op": {"type": "H"}},
        {"args": [["q", [1]]], "op": {"type": "H"}},
        {"args": [["q", [0]], ["c", [0]]], "op": {"type": "Measure"}},
        {"args": [["q", [1]], ["c", [1]]], "op": {"type": "Measure"}},
        {
            "args": [
                ["q", [0]],
                ["q", [1]],
                ["q", [2]],
                ["q", [3]],
                ["q", [4]],
                ["c", [0]],
                ["c", [1]],
                ["c", [2]],
                ["c", [3]],
                ["c", [4]],
            ],
            "op": {
                "signature": ["Q", "Q", "Q", "Q", "Q", "C", "C", "C", "C", "C"],
                "type": "Barrier",
            },
        },
        {"args": [["q", [0]]], "op": {"type": "Reset"}},
        {"args": [["q", [1]]], "op": {"type": "Reset"}},
        {"args": [["q", [4]]], "op": {"params": ["-0.25"], "type": "Rz"}},
        {
            "args": [["c", [1]], ["tk_SCRATCH_BIT", [0]]],
            "op": {
                "box": {
                    "exp": {"args": [["c", [1]], False], "op": "BitWiseOp.XOR"},
                    "id": "12c10add-5033-437b-b911-f939f97203ed",
                    "n_i": 1,
                    "n_io": 0,
                    "n_o": 1,
                    "type": "ClassicalExpBox",
                },
                "type": "ClassicalExpBox",
            },
        },
        {
            "args": [["c", [0]], ["tk_SCRATCH_BIT", [1]]],
            "op": {
                "box": {
                    "exp": {"args": [["c", [0]], False], "op": "BitWiseOp.XOR"},
                    "id": "7d9e1fc7-dac1-4c52-8202-c480ef1897e0",
                    "n_i": 1,
                    "n_io": 0,
                    "n_o": 1,
                    "type": "ClassicalExpBox",
                },
                "type": "ClassicalExpBox",
            },
        },
        {
            "args": [["c", [0]], ["tk_SCRATCH_BIT", [2]]],
            "op": {
                "box": {
                    "exp": {"args": [["c", [0]], False], "op": "BitWiseOp.XOR"},
                    "id": "319a085c-42b6-4aa7-8348-cd588f6aa3f5",
                    "n_i": 1,
                    "n_io": 0,
                    "n_o": 1,
                    "type": "ClassicalExpBox",
                },
                "type": "ClassicalExpBox",
            },
        },
        {
            "args": [["tk_SCRATCH_BIT", [0]], ["q", [2]]],
            "op": {
                "conditional": {"op": {"type": "X"}, "value": 1, "width": 1},
                "type": "Conditional",
            },
        },
        {
            "args": [["tk_SCRATCH_BIT", [2]], ["q", [3]]],
            "op": {
                "conditional": {"op": {"type": "Z"}, "value": 1, "width": 1},
                "type": "Conditional",
            },
        },
        {
            "args": [["tk_SCRATCH_BIT", [1]], ["q", [2]]],
            "op": {
                "conditional": {"op": {"type": "Z"}, "value": 1, "width": 1},
                "type": "Conditional",
            },
        },
        {"args": [["q", [3]]], "op": {"params": ["-0.25"], "type": "Rz"}},
        {"args": [["q", [2]]], "op": {"params": ["-0.25"], "type": "Rz"}},
    ],
    "implicit_permutation": [
        [["q", [0]], ["q", [0]]],
        [["q", [1]], ["q", [1]]],
        [["q", [2]], ["q", [2]]],
        [["q", [3]], ["q", [3]]],
        [["q", [4]], ["q", [4]]],
        [["q", [5]], ["q", [5]]],
    ],
    "phase": "0.0",
    "qubits": [["q", [0]], ["q", [1]], ["q", [2]], ["q", [3]], ["q", [4]], ["q", [5]]],
}


circ = Circuit.from_dict(circ_dict)

cu = CompilationUnit(circ)
print(cu)

from pytket.passes import FullMappingPass, RoutingPass, DefaultMappingPass
from pytket.architecture import SquareGrid

DefaultMappingPass(SquareGrid(4, 4)).apply(cu)
print(cu.circuit)
