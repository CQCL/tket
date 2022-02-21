

from pytket import Circuit
from pytket.predicates import CompilationUnit

circ_dict = {'bits': [['c', [0]], ['c', [1]], ['c', [2]], ['c', [3]], ['c', [4]], ['c', [5]], ['c', [6]], ['c', [7]], ['c', [8]], ['c', [9]], ['c', [10]], ['c', [11]], ['tk_SCRATCH_BIT', [0]], ['tk_SCRATCH_BIT', [1]], ['tk_SCRATCH_BIT', [2]], ['tk_SCRATCH_BIT', [3]], ['tk_SCRATCH_BIT', [4]], ['tk_SCRATCH_BIT', [5]], ['tk_SCRATCH_BIT', [6]], ['tk_SCRATCH_BIT', [7]], ['tk_SCRATCH_BIT', [8]], ['tk_SCRATCH_BIT', [9]], ['tk_SCRATCH_BIT', [10]], ['tk_SCRATCH_BIT', [11]], ['tk_SCRATCH_BIT', [12]], ['tk_SCRATCH_BIT', [13]], ['tk_SCRATCH_BIT', [14]], ['tk_SCRATCH_BIT', [15]], ['tk_SCRATCH_BIT', [16]], ['tk_SCRATCH_BIT', [17]], ['tk_SCRATCH_BIT', [18]], ['tk_SCRATCH_BIT', [19]], ['tk_SCRATCH_BIT', [20]], ['tk_SCRATCH_BIT', [21]], ['tk_SCRATCH_BIT', [22]], ['tk_SCRATCH_BIT', [23]]], 'commands': [{'args': [['q', [9]], ['q', [5]]], 'op': {'type': 'CZ'}}, {'args': [['q', [8]], ['q', [6]]], 'op': {'type': 'CZ'}}, {'args': [['q', [7]], ['q', [10]]], 'op': {'type': 'CZ'}}, {'args': [['q', [6]], ['q', [0]]], 'op': {'type': 'CZ'}}, {'args': [['q', [2]], ['q', [4]]], 'op': {'type': 'CZ'}},  {'args': [['q', [9]], ['q', [8]]], 'op': {'type': 'CZ'}}, {'args': [['q', [2]], ['q', [1]]], 'op': {'type': 'CZ'}}, {'args': [['q', [5]], ['q', [3]]], 'op': {'type': 'CZ'}}, {'args': [['q', [6]], ['q', [4]]], 'op': {'type': 'CZ'}}, {'args': [['q', [9]], ['q', [11]]], 'op': {'type': 'CZ'}}, {'args': [['q', [1]], ['q', [0]]], 'op': {'type': 'CZ'}}, {'args': [['q', [3]], ['q', [2]]], 'op': {'type': 'CZ'}},{'args': [['q', [1]], ['c', [1]]], 'op': {'type': 'Measure'}}, {'args': [['c', [1]], ['c', [4]], ['tk_SCRATCH_BIT', [6]]], 'op': {'box': {'exp': {'args': [{'args': [['c', [1]], False], 'op': 'BitWiseOp.XOR'}, ['c', [4]]], 'op': 'BitWiseOp.XOR'}, 'id': '00be1bef-0d05-485f-b53f-470111bb26e8', 'n_i': 2, 'n_io': 0, 'n_o': 1, 'type': 'ClassicalExpBox'}, 'type': 'ClassicalExpBox'}},{'args': [['tk_SCRATCH_BIT', [6]], ['q', [3]]], 'op': {'conditional': {'op': {'type': 'X'}, 'value': 1, 'width': 1}, 'type': 'Conditional'}}], 'implicit_permutation': [[['q', [0]], ['q', [0]]], [['q', [1]], ['q', [1]]], [['q', [2]], ['q', [2]]], [['q', [3]], ['q', [3]]], [['q', [4]], ['q', [4]]], [['q', [5]], ['q', [5]]], [['q', [6]], ['q', [6]]], [['q', [7]], ['q', [7]]], [['q', [8]], ['q', [8]]], [['q', [9]], ['q', [9]]], [['q', [10]], ['q', [10]]], [['q', [11]], ['q', [11]]]], 'phase': '0.0', 'qubits': [['q', [0]], ['q', [1]], ['q', [2]], ['q', [3]], ['q', [4]], ['q', [5]], ['q', [6]], ['q', [7]],['q', [8]], ['q', [9]], ['q', [10]], ['q', [11]]]}

circ = Circuit.from_dict(circ_dict)



from pytket.passes import FullMappingPass, RoutingPass, DefaultMappingPass
from pytket.architecture import SquareGrid
from pytket.placement import place_with_map

arc = SquareGrid(4,4)
nodes = arc.nodes
qubits = circ.qubits
# print(circ)

qmap = {qubits[0]: nodes[5], 
        qubits[1]: nodes[8], 
        qubits[2]: nodes[0], 
        qubits[3]: nodes[9], 
        qubits[4]: nodes[4],
        qubits[5]: nodes[10],
        qubits[6]: nodes[1],
        qubits[7]: nodes[2],
        qubits[8]: nodes[3],
        qubits[9]: nodes[7],
        qubits[10]: nodes[6],
        qubits[11]: nodes[11]}

place_with_map(circ, qmap)
cu = CompilationUnit(circ)

RoutingPass(SquareGrid(4,4)).apply(cu)

# DefaultMappingPass(SquareGrid(4,4)).apply(cu)

# FullMappingPass(SquareGrid(4,4)).apply(cu)
# print(cu.circuit)
# for x in cu.initial_map:
#   print(x, cu.initial_map[x])

# print(cu.final_map)
print("Finish Routing")
print(cu.circuit.get_commands())
# print(cu.circuit.qubits)








# q[0] gridNode[1, 1, 0]
# q[1] gridNode[2, 0, 0]
# q[2] gridNode[0, 0, 0]
# q[3] gridNode[2, 1, 0]
# q[4] gridNode[1, 0, 0]
# q[5] gridNode[2, 2, 0]
# q[6] gridNode[0, 1, 0]
# q[7] gridNode[0, 2, 0]
# q[8] gridNode[0, 3, 0]
# q[9] gridNode[1, 3, 0]
# q[10] gridNode[1, 2, 0]
# q[11] gridNode[2, 3, 0]

# 0 gridNode[0, 0, 0]
# 1 gridNode[0, 1, 0]
# 2 gridNode[0, 2, 0]
# 3 gridNode[0, 3, 0]
# 4 gridNode[1, 0, 0]
# 5 gridNode[1, 1, 0]
# 6 gridNode[1, 2, 0]
# 7 gridNode[1, 3, 0]
# 8 gridNode[2, 0, 0]
# 9 gridNode[2, 1, 0]
# 10 gridNode[2, 2, 0]
# 11 gridNode[2, 3, 0]
# 12 gridNode[3, 0, 0]
# 13 gridNode[3, 1, 0]
# 14 gridNode[3, 2, 0]
# 15 gridNode[3, 3, 0]
