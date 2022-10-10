from pytket import Circuit, OpType
from pytket.qasm import circuit_from_qasm_str, circuit_to_qasm_str  

c = Circuit(3)
# c.add_gate(OpType.PhasedISWAP, [0.2, 0.3], [0,1 ])
# c.add_gate(OpType.ECR, [0,1])
c.add_gate(OpType.CSWAP, [0,1,2])
print(c.get_commands())
qs = circuit_to_qasm_str(c)
print(qs)

print("\n\n\nTHE RETURN")
cs = circuit_from_qasm_str(qs)
print(cs.get_commands())