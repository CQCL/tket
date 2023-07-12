from pytket import Circuit
from pytket.circuit.display import get_circuit_renderer

def entanglement_circuit():
    circuit = Circuit(2)
    circuit.H(0)
    circuit.CX(0, 1)
    circuit.measure_all()
    renderer = get_circuit_renderer()
    renderer.view_browser(circuit)

def teleport_circuit():
    circuit = Circuit(3)
    circuit.H(0)
    circuit.CX(0, 1)
    circuit.CX(1, 2)
    circuit.H(0)
    circuit.measure_all()
    renderer = get_circuit_renderer()
    renderer.view_browser(circuit)

