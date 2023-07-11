from pytket import Circuit, OpType, Qubit
from pytket.pauli import Pauli
from pytket.transform import Transform
from pytket.passes import auto_rebase_pass
from typing import List, Tuple, Dict


def update_pauli_measurements(
    clifford_circuit: Circuit,
) -> Tuple[List[Dict[Qubit, Pauli]], List[bool]]:
    """
    Assumes input is n measurements, with a Z measurement on each qubit.
    For each gate in clifford_circuit (dag), updates each measurement.

    Returns a new n element List of QubitPauliString corresponding to the
    n measurements required.
    """
    auto_rebase_pass({OpType.H, OpType.S, OpType.CX, OpType.Rz}).apply(clifford_circuit)
    Transform.RebaseToCliffordSingles().apply(clifford_circuit)

    # Set up a Z measurement for each Qubit
    all_qubit_pauli_string: List[Dict[Qubit, Pauli]] = []
    for q in clifford_circuit.qubits:
        qubit_pauli_string: Dict[Qubit, Pauli] = {
            q: Pauli.I for q in clifford_circuit.qubits
        }
        qubit_pauli_string[q] = Pauli.Z
        all_qubit_pauli_string.append(qubit_pauli_string)

    phases: List[bool] = [False] * clifford_circuit.n_qubits
    # update each measurement for each Qubit
    # identites as on page 31 of:
    # https://www.cs.ox.ac.uk/people/aleks.kissinger/theses/cole-thesis.pdf
    for com in reversed(clifford_circuit.get_commands()):
        for index, measurement in enumerate(all_qubit_pauli_string):
            match com.op.type:
                # updates for H.P.H^
                case OpType.H:
                    q: Qubit = com.qubits[0]
                    match measurement[q]:
                        case Pauli.X:
                            measurement[q] = Pauli.Z
                        case Pauli.Y:
                            phases[index] = not phases[index]
                        case Pauli.Z:
                            measurement[q] = Pauli.X
                # updates for S.P.S^
                case OpType.S:
                    q: Qubit = com.qubits[0]
                    match measurement[q]:
                        case Pauli.X:
                            measurement[q] = Pauli.Y
                        case Pauli.Y:
                            measurement[q] = Pauli.X
                            phases[index] = not phases[index]
                # Updates for Z.P.Z^ = SS.P.S^S^
                case OpType.Z:
                    match measurement[com.qubits[0]]:
                        case Pauli.X | Pauli.Y:
                            phases[index] = not phases[index]
                # Updates for V.P.V^ = HSH.P.H^S^H^
                case OpType.V:
                    q: Qubit = com.qubits[0]
                    match measurement[q]:
                        # N.B For Pauli.X:
                        # H => Pauli.Z
                        # S => Pauli.Z
                        # H => Pauli.X
                        case Pauli.Y:
                            # H => flip phase
                            # S => Pauli.X, flip phase
                            # H => Pauli.Z
                            measurement[q] = Pauli.Z
                        case Pauli.Z:
                            # H => Pauli.X
                            # S => Pauli.Y
                            # H => flip phase
                            measurement[q] = Pauli.Y
                            phases[index] = not phases[index]
                # Updates for X.P.X^
                case OpType.X:
                    match measurement[com.qubits[0]]:
                        case Pauli.Y | Pauli.Z:
                            phases[index] = not phases[index]
                # updates for CX.P.CX^
                case OpType.CX:
                    control: Qubit = com.qubits[0]
                    target: Qubit = com.qubits[1]
                    match measurement[control]:
                        case Pauli.I:
                            match measurement[target]:
                                case Pauli.Z | Pauli.Y:
                                    measurement[control] = Pauli.Z
                        case Pauli.X:
                            match measurement[target]:
                                case Pauli.I:
                                    measurement[target] = Pauli.X
                                case Pauli.X:
                                    measurement[target] = Pauli.I
                                case Pauli.Y:
                                    measurement[control] = Pauli.Y
                                    measurement[target] = Pauli.Z
                                case Pauli.Z:
                                    measurement[control] = Pauli.Y
                                    measurement[target] = Pauli.Y
                                    phases[index] = not phases[index]
                        case Pauli.Y:
                            match measurement[target]:
                                case Pauli.I:
                                    measurement[target] = Pauli.X
                                case Pauli.X:
                                    measurement[target] = Pauli.I
                                case Pauli.Y:
                                    measurement[control] = Pauli.X
                                    measurement[target] = Pauli.Z
                                    phases[index] = not phases[index]
                                case Pauli.Z:
                                    measurement[control] = Pauli.X
                                    measurement[target] = Pauli.Y
                        case Pauli.Z:
                            match measurement[target]:
                                case Pauli.Y | Pauli.Z:
                                    measurement[control] = Pauli.I
                case _:
                    # {H,S,CX} is a universal gateset for Clifford Circuit
                    # Rebase before applying
                    assert False
            # make sure list is updated suitably
            all_qubit_pauli_string[index] = measurement
    return all_qubit_pauli_string, phases
