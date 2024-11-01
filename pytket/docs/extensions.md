# Extension packages

The pytket extensions are separate python modules which allow pytket to interface with backends from a range of providers including quantum devices from Quantinuum and IBM.
In pytket a `Backend` represents a connection to a QPU (Quantum Processing Unit) or simulator for processing quantum circuits. One can also access additional quantum devices and simulators via the cloud through the extensions for [Braket](inv:pytket-braket:std:doc#index) .

Additionally, the extensions allow pytket to cross-compile circuits from different quantum computing libraries with the extensions for [qiskit](inv:pytket-qiskit:std:doc#index), [cirq](inv:pytket-cirq:std:doc#index) and [pennylane](inv:pytket-pennylane:std:doc#index) . This enables pytket's compilation features to be used in conjunction with other software tools.

The additional modules can be installed adding the extension name to the installation command for pytket. For example pytket-quantinuum can be installed by running

```
pip install pytket-quantinuum
```

The types of `Backend` available in pytket are the following

## Types of Backend

- **QPUs** - These are real quantum computers that return shots based results. E.g the [QuantinuumBackend](inv:#*.QuantinuumBackend).
- **Cloud Access** - Cloud backends allow pytket to interface with cloud platforms to access additional QPUs and simulators. E.g [BraketBackend](inv:#*.BraketBackend).
- **Emulators** - These classically simulate a circuit and produce shots based results. Sometimes emulators use a noise model and have connectivity constraints to emulate real QPUs. E.g. [IBMQEmulatorBackend](inv:#*.IBMQEmulatorBackend)
- **Statevector Simulators** - Calculates the pure quantum state prepared by a circuit returning a vector/ndarray. Examples of statevector simulators are the [ForestStateBackend](inv:#*.ForestStateBackend) and the [AerStateBackend](inv:#*.AerStateBackend).
- **Unitary Simulators** - Unitary simulators calculate the unitary operator that is applied by a circuit. A unitary matrix/ndarray is returned [AerUnitaryBackend](inv:#*.AerUnitaryBackend) is an example of such a simulator.
- **Density Matrix Simulators** - These simulators compute the density matrix prepared by a circuit. The result can be a statistical mixture of states in contrast to statevector simulation. E.g. [CirqDensityMatrixSampleBackend](inv:#*.CirqDensityMatrixSampleBackend)
- **Other specialised simulators** - There are extensions for simulating specific types of circuit. For instance the [SimplexBackend](inv:#*.SimplexBackend) is designed to simulate Clifford circuits.

A full list of available pytket backends is shown below.

## QPUs

[QuantinuumBackend](inv:#*.QuantinuumBackend)
\- Interface to a remote Quantinuum device or simulator. There are currently two Quantinuum devices offered (H1-1 and H2-1).

[IBMQBackend](inv:#*.IBMQBackend)
\- A backend for running circuits on remote IBMQ devices.

[ForestBackend](inv:#*.ForestBackend)
\- A backend for running circuits on remote Rigetti devices.

[AQTBackend](https://cqcl.github.io/pytket-aqt/api/api.html#pytket.extensions.aqt.AQTBackend)
\- Interface to an AQT device or simulator.

[IQMBackend](inv:#*.IQMBackend)
\- Interface to an IQM device or simulator.

## Cloud Access

[BraketBackend](inv:#*.BraketBackend)
\- Interface to Amazon Braket service.

## Emulators

[IBMQEmulatorBackend](inv:#*.IBMQEmulatorBackend) - A backend which uses the [AerBackend](inv:#*.AerBackend) to emulate the behavior of IBMQBackend.

[QuantinuumBackend](inv:#*.QuantinuumBackend)
\- The QuantinuumBackend has two available emulators namely H1-1E and H2-1E. These are device specific emulators for the H1-1 and H2-1 devices. These emulators run remotely on a server.

## Statevector Simulators

[CirqStateSampleBackend](inv:#*.CirqStateSampleBackend)
\- Backend for Cirq statevector simulator sampling.

[CirqStateSimBackend](inv:#*.CirqStateSimBackend)
\- Backend for Cirq statevector simulator state return.

[AerStateBackend](inv:#*.AerStateBackend) - Backend for running simulations on the Qiskit Aer Statevector simulator.

[ForestStateBackend](inv:#*.ForestStateBackend) - State-based interface to a Rigetti device.

[ProjectQBackend](inv:#*.ProjectQBackend)
\- Backend for running statevector simulations on the ProjectQ simulator.

[QuESTBackend](https://docs.quantinuum.com/tket/extensions/pytket-quest/api.html#pytket.extensions.quest.QuESTBackend) Interface to the [QUEST simulator](https://quest.qtechtheory.org/docs/).

## Unitary Simulators

[AerUnitaryBackend](inv:#*.AerUnitaryBackend) - Backend for running simulations on the Qiskit Aer unitary simulator.

## Density Matrix Simulators

[AerDensityMatrixBackend](inv:#*qiskit.AerDensityMatrixBackend) - Backend for density matrix simulation using qiskit Aer. Can take a `NoiseModel` as an optional argument.

[CirqDensityMatrixSampleBackend](inv:#*.CirqDensityMatrixSampleBackend)
\- Backend for Cirq density matrix simulator sampling.

[CirqDensityMatrixSimBackend](inv:#*.CirqDensityMatrixSimBackend)
\- Backend for Cirq density matrix simulator density_matrix return.

[QulacsBackend](inv:#*.QulacsBackend) - This has a configurable density matrix simulation option.

Use `QulacsBackend(result_type="density_matrix")`.

## Clifford Simulators

[CirqCliffordSampleBackend](inv:#*.CirqCliffordSampleBackend)
\- Backend for Cirq Clifford simulator sampling.

[CirqCliffordSimBackend](inv:#*.CirqCliffordSimBackend)
\- Backend for Cirq Clifford simulator state return.

[SimplexBackend](inv:#*.SimplexBackend)- Backend for simulating Clifford circuits using pysimplex.

[StimBackend](inv:#*.StimBackend)
\- Backend for simulating Clifford circuits using Stim.

## Other

[AerBackend](inv:#*.AerBackend)
\- Backend for running simulations on the Qiskit Aer QASM simulator. This simulator is noiseless by default but can take a user defined `NoiseModel`.

[QulacsBackend](inv:#*.QulacsBackend)
\- Backend for running simulations of variational quantum circuits on the Qulacs simulator.

```{toctree}
:caption: 'Extensions:'
:maxdepth: 0

pytket-aqt <https://cqcl.github.io/pytket-aqt/api/index.html>
pytket-azure <https://docs.quantinuum.com/tket/extensions/pytket-azure>
pytket-braket <https://docs.quantinuum.com/tket/extensions/pytket-braket>
pytket-cirq <https://docs.quantinuum.com/tket/extensions/pytket-cirq>
pytket-iqm <https://docs.quantinuum.com/tket/extensions/pytket-iqm>
pytket-pennylane <https://docs.quantinuum.com/tket/extensions/pytket-pennylane>
pytket-projectq <https://docs.quantinuum.com/tket/extensions/pytket-projectq>
pytket-pyquil <https://docs.quantinuum.com/tket/extensions/pytket-pyquil>
pytket-pysimplex <https://docs.quantinuum.com/tket/extensions/pytket-pysimplex>
pytket-pyzx <https://docs.quantinuum.com/tket/extensions/pytket-pyzx>
pytket-qir <https://docs.quantinuum.com/tket/extensions/pytket-qir>
pytket-qiskit <https://docs.quantinuum.com/tket/extensions/pytket-qiskit>
pytket-quantinuum <https://docs.quantinuum.com/tket/extensions/pytket-quantinuum>
pytket-quest <https://docs.quantinuum.com/tket/extensions/pytket-quest>
pytket-cutensornet <https://docs.quantinuum.com/tket/extensions/pytket-cutensornet>
pytket-qulacs <https://docs.quantinuum.com/tket/extensions/pytket-qulacs>
pytket-qujax <https://docs.quantinuum.com/tket/extensions/pytket-qujax>
pytket-stim <https://docs.quantinuum.com/tket/extensions/pytket-stim>
```

