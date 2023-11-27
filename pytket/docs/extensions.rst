pytket extensions
=================

The pytket extensions are separate python modules which allow pytket to interface with backends from a range of providers including quantum devices from Quantinuum and IBM.
In pytket a ``Backend`` represents a connection to a QPU (Quantum Processing Unit) or simulator for processing quantum circuits. One can also access additional quantum devices and simulators via the cloud through the extensions for `Azure <https://cqcl.github.io/pytket-qsharp/api/api.html#pytket.extensions.qsharp.AzureBackend>`_ and `Braket <https://tket.quantinuum.com/extensions/pytket-braket/api/api.html#pytket.extensions.braket.BraketBackend>`_ . 

Additionally, the extensions allow pytket to cross-compile circuits from different quantum computing libraries with the extensions for `qiskit <https://tket.quantinuum.com/extensions/pytket-qiskit/api/index.html>`_, `cirq <https://tket.quantinuum.com/extensions/pytket-cirq/api/index.html>`_ and `pennylane <https://tket.quantinuum.com/extensions/pytket-pennylane/api/index.html>`_ . This enables pytket's compilation features to be used in conjunction with other software tools.

The additional modules can be installed adding the extension name to the installation command for pytket. For example pytket-quantinuum can be installed by running

::

   pip install pytket-quantinuum

The types of ``Backend`` available in pytket are the following

Types of Backend
----------------

* **QPUs** - These are real quantum computers that return shots based results. E.g the `QuantinuumBackend <https://tket.quantinuum.com/extensions/pytket-quantinuum/api.html#pytket.extensions.quantinuum.QuantinuumBackend>`_ .
* **Cloud Access** - Cloud backends allow pytket to interface with cloud platforms to access additional QPUs and simulators. E.g `BraketBackend <https://tket.quantinuum.com/extensions/pytket-braket/api.html#pytket.extensions.braket.BraketBackend>`_ .
* **Emulators** - These classically simulate a circuit and produce shots based results. Sometimes emulators use a noise model and have connectivity constraints to emulate real QPUs. E.g. `IBMQEmulatorBackend`_ 
* **Statevector Simulators** - Calculates the pure quantum state prepared by a circuit returning a vector/ndarray. Examples of statevector simulators are the `ForestStateBackend`_ and the `AerStateBackend`_. 
* **Unitary Simulators** - Unitary simulators calculate the unitary operator that is applied by a circuit. A unitary matrix/ndarray is returned `AerUnitaryBackend`_ is an example of such a simulator.
* **Density Matrix Simulators** - These simulators compute the density matrix prepared by a circuit. The result can be a statistical mixture of states in contrast to statevector simulation. E.g. `CirqDensityMatrixSampleBackend`_
* **Other specialised simulators** - There are extensions for simulating specific types of circuit. For instance the `SimplexBackend`_ is designed to simulate Clifford circuits. 

A full list of available pytket backends is shown below.

QPUs
----

`QuantinuumBackend <https://tket.quantinuum.com/extensions/pytket-quantinuum/api.html#pytket.extensions.quantinuum.QuantinuumBackend>`_
- Interface to a remote Quantinuum device or simulator. There are currently two Quantinuum devices offered (H1-1 and H2-1).

`IBMQBackend <https://tket.quantinuum.com/extensions/pytket-qiskit/api.html#pytket.extensions.qiskit.IBMQBackend>`_
- A backend for running circuits on remote IBMQ devices.

`IonQBackend <https://cqcl.github.io/pytket-ionq/api/api.html#pytket.extensions.ionq.IonQBackend>`_
- A backend for running circuits on remote IONQ devices.

`ForestBackend <https://tket.quantinuum.com/extensions/pytket-pyquil/api.html#pytket.extensions.pyquil.ForestBackend>`_
- A backend for running circuits on remote Rigetti devices.

`AQTBackend <https://cqcl.github.io/pytket-aqt/api/api.html#pytket.extensions.aqt.AQTBackend>`_
- Interface to an AQT device or simulator.

`IQMBackend <https://tket.quantinuum.com/extensions/pytket-iqm/api.html#pytket.extensions.iqm.IQMBackend>`_
- Interface to an IQM device or simulator.

Cloud Access
------------

`AzureBackend <https://cqcl.github.io/pytket-qsharp/api/api.html#pytket.extensions.qsharp.AzureBackend>`_
- Backend for running circuits remotely using Azure Quantum devices and simulators.

`BraketBackend <https://tket.quantinuum.com/extensions/pytket-braket/api.html#pytket.extensions.braket.BraketBackend>`_
- Interface to Amazon Braket service.

Emulators
---------

`IBMQEmulatorBackend`_ - A backend which uses the `AerBackend <https://tket.quantinuum.com/extensions/pytket-qiskit/api.html#pytket.extensions.qiskit.AerBackend>`_ to emulate the behavior of IBMQBackend.

`QuantinuumBackend <https://tket.quantinuum.com/extensions/pytket-quantinuum/api.html#pytket.extensions.quantinuum.QuantinuumBackend>`_
- The QuantinuumBackend has two available emulators namely H1-1E and H2-1E. These are device specific emulators for the H1-1 and H2-1 devices. These emulators run remotely on a server.

Statevector Simulators
----------------------

`CirqStateSampleBackend <https://tket.quantinuum.com/extensions/pytket-cirq/api.html#pytket.extensions.cirq.CirqStateSampleBackend>`_
- Backend for Cirq statevector simulator sampling.

`CirqStateSimBackend <https://tket.quantinuum.com/extensions/pytket-cirq/api.html#pytket.extensions.cirq.CirqStateSimBackend>`_
- Backend for Cirq statevector simulator state return.

`AerStateBackend`_ - Backend for running simulations on the Qiskit Aer Statevector simulator.

`ForestStateBackend`_ - State-based interface to a Rigetti device.

`ProjectQBackend <https://tket.quantinuum.com/extensions/pytket-projectq/api.html#pytket.extensions.projectq.ProjectQBackend>`_
- Backend for running statevector simulations on the ProjectQ simulator.

Unitary Simulators
------------------

`AerUnitaryBackend`_ - Backend for running simulations on the Qiskit Aer unitary simulator.

Density Matrix Simulators
-------------------------

`CirqDensityMatrixSampleBackend`_
- Backend for Cirq density matrix simulator sampling.

`CirqDensityMatrixSimBackend <https://tket.quantinuum.com/extensions/pytket-cirq/api.html#pytket.extensions.cirq.CirqDensityMatrixSimBackend>`_
- Backend for Cirq density matrix simulator density_matrix return.

Clifford Simulators
-------------------

`CirqCliffordSampleBackend <https://tket.quantinuum.com/extensions/pytket-cirq/api.html#pytket.extensions.cirq.CirqCliffordSampleBackend>`_
- Backend for Cirq Clifford simulator sampling.

`CirqCliffordSimBackend <https://tket.quantinuum.com/extensions/pytket-cirq/api.html#pytket.extensions.cirq.CirqCliffordSimBackend>`_
- Backend for Cirq Clifford simulator state return.

`SimplexBackend`_- Backend for simulating Clifford circuits using pysimplex.

`StimBackend <https://tket.quantinuum.com/extensions/pytket-stim/api.html#pytket.extensions.stim.StimBackend>`_
- Backend for simulating Clifford circuits using Stim.

Other
-----

`AerBackend <https://tket.quantinuum.com/extensions/pytket-qiskit/api.html#pytket.extensions.qiskit.AerBackend>`_
- Backend for running simulations on the Qiskit Aer QASM simulator. This simulator is noiseless by default but can take a user defined ``NoiseModel``.

`QulacsBackend <https://tket.quantinuum.com/extensions/pytket-qulacs/api.html#pytket.extensions.qulacs.QulacsBackend>`_
- Backend for running simulations of variational quantum circuits on the Qulacs simulator.

`QsharpSimulatorBackend <https://cqcl.github.io/pytket-qsharp/api/api.html#pytket.extensions.qsharp.QsharpSimulatorBackend>`_
- Backend for simulating a circuit using the QDK.

`QsharpToffoliSimulatorBackend <https://cqcl.github.io/pytket-qsharp/api/api.html#pytket.extensions.qsharp.QsharpToffoliSimulatorBackend>`_
- Backend for simulating a Toffoli circuit using the QDK.


.. toctree::
   :caption: Extensions:
   :maxdepth: 0

   pytket-aqt <https://cqcl.github.io/pytket-aqt/api/index.html>
   pytket-braket <https://tket.quantinuum.com/extensions/pytket-braket>
   pytket-cirq <https://tket.quantinuum.com/extensions/pytket-cirq>
   pytket-ionq <https://cqcl.github.io/pytket-ionq/api/index.html>
   pytket-iqm <https://tket.quantinuum.com/extensions/pytket-iqm>
   pytket-pennylane <https://tket.quantinuum.com/extensions/pytket-pennylane>
   pytket-projectq <https://tket.quantinuum.com/extensions/pytket-projectq>
   pytket-pyquil <https://tket.quantinuum.com/extensions/pytket-pyquil>
   pytket-pysimplex <https://tket.quantinuum.com/extensions/pytket-pysimplex>
   pytket-pyzx <https://tket.quantinuum.com/extensions/pytket-pyzx>
   pytket-qir <https://tket.quantinuum.com/extensions/pytket-qir>
   pytket-qiskit <https://tket.quantinuum.com/extensions/pytket-qiskit>
   pytket-qsharp <https://cqcl.github.io/pytket-qsharp/api/index.html>
   pytket-quantinuum <https://tket.quantinuum.com/extensions/pytket-quantinuum>
   pytket-cutensornet <https://tket.quantinuum.com/extensions/pytket-cutensornet> 
   pytket-qulacs <https://tket.quantinuum.com/extensions/pytket-qulacs>
   pytket-qujax <https://cqcl.github.io/pytket-qujax/api/index.html>
   pytket-stim <https://tket.quantinuum.com/extensions/pytket-stim>


.. _pytket: https://tket.quantinuum.com/tket/pytket/api/
.. _Quantinuum: https://quantinuum.com
.. _IBMQEmulatorBackend:  https://tket.quantinuum.com/extensions/pytket-qiskit/api/api.html#pytket.extensions.qiskit.IBMQEmulatorBackend
.. _AerStateBackend: https://tket.quantinuum.com/extensions/pytket-qiskit/api.html#pytket.extensions.qiskit.AerStateBackend
.. _ForestStateBackend: https://tket.quantinuum.com/extensions/pytket-pyquil/api/api.html#pytket.extensions.pyquil.ForestStateBackend
.. _AerUnitaryBackend: https://tket.quantinuum.com/extensions/pytket-qiskit/api/api.html#pytket.extensions.qiskit.AerUnitaryBackend
.. _CirqDensityMatrixSampleBackend: https://tket.quantinuum.com/extensions/pytket-cirq/api/api.html#pytket.extensions.cirq.CirqDensityMatrixSampleBackend
.. _SimplexBackend: https://tket.quantinuum.com/extensions/pytket-simplex/api.html#pytket.extensions.pysimplex.SimplexBackend