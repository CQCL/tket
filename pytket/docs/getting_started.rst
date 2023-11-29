Getting Started
===============

The tket compiler is a powerful tool for optimizing and manipulating
platform-agnostic quantum circuits, focused on enabling superior performance on
NISQ (noisy intermediate-scale quantum) devices. The pytket package provides an
API for interacting with tket and transpiling to and from other popular quantum
circuit specifications.

Pytket is compatible with 64-bit Python 3.9, 3.10 and 3.11, on Linux, MacOS
(11.0 or later) and Windows. Install pytket from PyPI using:

::

    pip install pytket

This will install the tket compiler binaries as well as the pytket package. For
those using an older version of pytket, keep up to date by installing with the
``--upgrade`` flag for additional features and bug fixes.

There are separate packages for managing the interoperability between pytket and
other quantum software packages which can also be installed via PyPI. For
details of these, see the
`pytket extensions <https://tket.quantinuum.com/api-docs/extensions>`_ documentation.


The quantum circuit is an abstraction of computation using quantum resources,
designed by initializing a system into a fixed state, then mutating it via
sequences of instructions (gates).

The native circuit interface built into pytket allows us to build circuits and
use them directly.

::

    from pytket import Circuit

    c = Circuit(2, 2) # define a circuit with 2 qubits and 2 bits
    c.H(0)            # add a Hadamard gate to qubit 0
    c.Rz(0.25, 0)     # add an Rz gate of angle 0.25*pi to qubit 0
    c.CX(1,0)         # add a CX gate with control qubit 1 and target qubit 0
    c.measure_all()   # measure qubits 0 and 1, recording the results in bits 0 and 1

Pytket provides many handy shortcuts and higher-level components for building
circuits, including custom gate definitions, circuit composition, gates with
symbolic parameters, and conditional gates.

On the other hand, pytket's flexibile interface allows you to take circuits
defined in a number of languages, including raw source code languages such as
OpenQASM and Quipper, or embedded python frameworks such as Qiskit and Cirq.

::

    from pytket.qasm import circuit_from_qasm

    c = circuit_from_qasm("my_qasm_file.qasm")

Or, if an extension module like ``pytket-qiskit`` is installed:

::

    from qiskit import QuantumCircuit
    from pytket.extensions.qiskit import qiskit_to_tk

    qc = QuantumCircuit()
    # ...
    c = qiskit_to_tk(qc)

See the
`pytket user manual <https://tket.quantinuum.com/user-manual/index.html>`_
for an extensive tutorial on pytket, providing a gentle introduction to its
features and how to run circuits on backend devices, with worked examples.

In pytket there is also a generic :py:class:`Backend` interface. This represents a connection to a quantum device or simulator.
It's possible to run circuits on platforms from different providers through the `extension modules <https://tket.quantinuum.com/api-docs/extensions>`_.

::

    from pytket import Circuit
    from pytket.extensions.qiskit import AerBackend

    # Define Circuit
    circ = Circuit(3)
    circ.H(0).CX(0, 1).CX(0, 2)

    # Initialise Backend and execute the Circuit
    backend = AerBackend()
    result = backend.run_circuit(circ, n_shots=1000)
    print(result.get_counts())


This prints out a summary of readouts (the final values of the classical bits) and their frequencies.

Each pytket :py:class:`Backend` comes with its own default compilation method. This is a recommended sequence of optimisation passes to meet the requirements of the specific :py:class:`Backend`. 

The following code snippet will show how to compile a circuit to run on an IBM device. This requires setting up IBM credentials (see the `credentials guide <https://cqcl.github.io/pytket-qiskit/api/index.html#access-and-credentials>`_).

::

    from pytket.extensions.qiskit import IBMQBackend

    circ = Circuit(3).X(0).CCX(0, 1, 2)
    nairobi_device = IBMQBackend('ibm_nairobi')

    # Compile Circuit to use supported gates of IBMQ Nairobi
    compiled_circ = nairobi_device.get_compiled_circuit(circ)
    result = backend.run_circuit(compiled_circ, n_shots=100)

Here the default compilation pass is applied by :py:meth:`IBMQBackend.get_compiled_circuit`. See `this page <https://cqcl.github.io/pytket-qiskit/api/index.html#default-compilation>`_ for more details.

As an alternative, We can experiment with constructing our own circuit compilation routines in pytket. Passes from the :py:mod:`pytket.passes` module can be applied individually or composed in sequence. 
See the section of the user manual on `circuit compilation <https://tket.quantinuum.com/user-manual/manual_compiler.html>`_ and the corresponding `notebook example <https://github.com/CQCL/pytket/blob/main/examples/compilation_example.ipynb>`_ for more.
