pytket
======

``pytket`` is a python module for interfacing with tket, a quantum computing toolkit and optimising compiler developed by `Quantinuum`_. We currently support circuits and device architectures from
`numerous providers <https://cqcl.github.io/pytket-extensions/api/index.html>`_, allowing the
tket tools to be used in conjunction with projects on their platforms.

``pytket`` is available for Python 3.9, 3.10 and 3.11, on Linux, MacOS and
Windows. To install, run

::

    pip install pytket

.. note::
    On M1-based Macs running in native (arm64) mode, this command may fail
    because of an issue installing ``scipy``. To fix this:

    1. Install `brew <https://brew.sh/>`_ (if you haven't already);
    2. ``brew install openblas``;
    3. ``pip install -U pip wheel``;
    4. ``OPENBLAS="$(brew --prefix openblas)" pip install scipy``;
    5. ``pip install pytket``.

If you have issues installing ``pytket`` please visit the `installation troubleshooting <https://cqcl.github.io/tket/pytket/api/install.html>`_ page.

To use ``pytket``, you can simply import the appropriate modules into your python code or in an interactive Python notebook. We can build circuits directly using the ``pytket`` interface by creating a blank circuit and adding gates in the order we want to apply them.

::

    from pytket import Circuit

    circ = Circuit(2,2) # define a circuit with 2 qubits and 2 bits
    circ.H(0)           # add a Hadamard gate to qubit 0
    circ.Rz(0.25, 0)    # add an Rz gate of angle 0.25*pi to qubit 0
    circ.CX(1,0)        # add a CX gate with control qubit 1 and target qubit 0
    circ.measure_all()  # measure qubits 0 and 1, recording the results in bits 0 and 1

Some of the extension modules define :py:class:`Backend` s, allowing the circuits to be run on simulators or real quantum hardware. For example, ``pytket-qiskit`` grants access to the :py:class:`AerBackend` simulator which can sample from measurements.

::

    from pytket.extensions.qiskit import AerBackend

    backend = AerBackend()                                 # connect to the backend
    compiled_circ = backend.get_compiled_circuit(circ)     # compile the circuit to satisfy the backend's requirements
    handle = backend.process_circuit(compiled_circ, 100)   # submit the job to run the circuit 100 times
    counts = backend.get_result(handle).get_counts()       # retrieve and summarise the results
    print(counts)

This prints out a summary of readouts (the final values of the classical bits) and their frequencies.

::

    {(0, 0): 49, (1, 0): 51}

See the `Getting Started`_ page for a basic tutorial on using
``pytket``. To get more in depth on features, see the `examples`_. See the `pytket user manual <https://cqcl.github.io/pytket/manual/index.html>`_ for an extensive introduction to ``pytket`` functionality and how to use it.

Extensions
~~~~~~~~~~

To use pytket in conjunction with other platforms you must download an
additional separate module for each. Each one of these adds either some new
methods to the ``pytket`` package to convert between the circuit
representations, or some new backends to submit circuits to within ``pytket``.

Extensions are separate python packages can be installed using ``pip``. The installation command is ``pip install pytket-X`` where ``X`` is the name of the extension.

To install the ``pytket-quantinuum`` package use the following command.
::
    
    pip install pytket-quantinuum

The extensions supported by tket are described
`here <https://cqcl.github.io/pytket-extensions/api/index.html>`_.

How to cite
~~~~~~~~~~~

If you wish to cite tket in any academic publications, we generally recommend citing our `software overview paper <https://doi.org/10.1088/2058-9565/ab8e92>`_ for most cases.

If your work is on the topic of specific compilation tasks, it may be more appropriate to cite one of our other papers:

- `"On the qubit routing problem" <https://doi.org/10.4230/LIPIcs.TQC.2019.5>`_ for qubit placement (aka allocation, mapping) and routing (aka swap network insertion, connectivity solving).
- `"Phase Gadget Synthesis for Shallow Circuits" <https://doi.org/10.4204/EPTCS.318.13>`_ for representing exponentiated Pauli operators in the ZX calculus and their circuit decompositions.
- `"A Generic Compilation Strategy for the Unitary Coupled Cluster Ansatz" <https://arxiv.org/abs/2007.10515>`_ for sequencing of terms in Trotterisation and Pauli diagonalisation.

We are also keen for others to benchmark their compilation techniques against us. We recommend checking our `benchmark repository <https://github.com/CQCL/tket_benchmarking>`_ for examples on how to run basic benchmarks with the latest version of ``pytket``. Please list the release version of ``pytket`` with any benchmarks you give, and feel free to get in touch for any assistance needed in setting up fair and representative tests.

User Support
~~~~~~~~~~~~

If you have problems with the use of tket or you think that you have found a bug there are several ways to contact us:

- You can write an issue on `github <https://github.com/CQCL/tket/issues>`_ with details of the problem and we will pick that up. Github issues are the preferred way to report bugs with tket or request features. You can also have a look on that page so see if your problem has already been reported by someone else.
- We have a slack channel for community discussion and support. You can join by following `this link <https://tketusers.slack.com/join/shared_invite/zt-18qmsamj9-UqQFVdkRzxnXCcKtcarLRA#/shared-invite/email>`_
- Write an email to tket-support@quantinuum.com and ask for help with your problem.
- There is also a tag on `quantum computing stack exchange <https://quantumcomputing.stackexchange.com/questions/tagged/pytket>`_ for questions relating to pytket.

We are really thankful for all help to fix bugs in tket. Usually you will get an answer from someone in the development team of tket soon.

LICENCE
~~~~~~~

Licensed under the `Apache 2 License <http://www.apache.org/licenses/LICENSE-2.0>`_.

.. _Getting Started: getting_started.html
.. _examples: https://github.com/CQCL/pytket/tree/main/examples
.. _Quantinuum: https://www.quantinuum.com/

.. toctree::
    :caption: Introduction:
    :maxdepth: 1

    getting_started.rst
    changelog.rst
    install.rst
    opensource.rst

.. toctree::
    :caption: More Documentation:
    :maxdepth: 1
    
    Manual <https://cqcl.github.io/pytket/manual/index.html>
    Extensions <https://cqcl.github.io/pytket-extensions/api/index.html>
    Example notebooks <https://github.com/CQCL/pytket/tree/main/examples>

.. toctree::
    :caption: API Reference:
    :maxdepth: 2

    backends.rst
    circuit.rst
    pauli.rst
    passes.rst
    predicates.rst
    partition.rst
    qasm.rst
    quipper.rst
    architecture.rst
    placement.rst
    mapping.rst
    tableau.rst
    transform.rst
    tailoring.rst
    wasm.rst
    zx.rst
    utils.rst
    logging.rst
    config.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
