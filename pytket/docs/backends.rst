pytket.backends
===============

Contains :py:class:`Backend` abstract class and associated methods. In pytket a :py:class:`Backend` represents an interface between pytket and a quantum device or simulator. Different backends are defined in the various pytket extension modules and inherit from the core pytket :py:class:`Backend` class.

Notebook tutorials on pytket backends

1. `Backends example <https://github.com/CQCL/pytket/blob/main/examples/backends_example.ipynb>`_ - Covers the different backends from the ``pytket-qiskit`` extension as well as compiling circuits, noise models and calculating expectation values.  
2. `Comparing simulators <https://github.com/CQCL/pytket/blob/main/examples/comparing_simulators.ipynb>`_ - Discusses capabilities of simulator backends from ``pytket-qiskit``, ``pytket-pyquil``, ``pytket-qsharp``, ``pytket-projectq`` and ``pytket-qulacs``.
3. `Qiskit integration <https://github.com/CQCL/pytket/blob/main/examples/qiskit_integration.ipynb>`_ - Covers compiling to IBM backends as well as using the `TketBackend <https://cqcl.github.io/pytket-qiskit/api/api.html#pytket.extensions.qiskit.tket_backend.TketBackend>`_ to run qiskit circuits.
4. `Creating Backends <https://github.com/CQCL/pytket/blob/main/examples/creating_backends.ipynb>`_ - Tutorial showing how to define your own pytket :py:class:`Backend`. This could be used as a guide for developing new pytket extensions. 

See also the `Running on backends <https://cqcl.github.io/pytket/manual/manual_backend.html>`_ section of the pytket user manual.

.. automodule:: pytket.backends
    :members: backend

pytket.backends.backend
~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: Backend
   :special-members: __init__
   :members:

pytket.backends.resulthandle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.backends.resulthandle
    :members: ResultHandle

pytket.backends.backendresult
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.backends.backendresult
    :members: BackendResult, StoredResult

pytket.backends.status
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.backends.status
    :members: CircuitStatus, StatusEnum

pytket.backends.backendinfo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: pytket.backends.backendinfo
    :members: BackendInfo
