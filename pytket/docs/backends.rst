pytket.backends
===============

Contains :py:class:`Backend` abstract class and associated methods. In pytket a :py:class:`~pytket.backends.Backend` represents an interface between pytket and a quantum device or simulator. Different backends are defined in the various pytket extension modules and inherit from the core pytket :py:class:`~pytket.backends.Backend` class.

If you are interested in developing your own :py:class:`Backend` or pytket extension then see the `creating backends <https://docs.quantinuum.com/tket/user-guide/examples/backends/creating_backends.html>`_ tutorial.

Documentation relating to Quantinuum systems device and emulator access can be found at https://docs.quantinuum.com/systems/.

See also the `Running on backends <https://docs.quantinuum.com/tket/user-guide/manual/manual_backend.html>`_ section of the pytket user manual.

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
