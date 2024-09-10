pytket.backends
===============

Contains :py:class:`Backend` abstract class and associated methods. In pytket a :py:class:`Backend` represents an interface between pytket and a quantum device or simulator. Different backends are defined in the various pytket extension modules and inherit from the core pytket :py:class:`Backend` class.

There are several `example notebooks <https://tket.quantinuum.com/examples>`_ on pytket :py:class:`Backend`\s. If you are interested in developing your own :py:class:`Backend` or pytket extension then see the `creating backends <https://tket.quantinuum.com/user-guide/examples/creating_backends.html>`_ tutorial.

Notebook tutorials specific to the :py:class:`QuantinuumBackend` can be found `here <https://github.com/CQCL/pytket-quantinuum/tree/develop/examples>`_.

See also the `Running on backends <https://tket.quantinuum.com/user-guide/manual/manual_backend.html>`_ section of the pytket user manual.

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
