pytket.backends
===============

Contains :py:class:`Backend` abstract class and associated methods. In pytket a :py:class:`Backend` represents an interface between pytket and a quantum device or simulator. Different backends are defined in the various pytket extension modules and inherit from the core pytket :py:class:`Backend` class.

Notebook tutorials specific to the :py:class:`QuantinuumBackend` can be found `here <https://github.com/CQCL/pytket-quantinuum/tree/develop/examples>`_.

There are several example tutorials on pytket :py:class:`Backend`\s which can be found `here <https://github.com/CQCL/pytket/tree/main/examples#pytket-examples>`_.

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
