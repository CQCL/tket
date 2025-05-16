# pytket.backends

Contains {py:class}`~.Backend` abstract class and associated methods. In pytket a {py:class}`~.Backend` represents an interface between pytket and a quantum device or simulator. Different backends are defined in the various pytket extension modules and inherit from the core pytket {py:class}`~.Backend` class.

If you are interested in developing your own {py:class}`~.Backend` or pytket extension then see the [creating backends](https://docs.quantinuum.com/tket/user-guide/examples/backends/creating_backends.html) tutorial.

Documentation relating to Quantinuum Systems device and emulator access can be found at <https://docs.quantinuum.com/systems/>.

See also the [Running on backends](https://docs.quantinuum.com/tket/user-guide/manual/manual_backend.html) section of the pytket user manual.

```{eval-rst}
.. automodule:: pytket.backends
    :members: backend
```

## pytket.backends.backend

```{eval-rst}
.. automodule:: pytket.backends.backend
```

```{eval-rst}
.. autoclass:: Backend
   :special-members: __init__
   :members:
```

```{eval-rst}
.. autoclass:: ResultHandleTypeError
```

## pytket.backends.resulthandle

```{eval-rst}
.. automodule:: pytket.backends.resulthandle
    :members: ResultHandle
```

## pytket.backends.backendresult

```{eval-rst}
.. automodule:: pytket.backends.backendresult
    :members: BackendResult, StoredResult
```

## pytket.backends.status

```{eval-rst}
.. automodule:: pytket.backends.status
    :members: CircuitStatus, StatusEnum
```

## pytket.backends.backendinfo

```{eval-rst}
.. automodule:: pytket.backends.backendinfo
    :members: BackendInfo, fully_connected_backendinfo
```

## pytket.backends.backend_exceptions

```{eval-rst}
.. automodule:: pytket.backends.backend_exceptions
    :members:
```
