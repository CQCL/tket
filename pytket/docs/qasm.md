# pytket.qasm

{py:class}`~.Circuit` objects can be converted to and from OpenQASM, although we do not support all operations.

However, we do support symbolic parameters of gates, both on import and export.

Any pytket {py:class}`~.Circuit` that is exported to OpenQASM format with `pytket.qasm` should be valid for importing again as a {py:class}`~.Circuit`, making this a convenient file format
to save your {py:class}`~.Circuit` objects.

In addition to the default `qelib1` qasm header, the `hqslib1` header is also supported.
We can set the `header` argument in the qasm conversion functions as follows.

```
from pytket.qasm import circuit_to_qasm_str

qasm_str = circuit_to_qasm_str(circ, header="hqslib1")
```

:::{note}
Unlike pytket backends, the qasm converters do not handle [implicit qubit permutations](https://docs.quantinuum.com/tket/user-guide/manual/manual_circuit.html#implicit-qubit-permutations). In other words if a circuit containing an implicit qubit permutation is converted to a qasm file the implicit permutation will not be accounted for and the circuit will be missing this permutation when reimported.
:::

```{eval-rst}
.. automodule:: pytket.qasm
```

```{eval-rst}
.. automodule:: pytket.qasm.qasm
    :members: circuit_from_qasm, circuit_from_qasm_wasm, circuit_to_qasm, circuit_from_qasm_str, circuit_to_qasm_str, circuit_from_qasm_io, circuit_to_qasm_io, QASMParseError, QASMUnsupportedError
```

```{eval-rst}
.. automodule:: pytket.qasm.grammar
```
