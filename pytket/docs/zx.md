# pytket.zx

```{eval-rst}
.. currentmodule:: pytket._tket.zx
```

```{eval-rst}
.. automodule:: pytket.zx
```

```{eval-rst}
.. automodule:: pytket._tket.zx
```

```{eval-rst}
.. autoclass:: pytket.zx.CliffordGen

   .. automethod:: __init__
   .. autoproperty:: param
```

```{eval-rst}
.. autoclass:: pytket.zx.DirectedGen

   .. automethod:: __init__
   .. autoproperty:: n_ports
   .. autoproperty:: signature
```

```{eval-rst}
.. autoclass:: pytket.zx.Flow

   .. automethod:: __init__
   .. automethod:: c
   .. automethod:: d
   .. automethod:: focus
   .. automethod:: identify_causal_flow
   .. automethod:: identify_focussed_sets
   .. automethod:: identify_pauli_flow
   .. automethod:: odd
   .. autoproperty:: cmap
   .. autoproperty:: dmap
```

```{eval-rst}
.. autoclass:: pytket.zx.PhasedGen

   .. automethod:: __init__
   .. autoproperty:: param
```

```{eval-rst}
.. autoenum:: pytket.zx.QuantumType
```

```{eval-rst}
.. autoclass:: pytket.zx.Rewrite

   .. automethod:: __init__
   .. automethod:: apply
   .. automethod:: basic_wires
   .. automethod:: decompose_boxes
   .. automethod:: extend_at_boundary_paulis
   .. automethod:: extend_for_PX_outputs
   .. automethod:: gadgetise_interior_paulis
   .. automethod:: internalise_gadgets
   .. automethod:: io_extension
   .. automethod:: merge_gadgets
   .. automethod:: parallel_h_removal
   .. automethod:: rebase_to_mbqc
   .. automethod:: rebase_to_zx
   .. automethod:: red_to_green
   .. automethod:: reduce_graphlike_form
   .. automethod:: remove_interior_cliffords
   .. automethod:: remove_interior_paulis
   .. automethod:: repeat
   .. automethod:: self_loop_removal
   .. automethod:: separate_boundaries
   .. automethod:: sequence
   .. automethod:: spider_fusion
   .. automethod:: to_MBQC_diag
   .. automethod:: to_graphlike_form
```

```{eval-rst}
.. autoclass:: pytket.zx.ZXBox

   .. automethod:: __init__
   .. autoproperty:: diagram
   .. autoproperty:: n_ports
   .. autoproperty:: signature
```

```{eval-rst}
.. autoclass:: pytket.zx.ZXDiagram

   .. automethod:: __init__
   .. automethod:: add_vertex
   .. automethod:: add_wire
   .. automethod:: add_zxbox
   .. automethod:: adj_wires
   .. automethod:: check_validity
   .. automethod:: count_vertices
   .. automethod:: count_wires
   .. automethod:: degree
   .. automethod:: free_symbols
   .. automethod:: get_boundary
   .. automethod:: get_name
   .. automethod:: get_qtype
   .. automethod:: get_vertex_ZXGen
   .. automethod:: get_wire_ends
   .. automethod:: get_wire_qtype
   .. automethod:: get_wire_type
   .. automethod:: get_zxtype
   .. automethod:: is_symbolic
   .. automethod:: multiply_scalar
   .. automethod:: neighbours
   .. automethod:: other_end
   .. automethod:: remove_vertex
   .. automethod:: remove_wire
   .. automethod:: set_vertex_ZXGen
   .. automethod:: set_wire_qtype
   .. automethod:: set_wire_type
   .. automethod:: symbol_substitution
   .. automethod:: to_circuit
   .. automethod:: to_doubled_diagram
   .. automethod:: to_graphviz_str
   .. automethod:: wire_at_port
   .. automethod:: wire_between
   .. automethod:: wires_between
   .. autoproperty:: n_vertices
   .. autoproperty:: n_wires
   .. autoproperty:: scalar
   .. autoproperty:: vertices
   .. autoproperty:: wires
```

```{eval-rst}
.. autoclass:: pytket.zx.ZXGen

   .. automethod:: __init__
   .. automethod:: create
   .. autoproperty:: qtype
   .. autoproperty:: type
```

```{eval-rst}
.. autoenum:: pytket.zx.ZXType
```

```{eval-rst}
.. autoclass:: pytket.zx.ZXVert
```

```{eval-rst}
.. autoclass:: pytket.zx.ZXWire
```

```{eval-rst}
.. autoenum:: pytket.zx.ZXWireType
```

```{eval-rst}
.. automethod:: pytket.zx.circuit_to_zx
```

## pytket.zx.tensor_eval

```{eval-rst}
.. automodule:: pytket.zx.tensor_eval
    :members:
```
