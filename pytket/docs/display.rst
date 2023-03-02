pytket.circuit.display
==================================

Conatins several functions for rendering interactive circuit diagrams.

.. note:: Rendering circuits in jupyter notebooks with ``pytket.circuit.display`` requires an internet connection as the source code is in a `remote github repository <https://github.com/CQCL/pytket-circuit-renderer>`_ and is not installed with the pytket python package. Using the circuit renderer offline is can be done by installing the `pytket-offline_display extension <https://github.com/CQCL/pytket-offline-renderer>`_. 

.. jupyter-execute::

    from pytket import Circuit
    from pytket.circuit.display import render_circuit_jupyter

    circ = Circuit(2)
    circ.H(0).H(1).CX(0, 1).Rz(0.4, 1).CX(0, 1).H(0).H(1)
    render_circuit_jupyter(circ)


.. automodule:: pytket.circuit.display
    :members: render_circuit_jupyter, render_circuit_as_html, view_browser
