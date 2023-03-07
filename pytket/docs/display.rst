pytket.circuit.display
==================================

Contains several functions for rendering interactive circuit diagrams.

.. note:: Rendering circuits with ``pytket.circuit.display`` requires an internet connection. Using the pytket circuit renderer offline can be done by installing the `pytket-offline_display extension <https://github.com/CQCL/pytket-offline-renderer>`_.

    ::

        pip install pytket-offline_display 

.. automodule:: pytket.circuit.display
    :members: render_circuit_jupyter, render_circuit_as_html, view_browser

.. jupyter-execute::

    from pytket import Circuit
    from pytket.circuit.display import render_circuit_jupyter

    circ = Circuit(2) # Define Circuit
    circ.H(0).H(1).CX(0, 1).Rz(0.4, 1).CX(0, 1).H(0).H(1)
    render_circuit_jupyter(circ) # Render interactive display