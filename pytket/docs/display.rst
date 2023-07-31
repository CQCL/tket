pytket.circuit.display
==================================

Contains several functions for rendering interactive circuit diagrams.

.. note:: Rendering circuits with ``pytket.circuit.display`` requires an internet connection. Using the pytket circuit renderer offline can be done by installing the `pytket-offline-display extension <https://github.com/CQCL/pytket-offline-renderer>`_.

    ::

        pip install pytket-offline-display

.. automodule:: pytket.circuit.display
    :members: get_circuit_renderer, CircuitRenderer

Example usage:
--------------

.. jupyter-execute::

    from pytket import Circuit
    from pytket.circuit.display import get_circuit_renderer

    circuit_renderer = get_circuit_renderer() # Instantiate a circuit renderer
    circuit_renderer.set_render_options(zx_style=True) # Configure render options
    circuit_renderer.condense_c_bits = False # You can also set the properties on the instance directly
    print("Render options:")
    print(circuit_renderer.get_render_options()) # View currently set render options

    circuit_renderer.min_height = "300px" # Change the display height

    circ = Circuit(2,2) # Define Circuit
    circ.H(0).H(1).CX(0, 1).Rz(0.4, 1).CX(0, 1).H(0).H(1).measure_all()

    circuit_renderer.render_circuit_jupyter(circ) # Render interactive display

If you are happy with the default render options, you can import the render
functions directly:

.. jupyter-execute::

        from pytket import Circuit
        from pytket.circuit.display import render_circuit_jupyter

        circ = Circuit(2,2) # Define Circuit
        circ.H(0).H(1).CX(0, 1).Rz(0.4, 1).CX(0, 1).H(0).H(1).measure_all()

        render_circuit_jupyter(circ) # Render with default options

This same diagram can be rendered with the offline renderer as follows
 
    ::

        from pytket.extensions.offline_display import get_circuit_renderer, render_circuit_jupyter

        custom_renderer = get_circuit_renderer()
        custom_renderer.render_circuit_jupyter(circ) # Render configurable display as above
        render_circuit_jupyter(circ) # Render using default options

