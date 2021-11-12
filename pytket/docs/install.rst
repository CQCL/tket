Installation Troubleshooting
==================================

General Installation Instructions
------------------------------------------
To install the core ``pytket`` package, run

``pip install pytket``

If the ``pip`` command is linked to a Python 2 distribution, you will need to replace it with ``pip3`` and use ``python3`` for running your code, or use ``python -m pip``.

The additional modules ``pytket-X`` (where ``X`` is ``qiskit``, ``pyquil``, etc.) can be installed in the same way, as

``pip install pytket-X``

You can update your installation to the most recent version using

``pip install --upgrade pytket``


Frequently Asked Questions
--------------------------

Is there a build of ``pytket`` for my system?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The non-commercial version of ``pytket``, and most additional modules available through PyPI, are built to work on Linux, MacOS and Windows with Python versions 3.7, 3.8, or 3.9, and ``pip`` version 20.0.0+.


Do all versions of ``pytket`` work with Windows?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``pytket`` versions 0.6 and above are compatible with most recent versions of Windows. LTSC versions of Windows are not supported and an `issue <https://github.com/CQCL/pytket/issues/36>`_ has been reported with these.


``pytket`` installed but modules mentioned in the docs could not be found. Why?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Some of the modules concerning interaction with another software package (e.g. ``qiskit``, ``pyquil``, etc.) need to be installed separately.

Otherwise, this may be that the version of ``pytket`` you obtained is not the most recent version released. We only publish our docs for the most recent version of ``pytket``.


When I ran ``pip install pytket``, I could only get an old version. What gives?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
At a couple of points in the development of the software, we had to increase the system requirements. Obtaining an old version from PyPI is likely the result of that being the most recent version compliant with your system.

If you received version 0.3.0 or 0.4.2, it is likely that you are using an old version of ``pip`` that cannot accept the more recent Linux builds. Try running ``pip install --upgrade pip`` to upgrade it to the most recent version and upgrade ``pytket``.


I've tried the recommended actions here and it still won't work! What can I do?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Our `examples repository <http://github.com/CQCL/pytket>`_ on GitHub has an issue tracker for current issues. You might find others who have had similar problems there. If not, feel free to add an issue describing your problem and our dev team will try to diagnose it and get back to you as soon as possible.
