# Nix support for tket and pytket

Tket exposes a Nix flake that you can build, test and use within a Nix environment.
To build and test tket and pytket on your local system, you can run

```
$ nix flake check github:CQCL/tket
```

This will take some time, as it requires the building of a patched symengine,
the compilation of tket and tket tests, and the compilation of pytket's dependencies,
before testing.

You can enter a development environment with tket and pytket available
with use of `nix develop`. For example,

```
$ nix develop github:CQCL/tket

$ python3
Python 3.10.12 (main, Jun  6 2023, 22:43:10) [GCC 12.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from pytket import Circuit
>>> circuit = Circuit(5)
>>> circuit.X(1)
[X q[1]; ]
>>> circuit.X(2)
[X q[1]; X q[2]; ]
```

You can also use tket and pytket within your own flake project.
An example of such a project is given in the [example-flake-project](example-flake-project)
directory.
