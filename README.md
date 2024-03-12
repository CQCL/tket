# tket

## Introduction

This repository contains the full source code for tket, a quantum SDK.

If you just want to use tket via Python, the easiest way is to install it with
`pip`:

```shell
pip install pytket
```

For getting started using pytket, check out the user manual and examples notebooks.

User manual - https://tket.quantinuum.com/user-manual/
Notebook examples - https://tket.quantinuum.com/examples/

The source content for the manual and example notebooks can be found in the [pytket-docs repository](https://github.com/CQCL/pytket-docs).


In addition to the core pytket package there are pytket extension modules which allow pytket to interface with quantum devices and simulators. Some extensions also provide interoperability with other software libraries such as qiskit, cirq and pennylane.

For a list of available pytket extensions see the [extensions index page](https://tket.quantinuum.com/api-docs/extensions).

These extensions are installed as separate python packages and the source code for each extension lives in its own github repository.

If you would like to build tket yourself and help to improve it, read on!

The codebase is split into two main projects:
 - [tket](tket): the core functionality of tket, optimised for execution speed
   and implemented in C++.
 - [pytket](pytket): the Python interface of tket. This consists of
   binder modules to tket (written in C++ and making use of `pybind11` to link to the tket
   shared library) and pure Python code that defines abstract interfaces 
   used by the extension modules such as the `Backend` and `BackendResult` classes,
   as well as various other utilities.

## How to build tket and pytket

### Prerequisites

#### Build tools

The following compiler toolchains are used to build tket on the CI and are
therefore known to work:

* Linux: gcc-11
* MacOS: apple-clang 14
* Windows: MSVC 19

It is recommended that you use these versions to build locally, as code may
depend on the features they support. The compiler version can be controlled by
setting `CC` and `CXX` in your environment (e.g. `CC=gcc-11` and `CXX=g++-11`),
or on Debian-based Linux systems using `update-alternatives`.

You should also have Python (3.10, 3.11 or 3.12) and `pip` installed. We use
`cmake` and the package manager `conan` to build tket and pytket. The latter can
be installed with `pip`:

```shell
pip install conan
```

You will need at least cmake version 3.26, and conan version 2.


#### Set up `conan` profile

Generate a profile that matches your current machine, and add the required
remote where some dependencies are stored:

```shell
conan profile detect
conan remote add tket-libs https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs --index 0
```

#### Optional: use `ninja` and `ccache`

It is recommended that you also install `ninja` and `ccache` to speed up the
build process. For example with `apt` on Debian/Ubuntu:
```shell
apt install ninja-build ccache
```
Homebrew on MacOS/Linux:
```shell
brew install ninja ccache
```
Chocolatey on Windows:
```shell
choco install ninja ccache
```

On MacOS/Linux:

- If installed, `ccache` is used automatically
- `ninja` must either be set as the default Cmake generator using the following command:
  ```shell
  echo "tools.cmake.cmaketoolchain:generator = Ninja" >> $(conan config home)/global.conf
  ```
  or be specified on a command-by-command basis by providing the argument
  `-c tools.cmake.cmaketoolchain:generator=Ninja` to conan

On Windows:
- Set `ninja` as generator as described above (less reliable than the default `Visual Studio` generator)
- `ccache` will be used automatically *only* when using `Ninja` or `Makefile` as the Cmake generator. It can
  also be used with `Visual Studio` generators by setting the environment
  variable `TKET_VSGEN_CCACHE_EXE` to the path of the `ccache` executable. **Note: this
  must be the path to the actual binary, not a symlink or shim (as used by Chocolatey)**. If using Chocolatey
  to install `ccache`, you can find the path to the binary using `ccache --shimgen-help`


### Building and testing the utility libraries

See the [README](libs/README.md) in the `libs` directory for instructions on
building and testing the utility libraries used by tket (for logging,
random-number generation and so on). This is not necessary if you just want to
build tket or pytket since the recipes or binaries will be automatically
downloaded from the above conan remote.

### Building and testing the tket library

See the [README](tket/README.md) in the `tket` directory for instructions on
building and testing tket as a standalone C++ library.

### Building and testing pytket

See the [README](pytket/README.md) in the `pytket` directory for instructions on
building and testing pytket.

### Nix Support

Tket and pytket are available as a Nix flake, with support for Linux and Apple Silicon systems.
See the [README](nix-support/README.md) in the `nix-support` directory for instructions
on building and testing tket and pytket through Nix, and on how to use it within a Nix project.

## API documentation

The `tket` (C++) API documentation (generated with `doxygen`, and still rather
patchy) is available
[here](https://cqcl.github.io/tket/tket/api/index.html).

The `pytket` (Python) API documentation is available
[here](https://cqcl.github.io/tket/pytket/api/index.html).
