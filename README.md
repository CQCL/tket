# tket

## Introduction

This repository contains the full source code for tket, a quantum SDK.

If you just want to use tket via Python, the easiest way is to install it with
`pip`:

```shell
pip install pytket
```

For full API documentation, as well as a comprehensive user manual and a
selection of example notebooks, please follow the links from the
[pytket](https://github.com/CQCL/pytket) main page.

Note that the various pytket extensions (which allow pytket to interface with
other software packages and with quantum devices) live in the separate
[pytket-extensions](https://github.com/CQCL/pytket-extensions) repository.

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

You should also have Python (3.9, 3.10 or 3.11) and `pip` installed. We use
`cmake` and the package manager `conan` to build tket and pytket. The latter can
be installed with `pip`:

```shell
pip install conan
```

You will need at least cmake version 3.26, and conan version 2.

It is recommended that you also install `ninja` and `ccache` to speed up the
build process. For example on Debian/Ubuntu:

```shell
sudo apt install ninja-build ccache
```

#### Set up `conan` profile

Generate a profile that matches your current machine, and add the required
remote where some dependencies are stored:

```shell
conan profile detect
conan remote add tket-libs https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs
```

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

## API documentation

The `tket` (C++) API documentation (generated with `doxygen`, and still rather
patchy) is available
[here](https://cqcl.github.io/tket/tket/api/index.html).

The `pytket` (Python) API documentation is available
[here](https://cqcl.github.io/tket/pytket/api/index.html).
