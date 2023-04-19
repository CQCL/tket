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
`cmake` and the package manager `conan` to build tket. Both can be installed
with `pip`:

```shell
pip install cmake conan
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

(Adding the remote will save time when building for tket for the first time,
but it is possible to build everything locally by first building the recipes
in the [`libs` directory](libs/README).

#### Test dependencies

A few of the tket tests require a working LaTeX installation, including
`latexmk` and the `quantikz` package. By default these are only run on Linux.
Passing `~[latex]` to the test executable will disable them. To install the
Latex dependencies on (Debian flavours of) Linux you can do:

```shell
sudo apt-get install texlive texlive-latex-extra latexmk
mkdir -p ~/texmf/tex/latex
wget http://mirrors.ctan.org/graphics/pgf/contrib/quantikz/tikzlibraryquantikz.code.tex -P ~/texmf/tex/latex
```

The Python tests require a few more packages. These can be installed with:

```shell
pip install -r pytket/tests/requirements.txt
```

### Adding local `pybind11` and `pybind11_json`

There is a known [issue](https://github.com/conan-io/conan-center-index/issues/6605) with using `pybind11`
from the `conan-center` that can lead to a Python crash when importing `pytket`. To remedy this, 
`pybind11` must be installed from the local recipe:

```shell
conan remove -c "pybind11/*"
conan create recipes/pybind11
```

It is also currently necessary to use the local `pybind11_json` recipe, since
the recipe on the `conan-center` is not yet compatible with conan 2:

```shell
conan create recipes/pybind11_json/all --version 0.2.13
```

### TKET libraries and conan packages

Some TKET functionality has been separated out into self-contained libraries,
as a way to modularize and reduce average build times. These are in
subdirectories of the `libs` directory. We anticipate that their number will
increase as we work towards greater modularization.

Recipes and some binaries for these are stored in the `tket-libs` repository;
they can also be built locally. See the [README](libs/README.md) in the `libs`
directory for instructions.

### Building and testing tket

See the [README](tket/README.md) in the `tket` directory for instructions on
building the TKET library and unit tests (but not pytket) using conan.

The tests with a running time >=1 second (on a regular modern laptop) are marked
as hidden, tagged with `"[long]"`, and are not run by default. To run the long
tests use the `"[long]"` tag as an argument:

```shell
<package_folder>/bin/test-tket "[long]"
```

To run the full suite manually you need to include also the short tests, like:

```shell
<package_folder>/bin/test-tket "[long],~[long]"
```

A smaller selection of the compiled tests can also be run by passing a filter of the test file name:

```shell
<package_folder>/bin/test-tket -# "[#test_name]"
```

See the [README](pytket/README.md) in the `pytket` directory for instructions on
how to build pytket.

## Test coverage

### tket

The code coverage of the `tket` tests is reported
[here](https://cqcl.github.io/tket/tket/test-coverage/index.html). This report
is generated weekly from the `develop` branch.

The libraries' coverage (from their own unit tests) is also reported: for
example [tklog](https://cqcl.github.io/tket/tket/tklog-coverage/index.html).
(For other libraries, just replace "tklog" with the library name in the URL.)

In both cases, PRs to `develop` check that the coverage has not decreased, and
merging is blocked until the coverage is at least as good as before.

### pytket

The code coverage of the `pytket` tests is reported
[here](https://cqcl.github.io/tket/pytket/test-coverage/index.html). This report
reflects the coverage of the `develop` branch, and is updated with every push.
The same report can be found in XML format
[here](https://cqcl.github.io/tket/pytket/test-coverage/cov.xml).

Lines and branch coverage results are also checked with every PR to `develop`.

## API documentation

The `tket` (C++) API documentation (generated with `doxygen`, and still rather
patchy) is available
[here](https://cqcl.github.io/tket/tket/api/index.html).

The `pytket` (Python) API documentation is available
[here](https://cqcl.github.io/tket/pytket/api/index.html).
