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

* Linux: gcc-10
* MacOS: apple-clang 13
* Windows: MSVC 19

It is recommended that you use these versions to build locally, as code may
depend on the features they support. The compiler version can be controlled by
setting `CC` and `CXX` in your environment (e.g. `CC=gcc-10` and `CXX=g++-10`),
or on Debian-based Linux systems using `update-alternatives`.

You should also have Python (3.8, 3.9 or 3.10) and `pip` installed. We use
`cmake` and the package manager `conan` to build tket. Both can be installed
with `pip`:

```shell
pip install cmake conan
```

It is recommended that you also install `ninja` and `ccache` to speed up the
build process. For example on Debian/Ubuntu:

```shell
sudo apt install ninja-build ccache
```

#### Set up `conan` profile

Generate a profile that matches your current machine. This profile does not have
to be called `tket`, but if you give it another name you will have to set
`CONAN_TKET_PROFILE` to its name in your environment when you build the Python
module.

```shell
conan profile new tket --detect
```

If this prints a warning about `gcc` ABI compatibility (as it probably will on
Linux), adjust the profile compiler settings with the following command, as
recommended in the warning message:

```shell
conan profile update settings.compiler.libcxx=libstdc++11 tket
```

Add the `tket.conan` repository to your remotes:

```shell
conan remote add tket-conan https://tket.jfrog.io/artifactory/api/conan/tket-conan
```

Enable revisions:

```shell
conan config set general.revisions_enabled=1
```

We want to build shared rather than static libraries, so set this in the
profile:

```shell
conan profile update options.tket:shared=True tket
```

If you wish you can set your profile to Debug mode:

```shell
conan profile update settings.build_type=Debug tket
```

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

### Adding local `pybind11`

There is a known [issue](https://github.com/conan-io/conan-center-index/issues/6605) with using `pybind11`
from the `conan-center` that can lead to a Python crash when importing `pytket`. To remedy this, 
`pybind11` must be installed from the local recipe:

```shell
conan remove -f pybind11/*
conan create --profile=tket recipes/pybind11
```

where the first line serves to remove any version already installed.

### Building symengine

The `symengine` dependency is built from a local conan recipe. Run:

```shell
conan create --profile=tket recipes/symengine
```

to build it. If you are using a conan configuration supported by the CI
(see above under "Build tools"), this is unnecessary as a pre-built package
will be downloaded from the `tket-conan` repository when you build `tket`.

### Building tket

#### Method 1

At this point you can run:

```shell
conan create --profile=tket recipes/tket
```

to build the tket library.

To build and run the tket tests:

```shell
conan create --profile=tket recipes/tket-tests
```

The tests with a running time >=1 second (on a regular modern laptop) are marked as hidden,
tagged with `"[long]"`, and are not run by default. To run the full suite of tests,
add `-o tket-tests:full=True` to the above `conan create` command (or to the tket profile).
The option `-o tket-tests:long=True` can also be used to run only the long tests.

If you want to build the tests without running them, pass `--test-folder None` to the
`conan` command. Then, you can manually run the binary.
If no arguments are provided only the default (short) tests are run.
To run the long tests use the `"[long]"` tag as an argument:

```shell
<package_folder>/bin/test_tket "[long]"
```

To run the full suite manually you need to include also the short tests, like:

```shell
<package_folder>/bin/test_tket "[long],~[long]"
```

A smaller selection of the compiled tests can also be run by passing a filter of the test file name:

```shell
<package_folder>/bin/test_tket -# "[#test_name]"
```

There is also a small set of property-based tests which you can build and run
with:

```shell
conan create --profile=tket recipes/tket-proptests
```

Now to build pytket, first install the `pybind11` headers:

```shell
conan create --profile=tket recipes/pybind11
```

Then build the pytket module:

```shell
cd pytket
pip install -e .
```

And then to run the Python tests:

```shell
cd tests
pytest
```

#### Method 2

In a development cycle, it may save time to break down the `conan create`
command from above into separate build and export commands.

First create a `build` folder in the project root. Then proceed as follows.

1. To install dependencies:

   ```shell
   conan install recipes/tket --install-folder=build --profile=tket --build=missing
   ```
2. To configure the build:

   ```shell
   conan build recipes/tket --configure --build-folder=build --source-folder=tket/src
   ```
3. To build:

   ```shell
   conan build recipes/tket --build --build-folder=build
   ```
4. To export to `conan` cache (necessary to build pytket):

   ```shell
   conan export-pkg recipes/tket -f --build-folder=build --source-folder=tket/src
   ```

## Test coverage

The code coverage of the `tket` tests is reported
[here](https://cqcl.github.io/tket/tket/test-coverage/index.html). This report
is generated weekly from the `develop` branch.

## API documentation

The `tket` (C++) API documentation (generated with `doxygen`, and still rather
patchy) is available
[here](https://cqcl.github.io/tket/tket/api/index.html).

The `pytket` (Python) API documentation is available
[here](https://cqcl.github.io/tket/pytket/api/index.html).
