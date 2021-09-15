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

The core functionality of tket is implemented in C++. The source for this
resides in the `bubble` directory.

The source code for pytket (which includes the abstract `Backend` interface used
by the extension modules) resides in the `pytket` directory. This consists of
binder modules (written in C++ and making use of `pybind11` to link to the tket
shared library) and pure Python code for the `Backend` interface and various
utilities.

## How to build tket and pytket

### Prerequisites

#### Build tools

The following compiler toolchains are used to build tket on the CI and are
therefore known to work:

* Linux: gcc-10
* MacOS: apple-clang 12
* Windows: MSVC 19

It is recommended that you use these versions to build locally, as code may
depend on the features they support. The compiler version can be controlled by
setting `CC` and `CXX` in your environment (e.g. `CC=gcc-10` and `CXX=g++-10`),
or on Debian-based Linux systems using `update-alternatives`.

You should also have Python (3.7, 3.8 or 3.9) and `pip` installed. We use
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

If you wish you can set your profile to Debug mode:

```shell
conan profile update settings.build_type=Debug tket
```

#### Enable revisions

In order to pick up the proper revision of the `pybind11` package, it is
currently necessary to do the following (or equivalent):

```shell
conan config set general.revisions_enabled=1
```

#### Test dependencies

A few of the bubble tests require a working LaTeX installation, including
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
### Building tket

#### Method 1

At this point you can run:

```shell
conan create --profile=tket recipes/tket
```

to build the tket library.

To build and run the bubble tests:

```shell
conan create --profile=tket recipes/tket-tests
```

If you want to build them without running them, pass `--test-folder None` to the
`conan` command. (You can still run them manually afterwards.)

There is also a small set of property-based tests which you can build and run
with:

```shell
conan create --profile=tket recipes/tket-proptests
```

Now to build pytket:

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
   conan build recipes/tket --configure --build-folder=build --source-folder=bubble/src
   ```
3. To build:

   ```shell
   conan build recipes/tket --build --build-folder=build
   ```
4. To export to `conan` cache (necessary to build pytket):

   ```shell
   conan export-pkg recipes/tket -f --build-folder=build --source-folder=bubble/src
   ```

## API documentation

The `bubble` (C++) API documentation is available
[here](https://cqcl.github.io/tket/bubble/doc/html/index.html).

The `pytket` (Python) API documentation is available
[here](https://cqcl.github.io/pytket/build/html/index.html).
