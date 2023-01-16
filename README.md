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
pip install cmake conan~=1.53
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

Set the `tket-libs` repository as your remote. (Note that the following commands
affect your conan configuration across all projects, so if you are working on
other projects with conan you will want to revert them afterwards. A simple way
is to back up the file `~/.conan/remotes.json`. You can view your current
remotes list with `conan remote list`.)

```shell
conan remote clean
conan remote add tket-libs https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs
```

Enable revisions:

```shell
conan config set general.revisions_enabled=1
```

We want to build tket and tklog as shared rather than static libraries, so set
this in the profile:

```shell
conan profile update options.tket:shared=True tket
conan profile update options.tklog:shared=True tket
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
conan remove -f "pybind11/*"
conan create --profile=tket recipes/pybind11
```

where the first line serves to remove any version already installed.

### TKET libraries and conan packages

Some TKET functionality has been separated out into self-contained libraries,
as a way to modularize and reduce average build times. These are in
subdirectories of the `libs` directory. We anticipate that their number will
increase as we work towards greater modularization.

If you are using a supported conan configuration (see above under "Build
tools"), you do not need to worry about these, unless you are modifying the code
in them. The main build of TKET will download a pre-built package for each of
them.

If you are using an unsupported configuration, or want to make changes to these
libraries, you will need to build them locally. For example:

```shell
conan create --profile=tket libs/tkrng tket/stable
```

If you make a change to one of these libraries, please increase the version
number and make a PR with that change only: the component will then be tested on
the CI, and on merge to `develop` the new version will be uploaded. Then it will
be possible to update conan requirements to use the new version.

A new version of TKET is uploaded to our conan repo with each push to `develop`
that changes the core library. This process is managed by CI workflows. If you
are making changes only to TKET tests or pytket, you do not need to build TKET
locally: the right version should be downloaded automatically from the conan
repo.

### Building tket

#### Method 1

At this point you can run:

```shell
conan create --profile=tket recipes/tket tket/stable
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
   conan export-pkg recipes/tket tket/${VERSION}@tket/stable -f --build-folder=build --source-folder=tket/src
   ```
   where `${VERSION}` is the tket library version, e.g. `1.0.3`.
 
#### Method 3: Makefile

An alternative build setup and development process is offered through a Makefile in the [`dev-tools`](dev-tools) directory.
This setup can also simplify integration with an IDE supporting `cmake` builds.
See the file [`dev-tools/README.md`](dev-tools/README.md) for instructions and more information.

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
