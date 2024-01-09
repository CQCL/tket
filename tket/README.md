# Building and testing tket

## One-off builds

This method is suitable for making local builds when you aren't going to be
editing the code. It is used for CI builds. See below for another method that is
more suitable for local development.

Assuming the prerequisite tools and configuration described in the top-level
README, to build `tket` (from the root directory of the repo):

```shell
conan create tket --user=tket --channel=stable --build=missing -o "boost/*":header_only=True
```

Note: The `create` method of conan does not play well with `ccache` because every
change to tket code results in a unique package revision which is stored in a
unique directory in the cache. To fix this, you can set the ccache `base_dir` to
your home directory by running:

```shell
ccache --set-config base_dir=${HOME}
```

## Building and testing

A few of the tests require a working LaTeX installation, including `latexmk` and
the `quantikz` package. Passing `~[latex]` to the test executable will disable
them. To install the Latex dependencies on (Debian flavours of) Linux you can
do:

```shell
sudo apt-get install texlive texlive-latex-extra latexmk
mkdir -p ~/texmf/tex/latex
wget http://mirrors.ctan.org/graphics/pgf/contrib/quantikz/tikzlibraryquantikz.code.tex -P ~/texmf/tex/latex
```

To build tket and its dependencies and run all unit and property tests, use:

```shell
conan build tket --build=missing -o "boost/*":header_only=True -o with_all_tests=True
```
or `with_test`/`with_proptest` to only build and run the unit tests or proptests, respectively.
 
To build with debug information, run the same `install` and `build` commands with the additional setting `-s build_type=Debug`.

## Test binaries

After building, the test binaries are located in the local
build directories `tket/build/Release` and/or `tket/build/Debug`. 
They can be run from those directories, e.g.:
```shell
cd tket/build/Release
./test/tket-test
./proptest/tket-proptest
```

The tests with a running time >=1 second (on a regular modern laptop) are marked
as hidden, tagged with `"[long]"`, and are not run by default. To run the long
tests use the `"[long]"` tag as an argument:
```shell
cd tket/build/Release
./test/test-tket "[long]"
```

To run the full suite manually you need to include also the short tests, like:
```shell
cd tket/build/Release
./test/test-tket "[long],~[long]"
```

A smaller selection of the compiled tests can also be run by passing a filter of
the test file name:
```shell
cd tket/build/Release
./test/test-tket -# "[#test_name]"
```

## Building without conan

It is possible to build tket without using conan at all: see
[here](../build-without-conan.md) for instructions.

## Generating a test coverage report

If using gcc and with gcovr installed (`pip install gcovr`), a test coverage
report can be generated using the following sequence of commands:

```shell
conan build tket --user=tket --channel=stable -s build_type=Debug --build=missing -o "boost/*":header_only=True -o "tket/*":profile_coverage=True -o "test-tket/*":with_coverage=True -o with_test=True -of build/tket

mkdir test-coverage
gcovr --print-summary --html --html-details -r ./tket --exclude-lines-by-pattern='.*\bTKET_ASSERT\(.*\);' --object-directory=${PWD}/build/tket/build/Debug/CMakeFiles/tket.dir/src -o test-coverage/index.html --decisions
```
