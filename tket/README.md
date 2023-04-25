# Building and testing tket

Assuming the prerequisite tools and configuration described in the top-level
README, to build `tket` (from the root directory of the repo):

```shell
conan create tket --user=tket --channel=stable --build=missing -o boost/*:header_only=True
```

(When using `zsh` shell, consider adding quotes around `boost/*`.)

To build the unit tests:

```shell
conan create tket/test --build=missing  -o boost/*:header_only=True --format json > test-tket.json
```

A few of the tests require a working LaTeX installation, including `latexmk` and
the `quantikz` package. Passing `~[latex]` to the test executable will disable
them. To install the Latex dependencies on (Debian flavours of) Linux you can
do:

```shell
sudo apt-get install texlive texlive-latex-extra latexmk
mkdir -p ~/texmf/tex/latex
wget http://mirrors.ctan.org/graphics/pgf/contrib/quantikz/tikzlibraryquantikz.code.tex -P ~/texmf/tex/latex
```

To extract the root package path to environment variable:

```shell
PKGPATH=`./rootpath test-tket.json test-tket`
```

To run the unit tests:

```shell
cd ${PKGPATH}/bin
./test-tket
cd -
```

The tests with a running time >=1 second (on a regular modern laptop) are marked
as hidden, tagged with `"[long]"`, and are not run by default. To run the long
tests use the `"[long]"` tag as an argument:

```shell
./test-tket "[long]"
```

To run the full suite manually you need to include also the short tests, like:

```shell
./test-tket "[long],~[long]"
```

A smaller selection of the compiled tests can also be run by passing a filter of
the test file name:

```shell
./test-tket -# "[#test_name]"
```


To build the property tests:

```shell
conan create tket/proptest --build=missing  -o boost/*:header_only=True --format json > proptest-tket.json
```

To extract the root package path to environment variable:

```shell
PKGPATH=`./rootpath proptest-tket.json proptest-tket`
```

To run the property tests:

```shell
cd ${PKGPATH}/bin
./proptest-tket
cd -
```

## Generating a test coverage report

If using gcc and with gcovr installed (`pip install gcovr`), a test coverage
report can be generated using the following sequence of commands:

```shell
conan install tket --user=tket --channel=stable -s build_type=Debug --build=missing -o boost/*:header_only=True -o tket/*:profile_coverage=True -of build/tket
conan build tket --user=tket --channel=stable -s build_type=Debug -o boost/*:header_only=True -o tket/*:profile_coverage=True -of build/tket
conan export-pkg tket --user=tket --channel=stable -s build_type=Debug -o boost/*:header_only=True -o tket/*:profile_coverage=True -of build/tket -tf ""

conan install tket/test -s build_type=Debug --build=missing -o boost/*:header_only=True -o test-tket/*:with_coverage=True -of build/tket-tests
conan build tket/test -s build_type=Debug --build=missing -o boost/*:header_only=True -o test-tket/*:with_coverage=True -of build/tket-tests
cp tket/test/src/test_circuits/*.json ./build/tket-tests/build/Debug
cp tket/test/src/test_architectures/*.json ./build/tket-tests/build/Debug

cd ./build/tket-tests/build/Debug
./test-tket
cd -
mkdir test-coverage
gcovr --print-summary --html --html-details -r ./tket --exclude-lines-by-pattern '.*\bTKET_ASSERT\(.*\);' --object-directory ${PWD}/build/tket/build/Debug/CMakeFiles/tket.dir/src -o test-coverage/index.html --decisions
```
