# Dependencies

Some TKET functionality has been separated out into self-contained libraries,
as a way to modularize and reduce average build times. These are in
subdirectories of the `libs` directory. We anticipate that their number will
increase as we work towards greater modularization.

Recipes and some binaries for these are stored in the `tket-libs` repository;
they can also be built locally.

The libraries in this directory have the following dependency graph (where a
downward line means "is required by"):

        tklog
          |
          |
      tkassert   tkrng
          |   \ /  |
          |    X   |
          |   / \  |
    tktokenswap tkwsm

## Building

Assuming the prerequisite tools and configuration described in the top-level
README, to build the library `tkxxx` (from the root directory of the repo):

```shell
conan create libs/tkxxx --build=missing -o boost/*:header_only=True
```

(When using `zsh` shell, consider adding quotes around `boost/*`.)

To build the unit tests:

```shell
conan create libs/tkxxx/test --build=missing -o boost/*:header_only=True --format json > test-tkxxx.json
```

To extract the root package path to environment variable:

```shell
PKGPATH=`./rootpath test-tkxxx.json test-tkxxx`
```

To run the unit tests:

```shell
cd ${PKGPATH}/bin
./test-tkxxx
cd -
```

To build the shared library, use the option `-o tkxxx/*:shared=True` in the
`conan create` command.

## Generating a test coverage report

If using gcc and with gcovr installed (`pip install gcovr`), a test coverage
report can be generated using the following sequence of commands:

```shell
conan install libs/tkxxx -s build_type=Debug --build=missing -o boost/*:header_only=True -o tkxxx/*:profile_coverage=True -of build/tkxxx
conan build libs/tkxxx -s build_type=Debug -o boost/*:header_only=True -o tkxxx/*:profile_coverage=True -of build/tkxxx
conan export-pkg libs/tkxxx -s build_type=Debug -o boost/*:header_only=True -o tkxxx/*:profile_coverage=True -of build/tkxxx -tf ""

conan install libs/tkxxx/test -s build_type=Debug --build=missing -o boost/*:header_only=True -o test-tkxxx/*:with_coverage=True -of build/tkxxx-tests
conan build libs/tkxxx/test -s build_type=Debug --build=missing -o boost/*:header_only=True -o test-tkxxx/*:with_coverage=True -of build/tkxxx-tests
cd ./build/tkxxx-tests/build/Debug
./test-tkxxx
cd -
mkdir tkxxx-coverage
gcovr --print-summary --html --html-details -r ./libs/tkxxx --exclude-lines-by-pattern '.*\bTKET_ASSERT\(.*\);' --object-directory ${PWD}/build/tkxxx/build/Debug/CMakeFiles/tkxxx.dir/src -o tkxxx-coverage/index.html --decisions
```
