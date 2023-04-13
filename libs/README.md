# Dependencies

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

# Building

We recommend building with conan 2, and only conan 2 builds are used on the CI.
However, the libraries may be built using either conan 1 or conan 2.

The build requires `cmake` version 3.23 or above.


## conan 1

```shell
pip install conan~=1.59
```

Create a new profile called `tket`:

```
conan profile new tket --detect
```

(On Linux, a warning will be shown about the `libcxx` configuration. Follow the
recommendation there.)

Add the remote:

```shell
conan remote add tket-libs https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs
```

To build the library `tkxxx`:

```shell
conan create --profile=tket libs/tkxxx --build=missing
```

To build the unit tests, writing package information to `test-tkxxx.json`:

```shell
conan create --profile=tket libs/tkxxx/test --build=missing --json test-tkxxx.json
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

To build the shared library, use the option `-o tkxxx:shared=True` in the
`conan create` command.

## conan 2

```shell
pip install conan~=2.0
```

Use the default profile:

```
conan profile detect
```

Add the remote:

```shell
conan remote add tket-libs https://quantinuumsw.jfrog.io/artifactory/api/conan/tket1-libs
```

To build the library `tkxxx`:

```shell
conan create libs/tkxxx --build=missing -o boost/*:header_only=True
```

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
