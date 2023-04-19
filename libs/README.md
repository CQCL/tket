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

Building requires conan version 2, and `cmake` version 3.26 or above.

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
