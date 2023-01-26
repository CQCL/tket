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

## conan 1

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

## conan 2

To build the library `tkxxx`:

```shell
conan create --profile=tket libs/tkxxx --build=missing
```

To build the unit tests:

```shell
conan create --profile=tket libs/tkxxx/test --build=missing --format json > test-tkxxx.json
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
