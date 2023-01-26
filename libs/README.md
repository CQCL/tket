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

To build the unit tests inside a `build` directory:

```shell
conan install libs/tkxxx/test --install-folder=build/test-tkxxx --output-folder=build/test-tkxxx --profile=tket
conan build libs/tkxxx/test --configure --build-folder=build/test-tkxxx --source-folder=libs/tkxxx/test
conan build libs/tkxxx/test --build --build-folder=build/test-tkxxx}
```

To run them:

```shell
./build/test-tkxxx/build/Release/test-tkxxx
```

## conan 2

To build the library `tkxxx`:

```shell
conan create --profile=tket libs/tkxxx --build=missing
```

To build the unit tests inside a `build` directory:

```
conan install libs/tkxxx/test --output-folder=build/test-tkxxx --profile=tket
conan build libs/tkxxx/test --output-folder=build/test-tkxxx
```

To run them:

```shell
./build/test-tkxxx/build/Release/test-tkxxx
```
