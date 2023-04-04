# Building

We recommend building with conan 2, and only conan 2 builds are used on the CI.
However, tket may be built using either conan 1 or conan 2.

The build requires `cmake` version 3.23 or above.

## conan 1

To build `tket`:

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

It is currently necessary to make a small modification to `tket/CMakeLists.txt`:
change `find_package(Eigen3 CONFIG REQUIRED)` to
`find_package(eigen CONFIG REQUIRED)`; and change
`target_link_libraries(tket PRIVATE Eigen3::Eigen)` to
`target_link_libraries(tket PRIVATE eigen::eigen)`.

Then, to build `tket`:

```shell
conan create --profile=tket tket/ tket/stable --build=missing
```

To build the unit tests, writing package information to `test-tket.json`:

```shell
conan create --profile=tket tket/test --build=missing --json test-tket.json
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

To build the property tests, writing package information to
`proptest-tket.json`:

```shell
conan create --profile=tket tket/proptest --build=missing --json proptest-tket.json
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

To build `tket`:

```shell
conan create tket --user=tket --channel=stable --build=missing -o boost/*:header_only=True
```

To build the unit tests:

```shell
conan create tket/test --build=missing  -o boost/*:header_only=True --format json > test-tket.json
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
