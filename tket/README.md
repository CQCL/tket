# Building and testing tket

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

To build `tket`:

```shell
conan create tket --user=tket --channel=stable --build=missing -o boost/*:header_only=True
```
(When using `zsh` shell, consider adding quotes around `boost/*`.)

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
