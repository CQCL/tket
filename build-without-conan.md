# Building tket and pytket without conan

It is possible to build tket (including tests) and pytket without using conan to
manage the packages. This may be useful in some contexts, but is not used for
official releases. The scripts below were developed for Linux and have not been
tested on other platforms.

For convenience, set a few environment variables:

* `TKET_DIR`: the root directory of this repo
* `INSTALL_DIR`: a directory where libraries will be installed
* `TMP_DIR`: a directory where external sources will be downloaded and built

## Installing external dependencies

The versions should match the current requirements as specified in the relevant
`conanfile.py`. Those specified below are correct for pytket 1.15.0.

(Some of these may be installed already in system directories.)

### `boost`

```
cd ${TMP_DIR}
wget -O boost_1_83_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.83.0/boost_1_83_0.tar.gz/download
tar xzvf boost_1_83_0.tar.gz
cd boost_1_83_0/
./bootstrap.sh --prefix=${INSTALL_DIR}
./b2
./b2 install
```

### `gmp`

```
cd ${TMP_DIR}
wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.bz2
bzip2 -dk gmp-6.3.0.tar.bz2
tar xvf gmp-6.3.0.tar
cd gmp-6.3.0/
./configure --prefix=${INSTALL_DIR} --enable-cxx=yes
make
make check
make install
```

### `symengine`

```
cd ${TMP_DIR}
wget https://github.com/symengine/symengine/releases/download/v0.11.2/symengine-0.11.2.tar.gz
tar xzvf symengine-0.11.2.tar.gz
cd symengine-0.11.2/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF ..
cmake --build .
cmake --install .
````

### `eigen`

```
cd ${TMP_DIR}
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
bzip2 -dk eigen-3.4.0.tar.bz2
tar xvf eigen-3.4.0.tar
cd eigen-3.4.0/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} ..
cmake --build .
cmake --install .
```

### `nlohmann_json`

```
cd ${TMP_DIR}
wget https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
tar xvf json.tar.xz
cd json/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DJSON_BuildTests=OFF ..
cmake --build .
cmake --install .
```

### `catch2` (needed for tests)

```
cd ${TMP_DIR}
wget https://github.com/catchorg/Catch2/archive/refs/tags/v3.5.0.tar.gz
tar xzvf v3.5.0.tar.gz
cd Catch2-3.5.0/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} ..
cmake --build .
cmake --install .
```

### `rapidcheck` (needed for property tests)

```
cd ${TMP_DIR}
wget https://github.com/emil-e/rapidcheck/archive/8fafda42e732164db58003e542196e94a28481f9.zip
unzip 8fafda42e732164db58003e542196e94a28481f9.zip
cd rapidcheck-8fafda42e732164db58003e542196e94a28481f9/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} ..
cmake --build .
cmake --install .
```

### `pybind11` (needed for pytket)

```
cd ${TMP_DIR}
wget https://github.com/pybind/pybind11/archive/refs/tags/v2.11.1.tar.gz
tar xzvf v2.11.1.tar.gz
cd pybind11-2.11.1/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DPYBIND11_TEST=OFF ..
cmake --build .
cmake --install .
```

### `pybind11_json` (needed for pytket)

```
cd ${TMP_DIR}
wget https://github.com/pybind/pybind11_json/archive/refs/tags/0.2.13.tar.gz
tar xzvf 0.2.13.tar.gz
cd pybind11_json-0.2.13/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} ..
cmake --build .
cmake --install .
```

## TKET utility libraries

### `tklog`

```
cd ${TKET_DIR}/libs/tklog/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} ..
cmake --build .
cmake --install .
```

### `tklog` (shared library, needed for pytket)

```
cd ${TKET_DIR}/libs/tklog/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DBUILD_SHARED_LIBS=1 ..
cmake --build .
cmake --install .
```

### `tkrng`

```
cd ${TKET_DIR}/libs/tkrng/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_POSITION_INDEPENDENT_CODE=ON ..
cmake --build .
cmake --install .
```

### `tkassert`

```
cd ${TKET_DIR}/libs/tkassert/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_POSITION_INDEPENDENT_CODE=ON ..
cmake --build .
cmake --install .
```

### `tkwsm`

```
cd ${TKET_DIR}/libs/tkwsm/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_POSITION_INDEPENDENT_CODE=ON ..
cmake --build .
cmake --install .
```

### `tktokenswap`

```
cd ${TKET_DIR}/libs/tktokenswap/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_POSITION_INDEPENDENT_CODE=ON ..
cmake --build .
cmake --install .
```

## `tket`

```
cd ${TKET_DIR}/tket
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} ..
cmake --build .
cmake --install .
```

## `tket` with tests and proptests

This needs the static `tklog` and `tket` libraries to be installed.

```
cd ${TKET_DIR}/tket/test
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DBUILD_TKET_TEST=ON -DBUILD_TKET_PROPTEST=ON ..
cmake --build .
cmake --install .
```

## `tket` (shared library, for pytket)


```
cd ${TKET_DIR}/tket
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DBUILD_SHARED_LIBS=1 -DCMAKE_BUILD_WITH_INSTALL_RPATH=1 -DCMAKE_INSTALL_RPATH="\${ORIGIN}" ..
cmake --build .
cmake --install .
```

## `pytket`

This needs the shared `tklog` and `tket` libraries to be installed.

```
cd ${TKET_DIR}/pytket
NO_CONAN=1 pip install -v -e .
```
