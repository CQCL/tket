# Building and testing pytket

Assuming the prerequisite tools and configuration described in the top-level
README, to build `tket` with the required configuration (from the root directory
of the repo):

```shell
conan create tket --user tket --channel stable --build=missing -o boost/*:header_only=True -o tklog/*:shared=True -o tket/*:shared=True -tf ""
```

There is a known
[issue](https://github.com/conan-io/conan-center-index/issues/6605) with using
`pybind11` from the `conan-center` that can lead to a Python crash when
importing `pytket`. To remedy this, `pybind11` must be installed from the local
recipe:

```shell
conan remove -c "pybind11/*"
conan create recipes/pybind11
```

It is also currently necessary to use the local `pybind11_json` recipe, since
the recipe on the `conan-center` is not yet compatible with conan 2:

```shell
conan create recipes/pybind11_json/all --version 0.2.13
```

Then build the pytket module:

```shell
cd pytket
pip install -e . -v
```

The Python tests require a few more packages. These can be installed with:

```shell
pip install -r tests/requirements.txt
```

And then to run the Python tests:

```shell
cd tests
pytest
```

To generate a test coverage report:

```shell
pytest --hypothesis-seed=1 --cov=../pytket --cov-branch --cov-report=html --cov-report=xml:htmlcov/cov.xml
```

## Building without conan

It is possible to build pytket without using conan at all: see
[here](../build-without-conan.md) for instructions.
