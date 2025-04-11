# Building and testing pytket

Assuming the prerequisite tools and configuration described in the top-level
README, to build `tket` with the required configuration (from the root directory
of the repo):

```shell
conan create tket --user=tket --channel=stable --build=missing -o "boost/*":header_only=True -o "tklog/*":shared=True -o "tket/*":shared=True -tf ""
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

## Stub generation

Python type stubs are generated from the `nanobind` modules using mypy's `stubgen`. Changes to the
binding code under [binders](binders) (or, in some cases, to `tket` itself) may require stub regeneration.
See [stub_generation/README.md](stub_generation/README.md) for more information.

## Building without conan

It is possible to build pytket without using conan at all: see
[here](../build-without-conan.md) for instructions.
