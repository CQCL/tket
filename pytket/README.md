# Building and testing pytket

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

Starting from the root directory of this repo, build tket with the required configuration:

```shell
conan create tket --user tket --channel stable --build=missing -o boost/*:header_only=True -o tklog/*:shared=True -o tket/*:shared=True -tf ""
```

Install the other pytket requirements:

```shell
conan create recipes/pybind11
conan create recipes/pybind11_json/all --version 0.2.13
```

Then build the pytket module:

```shell
cd pytket
pip install -e . -v
```

And then to run the Python tests:

```shell
cd tests
pip install -r requirements.txt
pytest
```
