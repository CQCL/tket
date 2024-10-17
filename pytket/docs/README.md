Pytket Docs

# Clone repository

```
git clone git@github.com:CQCL/tket.git --recurse-submodules
```

# Move to docs directory

```
cd pytket/docs
```

# Install python deps

```
poetry install
```

# Install pytket

The pytket package is not installed as a poetry dependency so needs to be installed seperately

```
poetry run pip install -U pytket
```
You can install a pypi version as above or an editable wheel.

# Build html

```
poetry run bash ./build-docs.sh
```
# Serve built html locally

```
npx serve build
```
