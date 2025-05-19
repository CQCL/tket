Pytket Docs

This builds the main pages of the pytket documentation site and the API reference. The user manual and examples are maintained in the [pytket-docs](https://github.com/CQCL/pytket-docs/) repository, and docs for each extension is maintained in the corresponding repository.

All documentation pages are written in [MyST Markdown](https://mystmd.org) and built with sphinx using the [MyST parser](https://myst-parser.readthedocs.io/en/latest/index.html), with code cells supported through [myst-nb](https://myst-nb.readthedocs.io/en/latest/).

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
