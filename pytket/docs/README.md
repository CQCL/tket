Pytket Docs

# Clone repository

```
git@github.com:CQCL/tket.git
```

```
git submodule update --recursive --init
```

```
cd pytket/docs
```

# Install python deps

```
poetry install
```
# Build html

```
poetry run make html
```
# Serve build

```
npx serve ./build/html
```
