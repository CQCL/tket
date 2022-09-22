# Contributing to TKET

Pull requests are welcome. To make a PR, first fork the repo, make your proposed
changes on the `develop` branch, and open a PR from your fork. If it passes
tests and is accepted after review, it will be merged in.

When adding a new feature, please add tests for it. When fixing a bug, please
add a test that demonstrates the fix.

## Code style

### C++

C++20 features may be used whenever they are supported by all the compilers
used on the CI (listed in the README).

All C++ code should be formatted with `clang-format` (v13) using the
configuration file `.clang-format` in the root directory. This is checked on
the CI. The script `do-clang-format` will run this over all C++ files in the
repository and fix them up.

You can mark sections of code with `// clang-format off` and
`// clang-format on` to tell the tool to skip them, e.g. when it is helpful to
have data laid out in a certain way.

In other matters of style, please try to follow the
[Google style guide](https://google.github.io/styleguide/cppguide.html).

Declarations in header files should have Doxygen-style documentation,
sufficient to make it clear what each object is and what each method does.

The macro `TKET_ASSERT(...)` is used to check logical correctness of the code
(_not_ to catch situations that can occur through user error). The condition is
always checked, even in Release builds, so please avoid using it with
conditions that incur significant computational overhead. Please place
`TKET_ASSERT(...);` on a line containing no other code; then the coverage
checker will ignore it when calculating line and branch coverage statistics.

### Python

#### Formatting

All Python code should be formatted using
[black](https://black.readthedocs.io/en/stable/), with default options. This is
checked on the CI.

#### Type annotation

On the CI, [mypy](https://mypy.readthedocs.io/en/stable/) is used as a static
type checker and all submissions must pass its checks. You should therefore run
`mypy` locally on any changed files before submitting a PR. The following
command will perform all required checks:

```shell
cd pytket
mypy --config-file=mypy.ini -p pytket -p tests
```

#### Linting

We use [pylint](https://pypi.org/project/pylint/) on the CI to check compliance
with a set of style requirements (listed in `pytket/.pylintrc`). You should run
`pylint` over any changed files before submitting a PR, to catch any issues.
