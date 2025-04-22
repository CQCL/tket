# Stub Generation

Python stub files are used to add python type information to the C++ binder
modules. The stub files are generated automatically using the
[`nanobind.stubgen`](https://github.com/wjakob/nanobind/blob/master/src/stubgen.py)
tool with some customizations within the `regenerate_stubs` script.

To regenerate the stubs after changes to the binder code, (re-)install `pytket`
locally (according to the [README](../README.md)). Then run the
`./regenerate_stubs` script.

To test the generated stubs run:
`mypy --config-file=../mypy.ini --no-incremental  -p pytket -p tests`.
