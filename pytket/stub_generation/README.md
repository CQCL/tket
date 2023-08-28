# Stub Generation

Python stub files are used to add python type information to the C++ binder modules. The stub files
are generated automatically using mypy's `stubgen` tool and some custom cleanup code
within the `regenerate_stubs` script.

To regenerate the stubs after changes to the binder code, (re-)install `pytket` locally
(according to the [README](../README.md)). Then, with `mypy` installed, run the `./regenerate_stubs` script.

To test the generated stubs run: `mypy --config-file=../mypy.ini --no-incremental  -p pytket -p tests`.

There are a number of caveats and workarounds needed to make the stubs generated from the `pybind11`
modules acceptable for `mypy` (see below).


## Binder code tips and tricks

### When defining `__eq__`, `__ne__`

> tldr: Use the custom equality checkers `py_equals` and `py_not_equals` from [py_operators.hpp](../binders/include/py_operators.hpp) to define `__eq__` and `__ne__` for
> a C++ class in the binder code. If you define `__eq__` for a class and don't define a corresponding `__hash__`
> function, define `__hash__` using the custom "deleted" hash function from [deleted_hash.hpp](../binders/include/deleted_hash.hpp), i.e. add
> ```c++
> .def("__hash__", deletedHash<MyClass>, deletedHashStringj)
> ```
> to your `pybind11` bound class. See examples in binding code.

In C++, an equality check between different types is typically a compile time error, while, in python, it
typically results in a runtime `False`. Thus, in python, the equality function of a class can be called
for any python object. `pybind11` offers the syntax `pybind11::self ==  pybind11::self` to automatically
generate this functionality, but the stubs generated are still too restrictive and `mypy` will complain. The
custom functions above implement the python type equality and make sure the stubs specify that the
correct parameter type, a general python `object`.

In python, if `__eq__` is defined explicitly but not `__hash__`, the `__hash__` function is "deleted".
In practice, `__hash__` will simply return a `TypeError`. `pybind11` does this correctly, but for some reason
the stubs generated indicate that `__hash__` has type `None`. `mypy` doesn't like this because it is not compatible
with the type signature of the parent function `object.__hash__`, which is `(self) -> int`. The custom deleted hash
function uses the correct signature and raises the same `TypeError` when called.


### Binding `nlohmann::json` 

> tldr: Don't use `nlohmann::json` as paramter or return types in the binding code. Convert explicitly from and to
> `py::dict`, `py::list`, and `py::object` instead.

If `pybind11_json` is good for casting between json and py::objects within the
binder code but leads to bad stubs if the custom `type_caster` is used.

### Docstrings of methods with no parameters 

Say you have a method that should have the signature `my_method(self) -> ReturnType`.
For some very weird reason, if you include the method name with empty parenthesis (`my_method()`) within
the docstring, e.g. within an example, extra stub overloads will be generated with signature `(self) -> Any`. `mypy` will
complain because it doesn't make sense. Probably best to avoid doing this. If you must, delete the extra stubs
within the stub generation script, as done currently for the `Circuit().depth` method.
