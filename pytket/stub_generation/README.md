# Stub Generation

Python stub files are used to add python type information to the C++ binder modules. The stub files
are generated automatically using the [`pybind11-stubgen`](https://github.com/sizmailov/pybind11-stubgen)
tool and some custom cleanup code within the `regenerate_stubs` script.

To regenerate the stubs after changes to the binder code, (re-)install `pytket` locally
(according to the [README](../README.md)). Then, with `pybind11-stubgen` installed, run the `./regenerate_stubs` script.

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
custom functions above implement the python type equality and make sure the stubs specify the
correct parameter type, a general python `object`.

In python, if `__eq__` is defined explicitly but not `__hash__`, the `__hash__` function is "deleted".
In practice, `__hash__` will simply return a `TypeError`. `pybind11` does this correctly, but for some reason
the stubs generated indicate that `__hash__` has type `None`. `mypy` doesn't like this because it is not compatible
with the type signature of the parent function `object.__hash__`, which is `(self) -> int`. The custom deleted hash
function uses the correct signature and raises the same `TypeError` when called.


### Binding `nlohmann::json` 

> tldr: Don't use `nlohmann::json` as paramter or return types in the binding code. Convert explicitly from and to
> `py::dict`, `py::list`, and `py::object` instead.

The `pybind11_json` package is good for casting between json and py::objects within the
binder code but leads to bad stubs if its custom `type_caster` is used.

### Special classes for binding vector and list to/from typing.Sequence and tuple

The `pybind11` typecasters for `std::vector<T>` and `std::list<T>` will accept any `typing.Sequence`
(i.e. both `list` and `tuple`), but the type name always appears as `list` in the typing information
generated. It is sometimes beneficial to have parameters typed as `typing.Sequence` or `tuple`. Some
examples include:

- mypy will complain about passing a `list[Qubit]` to a method accepting `list[UnitID]`, even though a
  `Qubit` is a `UnitID`. The reason is that the method could append a `UnitID` (which is not necessarily
  a `Qubit`) to the list which could cause problems down the line. If the method is typed as accepting 
  `typing.Sequence`, there is a "guarantee" that the method won't modify the object passed. So mypy is happy.
- In some places we use `std::maps<std::vector<T>, VALUE_TYPE>`, which pybind11 happily
  translates to `dict[list[T], VALUE_TYPE]`, a type that cannot be created in python (since lists are not
  hashable). The typecaster accepts a `tuple[T]`, but mypy will complain if you use a `tuple` here.

To solve these problems there a few types defined in [pytket/binders/include/typecast.hpp](../binders/include/typecast.hpp)
with special typecasters to `typing.Sequence` or `tuple`. One example is `py::tket_custom::SequenceVec<T>`, which
inherits from `std::vector<T>` and automatically casts to `std::vector<T>`. Using it as the parameter of a function
will cause the resulting stub to be typed with a `typing.Sequence` instead of `list`. `py::tket_custom::TupleVec<T>`
also inherits from `std::vector<T>` and can be cast to and from `Tuple[T, ...]`. 
