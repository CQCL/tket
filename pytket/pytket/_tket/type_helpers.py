from typing import Callable, Union

import numpy
from sympy import Symbol

json = object
Expression = Union[int, float, numpy.complex128, Symbol]
function = Callable
