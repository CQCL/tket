from typing import Callable, Union

import numpy
from sympy import Symbol

Expression = Union[int, float, numpy.complex128, Symbol]
function = Callable
