from __future__ import annotations
import typing
__all__ = ['complex_to_list', 'list_to_complex']
def complex_to_list(arg0: complex) -> typing.Any:
    """
    Convert complex number to serializable list [real, imag].
    """
def list_to_complex(arg0: typing.Any) -> complex:
    """
    Convert serializable list as output by `complex_to_list` to complex number.
    """
