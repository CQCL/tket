from __future__ import annotations
import typing
__all__ = ['level', 'set_level']
class level:
    """
    Members:
    
      trace : all logs
    
      debug : debug logs and above
    
      info : informational logs and above
    
      warn : warnings and above
    
      err : error and critical logs only (default)
    
      critical : critical logs only
    
      off : no logs
    """
    __members__: typing.ClassVar[dict[str, level]]  # value = {'trace': <level.trace: 0>, 'debug': <level.debug: 1>, 'info': <level.info: 2>, 'warn': <level.warn: 3>, 'err': <level.err: 4>, 'critical': <level.critical: 5>, 'off': <level.off: 6>}
    critical: typing.ClassVar[level]  # value = <level.critical: 5>
    debug: typing.ClassVar[level]  # value = <level.debug: 1>
    err: typing.ClassVar[level]  # value = <level.err: 4>
    info: typing.ClassVar[level]  # value = <level.info: 2>
    off: typing.ClassVar[level]  # value = <level.off: 6>
    trace: typing.ClassVar[level]  # value = <level.trace: 0>
    warn: typing.ClassVar[level]  # value = <level.warn: 3>
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
def set_level(log_level: level) -> None:
    """
    Set the global logging level.
    
    :param log_level: Desired logging level
    """
