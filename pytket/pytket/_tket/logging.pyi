import enum


class level(enum.Enum):
    trace = 0
    """all logs"""

    debug = 1
    """debug logs and above"""

    info = 2
    """informational logs and above"""

    warn = 3
    """warnings and above"""

    err = 4
    """error and critical logs only (default)"""

    critical = 5
    """critical logs only"""

    off = 6
    """no logs"""

def set_level(log_level: level) -> None:
    """
    Set the global logging level.

    :param log_level: Desired logging level
    """
