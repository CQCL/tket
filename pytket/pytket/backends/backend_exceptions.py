# Copyright 2019-2024 Cambridge Quantum Computing
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import Union, Optional
from .resulthandle import ResultHandle


class CircuitNotValidError(Exception):
    """Raised when a submitted circuit does not satisfy all predicates"""

    def __init__(self, message: Union[str, int], failed_pred: Optional[str] = None):
        if isinstance(message, int):
            message = (
                "Circuit with index {0} in submitted does not satisfy "
                "{1} (try compiling with backend.get_compiled_circuits first)."
            ).format(message, failed_pred or "all predicates")
        super().__init__(message)


class CircuitNotRunError(Exception):
    """Raised when a result is retrieved corresponding to a handle that has not been
    executed"""

    def __init__(self, handle: ResultHandle):
        super().__init__(
            "Circuit corresponding to {0!r} ".format(handle)
            + "has not been run by this backend instance."
        )


class InvalidResultType(Exception):
    """Raised when a BackendResult instance cannot produce the required result type."""

    def __init__(self, result_type: str):
        super().__init__(
            "BackendResult cannot produce result of type {}.".format(result_type)
        )
