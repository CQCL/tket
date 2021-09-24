# Copyright 2019-2021 Cambridge Quantum Computing
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

"""Status classes for circuits submitted to backends.
"""
from typing import Any, Dict, NamedTuple
from enum import Enum


class StatusEnum(Enum):
    """Enumeration for the possible status of a circuit submitted to a backend."""

    COMPLETED = "Circuit has completed. Results are ready."
    QUEUED = "Circuit is queued."
    SUBMITTED = "Circuit has been submitted."
    RUNNING = "Circuit is running."
    CANCELLED = "Circuit has been cancelled."
    ERROR = "Circuit has errored. Check CircuitStatus.message for error message."


class CircuitStatus(NamedTuple):
    """The status of a circuit along with optional long description, \
for example an error message."""

    status: StatusEnum
    message: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Return JSON serializable dictionary representation."""
        return {"status": self.status.name, "message": self.message}

    @classmethod
    def from_dict(cls, dic: Dict[str, Any]) -> "CircuitStatus":
        """Construct from JSON serializable dictionary."""
        invalid = ValueError(f"Dictionary invalid format for CircuitStatus: {dic}")
        if "message" not in dic or "status" not in dic:
            raise invalid
        try:
            status = next(s for s in StatusEnum if dic["status"] == s.name)
        except StopIteration as e:
            raise invalid from e
        return cls(status, dic["message"])


WAITING_STATUS = {StatusEnum.QUEUED, StatusEnum.SUBMITTED, StatusEnum.RUNNING}
