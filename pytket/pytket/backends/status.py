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

"""Status classes for circuits submitted to backends.
"""
from collections.abc import Callable
from datetime import datetime
from enum import Enum
from typing import Any, NamedTuple


class StatusEnum(Enum):
    """Enumeration for the possible status of a circuit submitted to a backend."""

    COMPLETED = "Circuit has completed. Results are ready."
    QUEUED = "Circuit is queued."
    SUBMITTED = "Circuit has been submitted."
    RUNNING = "Circuit is running."
    RETRYING = "Circuit is being retried."
    CANCELLING = "Cancellation has been requested."
    CANCELLED = "Circuit has been cancelled."
    ERROR = "Circuit has errored. Check CircuitStatus.message for error message."


class CircuitStatus(NamedTuple):
    """The status of a circuit along with an optional description.

    Optionally can also include extra fields such as:
    * Detailed error information.
    * Timestamps for changes in status.
    * Queue position.
    """

    status: StatusEnum
    message: str = ""
    error_detail: str | None = None

    # Timestamp for when a status was last entered.
    completed_time: datetime | None = None
    queued_time: datetime | None = None
    submitted_time: datetime | None = None
    running_time: datetime | None = None
    cancelled_time: datetime | None = None
    error_time: datetime | None = None

    queue_position: int | None = None

    def to_dict(self) -> dict[str, Any]:
        """Return JSON serializable dictionary representation."""
        circuit_status_dict: dict[str, Any] = {
            "status": self.status.name,
            "message": self.message,
        }
        if self.error_detail is not None:
            circuit_status_dict["error_detail"] = self.error_detail

        if self.completed_time is not None:
            circuit_status_dict["completed_time"] = self.completed_time.isoformat()
        if self.queued_time is not None:
            circuit_status_dict["queued_time"] = self.queued_time.isoformat()
        if self.submitted_time is not None:
            circuit_status_dict["submitted_time"] = self.submitted_time.isoformat()
        if self.running_time is not None:
            circuit_status_dict["running_time"] = self.running_time.isoformat()
        if self.cancelled_time is not None:
            circuit_status_dict["cancelled_time"] = self.cancelled_time.isoformat()
        if self.error_time is not None:
            circuit_status_dict["error_time"] = self.error_time.isoformat()

        if self.queue_position is not None:
            circuit_status_dict["queue_position"] = self.queue_position

        return circuit_status_dict

    @classmethod
    def from_dict(cls, dic: dict[str, Any]) -> "CircuitStatus":
        """Construct from JSON serializable dictionary."""
        invalid = ValueError(f"Dictionary invalid format for CircuitStatus: {dic}")
        if "message" not in dic or "status" not in dic:
            raise invalid

        try:
            status = next(s for s in StatusEnum if dic["status"] == s.name)
        except StopIteration as e:
            raise invalid from e

        error_detail = dic.get("error_detail")

        read_optional_datetime: Callable[[str], datetime | None] = lambda key: (
            datetime.fromisoformat(x) if (x := dic.get(key)) is not None else None
        )
        completed_time = read_optional_datetime("completed_time")
        queued_time = read_optional_datetime("queued_time")
        submitted_time = read_optional_datetime("submitted_time")
        running_time = read_optional_datetime("running_time")
        cancelled_time = read_optional_datetime("cancelled_time")
        error_time = read_optional_datetime("error_time")

        queue_position = dic.get("queue_position")

        return cls(
            status,
            dic["message"],
            error_detail,
            completed_time,
            queued_time,
            submitted_time,
            running_time,
            cancelled_time,
            error_time,
            queue_position,
        )


WAITING_STATUS = {StatusEnum.QUEUED, StatusEnum.SUBMITTED, StatusEnum.RUNNING}
