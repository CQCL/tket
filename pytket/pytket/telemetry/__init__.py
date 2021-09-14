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

"""The telemetry module defines how pytket installations can register themselves."""
from pathlib import Path
from urllib import request
from urllib.error import HTTPError, URLError
from urllib.request import Request
import json
from logging import getLogger
import os
import platform
import sys

from pytket.config import PytketConfig, get_config_file_path

logger = getLogger(__name__)


def _register(pytket_config_file: Path) -> None:
    pytket_version: str
    try:
        # For python 3.8 onwards
        from importlib.metadata import version  # type: ignore

        pytket_version = version("pytket")
    except ImportError:
        import pkg_resources

        pytket_version = pkg_resources.get_distribution("pytket").version

    data = {
        "version": pytket_version,
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "system": platform.system(),
        "system_version": platform.version(),
        "system_release": platform.release(),
        "c_api_version": sys.api_version,
    }

    headers = {"Content-Type": "application/json"}
    json_data = json.dumps(data).encode("utf8")

    try:
        resp = request.urlopen(
            Request(
                "https://telemetry.cambridgequantum.com/v3/register",
                json_data,
                headers,
            ),
            timeout=10,
        )
        if resp.status != 200:
            logger.error(
                "failed to register pytket with http status code: %s", resp.status
            )
        else:
            resp_body = json.loads(resp.read(64).decode("utf-8"))
            telemetry_id = resp_body["telemetry_id"]

            telemetry_config = PytketConfig.read_file(pytket_config_file)
            telemetry_config.telemetry_id = telemetry_id
            telemetry_config.write_file(pytket_config_file)

    except (URLError, HTTPError) as err:
        logger.error("failed to register pytket with exception: %s", err)


def _set_telemetry_preference(enabled: bool) -> None:
    pytket_config_file = get_config_file_path()
    if not pytket_config_file.exists():
        default_config = PytketConfig.default()
        default_config.write_file(pytket_config_file)

    telemetry_config = PytketConfig.read_file(pytket_config_file)
    telemetry_config.enable_telemetry = enabled
    telemetry_config.write_file(pytket_config_file)


def opt_in() -> None:
    """Opt into pytket telemetry"""
    _set_telemetry_preference(True)
    print("Successfully opted into telemetry")


def opt_out() -> None:
    """Opt out of pytket telemetry"""
    _set_telemetry_preference(False)
    print("Successfully opted out of telemetry")


def _on_module_load() -> None:
    config: PytketConfig
    pytket_config_file = get_config_file_path()
    if pytket_config_file.exists():
        config = PytketConfig.read_file(pytket_config_file)
    else:
        config = PytketConfig.default()
        config.write_file(pytket_config_file)

    if config.enable_telemetry and config.telemetry_id is None:
        _register(pytket_config_file)


_on_module_load()
