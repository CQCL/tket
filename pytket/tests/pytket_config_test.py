# Copyright Quantinuum
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

import json
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar

import pytest
from jsonschema import validate  # type: ignore

from pytket.config import (
    PytketConfig,
    PytketExtConfig,
    get_config_file_path,
    load_config_file,
    write_config_file,
)


@dataclass
class SampleExtConfig(PytketExtConfig):
    ext_dict_key: ClassVar[str] = "tests_sample"

    field1: str | None
    field2: int | None

    @classmethod
    def from_extension_dict(
        cls: type["SampleExtConfig"], ext_dict: dict[str, Any]
    ) -> "SampleExtConfig":
        return cls(ext_dict.get("field1"), ext_dict.get("field2"))


def test_pytket_config() -> None:
    config_file = get_config_file_path()
    assert config_file.exists()
    with open(config_file) as f:
        config_dict = json.load(f)

    curr_file_path = Path(__file__).resolve().parent

    with open(curr_file_path.parent.parent / "schemas/pytket_config_v1.json") as f:
        schema = json.load(f)

    validate(instance=config_dict, schema=schema)

    config = load_config_file()

    assert isinstance(config.extensions, dict)


@pytest.fixture()
def pytket_config() -> Iterator[PytketConfig]:
    config = load_config_file()
    yield config
    if "tests_sample" in config.extensions:
        del config.extensions["tests_sample"]
    write_config_file(config)


def test_sample_extconfig(pytket_config: PytketConfig) -> None:
    assert "tests_sample" not in pytket_config.extensions
    pytket_config.extensions["tests_sample"] = {"field1": "foo"}

    write_config_file(pytket_config)

    assert "tests_sample" in load_config_file().extensions

    ext_config = SampleExtConfig.from_default_config_file()
    assert ext_config.field1 == "foo"
    assert ext_config.field2 is None

    ext_config.field1 = "bar"
    ext_config.field2 = 2

    ext_config.update_default_config_file()

    new_ext_config = SampleExtConfig.from_default_config_file()
    assert new_ext_config.field1 == "bar" and new_ext_config.field2 == 2
