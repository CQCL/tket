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

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, ClassVar, Dict, Optional, Type, TypeVar
from uuid import UUID
from dataclasses import asdict
import json
import os


def get_config_file_path() -> Path:
    """Get a path to the config file on this machine."""
    config_dir: Path
    xdg_conifg_dir = os.environ.get("XDG_CONFIG_HOME")
    if xdg_conifg_dir is None:
        config_dir = Path.home() / ".config"
    else:
        config_dir = Path(xdg_conifg_dir)

    pytket_config_file = config_dir / "pytket" / "config.json"

    return pytket_config_file


class PytketConfig:
    """PytketConfig represents a loaded config file for
    pytket and extension packages."""

    enable_telemetry: bool
    telemetry_id: Optional[UUID]
    extensions: Dict[str, Any]

    def __init__(
        self,
        enable_telemetry: bool,
        telemetry_id: Optional[UUID],
        extensions: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Construct a PytketConfig object with inital config parameter values.

        :param enable_telemetry: Set pytket telemetery on.
        :type enable_telemetry: bool
        :param telemetry_id: UUID identifying this system for telemetery
        :type telemetry_id: Optional[UUID]
        :param extensions: Dictionary holding parameter values for extension packages,
            defaults to None
        :type extensions: Optional[Dict[str, Any]], optional
        """

        self.enable_telemetry = enable_telemetry
        self.telemetry_id = telemetry_id
        self.extensions = {} if extensions is None else extensions

    @classmethod
    def default(cls) -> "PytketConfig":
        """Construct a default PytketConfig"""
        return PytketConfig(enable_telemetry=False, telemetry_id=None)

    @classmethod
    def read_file(cls, config_file_path: Path) -> "PytketConfig":
        """Construct a PytketConfig from reading a file with a given Path."""
        with config_file_path.open("r", encoding="utf-8") as config_file:
            config = json.load(config_file)
            return PytketConfig(
                config.get("enable_telemetry", False),
                config.get("telemetry_id", None),
                config.get("extensions", dict()),
            )

    def write_file(self, config_file_path: Path) -> None:
        """Write a PytketConfig to a file with a given Path."""
        config_file_path.parent.mkdir(parents=True, exist_ok=True)
        with config_file_path.open("w", encoding="utf-8") as config_file:
            config = {
                "enable_telemetry": self.enable_telemetry,
                "telemetry_id": self.telemetry_id,
                "extensions": self.extensions,
            }
            json.dump(config, config_file, indent=2)


def load_config_file() -> PytketConfig:
    """Load config from default file path."""
    return PytketConfig.read_file(get_config_file_path())


def write_config_file(config: PytketConfig) -> None:
    """Write config to default file path."""
    config.write_file(get_config_file_path())


T_ext = TypeVar("T_ext", bound="PytketExtConfig")


class PytketExtConfig(ABC):
    """Abstract base class for pytket extension config classes."""

    ext_dict_key: ClassVar[str] = ""

    @classmethod
    @abstractmethod
    def from_extension_dict(cls: Type[T_ext], ext_dict: Dict[str, Any]) -> T_ext:
        """Abstract method to build PytketExtConfig from dictionary serialized form."""
        ...

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return asdict(self)

    @classmethod
    def from_pytketconfig(cls: Type[T_ext], p_config: PytketConfig) -> T_ext:
        """Build from PytketConfig instance."""
        if cls.ext_dict_key in p_config.extensions:
            return cls.from_extension_dict(p_config.extensions[cls.ext_dict_key])
        return cls.from_extension_dict({})

    @classmethod
    def from_default_config_file(cls: Type[T_ext]) -> T_ext:
        """Load from default config file."""
        return cls.from_pytketconfig(load_config_file())

    def update_pytket_config(self, pytket_config: PytketConfig) -> None:
        """Update a PytketConfig instance from this extension config."""
        pytket_config.extensions.update({self.ext_dict_key: self.to_dict()})

    def update_default_config_file(self) -> None:
        """Update default config file with current parameters
        in this extension config."""
        config = load_config_file()
        self.update_pytket_config(config)
        write_config_file(config)
