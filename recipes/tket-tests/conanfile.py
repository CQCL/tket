# Copyright 2019-2022 Cambridge Quantum Computing
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

from conans import ConanFile, CMake, tools
from conans.errors import ConanInvalidConfiguration
import os, platform


class TketTestsConan(ConanFile):
    name = "tket-tests"
    version = "0.6.2"
    license = "CQC Proprietary"
    author = "Alec Edgington <alec.edgington@cambridgequantum.com>"
    url = "https://github.com/CQCL/tket"
    description = "Unit tests for tket"
    topics = ("quantum", "computation", "compiler")
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "with_coverage": [True, False],
        "full": [True, False],
        "long": [True, False],
    }
    default_options = {"with_coverage": False, "full": False, "long": False}
    generators = "cmake"
    exports_sources = "../../tket/tests/*"
    requires = ("tket/1.0.1", "catch2/3.0.1")

    _cmake = None

    def validate(self):
        if self.options.with_coverage and self.settings.compiler != "gcc":
            raise ConanInvalidConfiguration(
                "`with_coverage` option only available with gcc"
            )

    def _configure_cmake(self):
        if self._cmake is None:
            self._cmake = CMake(self)
            self._cmake.configure()
        return self._cmake

    def configure(self):
        if self.options.with_coverage:
            self.options["tket"].profile_coverage = True

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        executable_filename = "bin/test_tket"
        if platform.system() == "Windows":
            executable_filename = executable_filename + ".exe"
        self.copy(executable_filename)
        self.copy(
            os.path.join("test_circuits", "lexiroute_circuit.json"),
            dst="bin",
            keep_path=False,
        )
